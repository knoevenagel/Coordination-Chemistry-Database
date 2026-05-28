#!/usr/bin/env python3
"""
ChemDB Step10 - DID layer merger

Read five input CSVs, map records to DID across layers L1-L5, and write a single CSV:
Columns: DID, L1, L2, L3, L4, L5
- Each Lx contains a JSON string of all fields from the first matching row in that file
  (with some field renames), or <Empty> if absent.

Input locations (same as step9.py conventions):
- ./tmp/complex_data.csv            (L1, key: did)
- ./tmp/ligand_data.csv             (L2, key: ligand_did)
- ./tmp/repaired_ligand_data.csv    (L3, key: ligand_new_did)
- ./tmp/fragments.csv               (L4, key: fragment_DID)
- ./tmp/IRL_filtered.csv            (L5, key: DID)

Field renames per layer (source -> target):
- L1: complex_smiles -> smiles
- L2: ligand_smiles  -> smiles
- L3: new_smiles     -> smiles, is_repaired (kept as-is)
- L4: fragment_smiles-> smiles
- L5: SMILES         -> smiles, GAC -> gac

Notes:
- Case-insensitive matching for column names when reading keys and rename sources.
- On duplicate keys within a single file, keep only the first occurrence.
- Uses chunked reading to limit memory usage.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, Optional, Set, Tuple, List
import time

import pandas as pd

# Default configuration variables (aligned with step9.py)
INPUT_DIR = "./tmp"
DATA_DIR = "./data"
OUTPUT_DIR = "./tmp"

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def _find_column_case_insensitive(columns: List[str], target: str) -> Optional[str]:
    """Return the actual column name matching target case-insensitively, or None."""
    tl = target.lower()
    for c in columns:
        if c.lower() == tl:
            return c
    return None


def _apply_renames(row_dict: Dict[str, object], renames: Dict[str, str]) -> Dict[str, object]:
    """Apply case-insensitive renames in-place: for each src->dst, move if present.
    - If src present (case-insensitive), set dst and remove the original key.
    - Do not overwrite an existing dst if already present with a non-null value.
    - Remove pandas NaN by converting to None.
    """
    # Normalize keys map: lower -> actual
    key_map = {k.lower(): k for k in list(row_dict.keys())}
    for src, dst in renames.items():
        src_l = src.lower()
        if src_l in key_map:
            actual_src = key_map[src_l]
            val = row_dict.get(actual_src)
            # Normalize NaN/NaT to None
            if isinstance(val, float) and pd.isna(val):
                val = None
            if val is None or (isinstance(val, str) and val.strip().lower() == 'nan'):
                # treat as missing
                del row_dict[actual_src]
                continue
            # Only set dst if not already present or present but null-like
            if dst not in row_dict or row_dict[dst] in (None, '', 'nan'):
                row_dict[dst] = val
            # Remove the original key if different
            if actual_src != dst and actual_src in row_dict:
                del row_dict[actual_src]
    # Normalize remaining NaN to None
    for k in list(row_dict.keys()):
        v = row_dict[k]
        if isinstance(v, float) and pd.isna(v):
            row_dict[k] = None
        elif isinstance(v, str) and v.strip().lower() == 'nan':
            row_dict[k] = None
    return row_dict


def _build_selected_record(cols: List[str], row_tuple: Tuple[object, ...], select_map: Dict[str, str]) -> Dict[str, object]:
    """Select only requested fields from a row (tuple) and apply renames.
    - Case-insensitive source matching based on provided cols list.
    - Skips fields that are missing or NaN/None.
    """
    result: Dict[str, object] = {}
    # Build lookup map: lowercased column -> (actual_name, index)
    col_index = {c.lower(): (c, i) for i, c in enumerate(cols)}
    for src, dst in select_map.items():
        key = src.lower()
        if key not in col_index:
            continue
        _, idx = col_index[key]
        val = row_tuple[idx]
        if val is None or (isinstance(val, float) and pd.isna(val)):
            continue
        if isinstance(val, str) and val.strip().lower() == 'nan':
            continue
        result[dst] = val
    return result


class Step10Processor:
    def __init__(self, record_limit: int = -1):
        self.input_dir = Path(INPUT_DIR)
        self.data_dir = Path(DATA_DIR)
        self.output_dir = Path(OUTPUT_DIR)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Output file
        self.output_file = self.output_dir / "step10_layers.csv"

                # Specs for each layer: (layer_name, file_path, key_column, select_map)
        # select_map: desired source column -> output name (case-insensitive source matching)
        self.layer_specs = [
            (
                "L1",
                self.input_dir / "complex_data.csv",
                "did",
                {"cid": "cid", "complex_smiles": "smiles", "ligand_count": "ligand_count"},
            ),
            (
                "L2",
                self.input_dir / "ligand_data.csv",
                "ligand_did",
                {"ligand_smiles": "smiles"},
            ),
            (
                "L3",
                self.input_dir / "repaired_ligand_data.csv",
                "ligand_new_did",
                {"new_smiles": "smiles", "is_repaired": "is_repaired"},
            ),
            (
                "L4",
                self.input_dir / "fragments.csv",
                "fragment_DID",
                {"fragment_smiles": "smiles"},
            ),
            (
                "L5",
                self.data_dir / "IRL_filtered.csv",
                "DID",
                {"SMILES": "smiles", "GAC": "gac"},
            ),
        ]
        self.record_limit = record_limit  # limit of total unique DIDs to include; -1 for no limit
        self.elapsed_seconds: float = 0.0

    def _read_first_hit_map(
        self, csv_path: Path, key_col_name: str, select_map: Dict[str, str]
    ) -> Tuple[Dict[str, Dict[str, object]], Set[str]]:
        """Stream through CSV and build DID -> row_dict (first occurrence only).
        Returns (mapping, dids_seen)
        """
        mapping: Dict[str, Dict[str, object]] = {}
        dids_seen: Set[str] = set()

        if not csv_path.exists():
            logger.warning(f"File not found, skip: {csv_path}")
            return mapping, dids_seen

        logger.info(f"Loading: {csv_path.name}")
        chunksize = 10000
        try:
            for chunk in pd.read_csv(csv_path, dtype=str, chunksize=chunksize):
                cols = list(chunk.columns)
                actual_key_col = _find_column_case_insensitive(cols, key_col_name)
                if not actual_key_col:
                    logger.warning(f"Key column '{key_col_name}' not found in {csv_path.name}, skipping chunk")
                    continue

                # Pre-compute column index for fast tuple access
                col_index = {c: i for i, c in enumerate(cols)}
                key_idx = col_index[actual_key_col]

                for row in chunk.itertuples(index=False, name=None):
                    key_val = row[key_idx]
                    if key_val is None or (isinstance(key_val, float) and pd.isna(key_val)):
                        continue
                    did = str(key_val).strip()
                    if not did or did.lower() == 'nan':
                        continue
                    if did in mapping:
                        continue  # keep first occurrence

                    # Build filtered/renamed record
                    record = _build_selected_record(cols, row, select_map)
                    if not record:
                        continue  # treat as absent if none of the selected fields are present
                    mapping[did] = record
                    dids_seen.add(did)

                    # Optional global limit based on union size; check externally
                if self.record_limit > 0 and len(dids_seen) >= self.record_limit:
                    logger.info("Record limit reached while reading file; stopping early")
                    break
        except Exception as e:
            logger.error(f"Failed reading {csv_path}: {e}")

        logger.info(f"Collected {len(mapping)} first-hit records from {csv_path.name}")
        return mapping, dids_seen

    def process(self) -> Path:
        # Load each layer map
        layer_maps: Dict[str, Dict[str, Dict[str, object]]] = {}
        all_dids: Set[str] = set()
        start_time = time.time()

        for layer_name, file_path, key_col, select_map in self.layer_specs:
            m, dids = self._read_first_hit_map(file_path, key_col, select_map)
            layer_maps[layer_name] = m
            all_dids.update(dids)
            if self.record_limit > 0 and len(all_dids) >= self.record_limit:
                logger.info(f"Global record limit reached at layer {layer_name}")
                break

        if not all_dids:
            logger.warning("No DIDs found across all inputs. Writing an empty CSV header only.")
            pd.DataFrame(columns=["DID", "L1", "L2", "L3", "L4", "L5"]).to_csv(
                self.output_file, index=False
            )
            return self.output_file

        # Build output rows
        rows = []
        for did in sorted(all_dids):
            out_row = {"DID": did}
            for layer_name, _, _, _ in self.layer_specs:
                record = layer_maps.get(layer_name, {}).get(did)
                if record is None:
                    out_row[layer_name] = None
                else:
                    try:
                        out_row[layer_name] = json.dumps(record, ensure_ascii=False, separators=(",", ":"))
                    except Exception:
                        # Fallback: convert non-serializable values to string
                        safe_record = {k: (None if v is None else str(v)) for k, v in record.items()}
                        out_row[layer_name] = json.dumps(safe_record, ensure_ascii=False, separators=(",", ":"))
            rows.append(out_row)

        df = pd.DataFrame(rows, columns=["DID", "L1", "L2", "L3", "L4", "L5"])
        # Write CSV with explicit <Empty> for missing values
        df.to_csv(self.output_file, index=False, na_rep='')
        self.elapsed_seconds = time.time() - start_time
        logger.info(f"Output written to: {self.output_file} (elapsed {self.elapsed_seconds:.2f}s)")
        return self.output_file


def main():
    import sys
    record_limit = -1
    if len(sys.argv) > 1:
        try:
            record_limit = int(sys.argv[1])
            print(f"Using record limit: {record_limit}")
        except ValueError:
            print("Invalid record limit, using default (-1)")
    processor = Step10Processor(record_limit=record_limit)
    out_path = processor.process()
    if getattr(processor, 'elapsed_seconds', None) is not None:
        print(f"Step10 completed. Output: {out_path} | Elapsed: {processor.elapsed_seconds:.2f}s")
    else:
        print(f"Step10 completed. Output: {out_path}")


if __name__ == "__main__":
    main()
