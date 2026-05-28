#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare TSV data under data/mysql/ with generated CSVs under tmp/ to find core row differences.

Outputs per-pair reports into log/:
- comparison_<name>_missing_in_right.csv (keys present in left but missing in right)
- comparison_<name>_missing_in_left.csv (keys present in right but missing in left)
- comparison_<name>_mismatches.csv (common keys with column value differences)
Also prints concise on-screen summary and first few rows (heads) from each file.

Usage:
  python src/comparison.py
"""
from __future__ import annotations

import ast
import csv
import json
import os
import sys
from dataclasses import dataclass
import argparse
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
DATA_DIR = os.path.join(PROJECT_ROOT, "data", "mysql")
TMP_DIR = os.path.join(PROJECT_ROOT, "tmp")
LOG_DIR = os.path.join(PROJECT_ROOT, "tmp", "comparison")
os.makedirs(LOG_DIR, exist_ok=True)


def _read_table(path: str, sep: Optional[str] = None, nrows: Optional[int] = None) -> pd.DataFrame:
    """Robust reader for TSV/CSV with encoding fallbacks.
    - sep: autodetect by file extension if not provided
    - nrows: limit rows for quick head preview
    """
    if sep is None:
        lower = path.lower()
        if lower.endswith(".tsv"):
            sep = "\t"
        elif lower.endswith(".csv"):
            sep = ","
        else:
            # best effort: try tab first then comma
            sep = None
    encodings = ["utf-8", "utf-8-sig", "latin-1"]
    last_err: Optional[Exception] = None
    for enc in encodings:
        try:
            return pd.read_csv(
                path,
                sep=sep,
                engine="python",  # tolerant parser
                dtype="object",   # keep as strings for stable comparison
                na_filter=False,   # keep literal strings like "NULL" and "[]"
                quoting=csv.QUOTE_MINIMAL,
                quotechar='"',
                escapechar='\\',
                on_bad_lines="skip",
                nrows=nrows,
                encoding=enc,
            )
        except Exception as e:  # try next encoding
            last_err = e
            continue
    raise RuntimeError(f"Failed reading {path!r}: {last_err}")


def _to_bool_like(x: str) -> Optional[bool]:
    if x is None:
        return None
    s = str(x).strip().lower()
    if s in {"1", "true", "t", "yes", "y"}:
        return True
    if s in {"0", "false", "f", "no", "n"}:
        return False
    # allow numeric 0/1
    try:
        if float(s) == 1.0:
            return True
        if float(s) == 0.0:
            return False
    except Exception:
        pass
    return None  # not a boolean


def _parse_list_like(x: str):
    """Attempt to parse list/dict-like strings (e.g., "[]", "[{'Ca': '0'}]")."""
    if x is None:
        return []
    s = str(x).strip()
    if s == "" or s.upper() == "NULL":
        return []
    # normalize quotes for ast
    try:
        return ast.literal_eval(s)
    except Exception:
        return s  # fallback to raw string


def _normalize_series(name: str, s: pd.Series) -> pd.Series:
    lname = name.lower()
    if lname in {"from_pubchem", "inactive", "ligand_inactive"}:
        return s.apply(_to_bool_like).astype("object")
    if lname in {"label_list", "labels", "metal_info"}:
        return s.apply(_parse_list_like).astype("object")
    # default: strip whitespace
    return s.astype("string").str.strip().astype("object")


@dataclass
class PairSpec:
    name: str
    left_path: str
    right_path: str
    left_key: Union[str, Sequence[str]]
    right_key: Union[str, Sequence[str]]
    # Optional explicit column name mapping from left/right to unified names
    left_renames: Optional[Dict[str, str]] = None
    right_renames: Optional[Dict[str, str]] = None
    # If provided, only these unified columns will be compared; otherwise uses intersection
    compare_columns: Optional[Sequence[str]] = None
    # Columns to ignore in auto-intersection
    ignore_columns: Optional[Sequence[str]] = None


def _prep_df(path: str, key_col: Union[str, Sequence[str]], renames: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    df = _read_table(path)
    # lower-case columns for robust matching
    df.columns = [c.strip() for c in df.columns]
    lower_map = {c: c.lower() for c in df.columns}
    df.rename(columns=lower_map, inplace=True)
    # apply additional renames (using lowered keys)
    ren_low: Dict[str, str] = {}
    if renames:
        # ensure both mapping keys and values are lowered for consistency
        ren_low = {k.lower(): v.lower() for k, v in renames.items()}
        df.rename(columns=ren_low, inplace=True)
    # resolve effective key name(s) (after renames)
    if isinstance(key_col, (list, tuple)):
        key_names = list(key_col)
    else:
        key_names = [key_col]

    def _resolve_one(kname: str) -> str:
        k_low = str(kname).lower()
        effective = ren_low.get(k_low, k_low)
        if effective not in df.columns:
            raise KeyError(
                f"Key column {kname!r} (effective {effective!r}) not found in {path}. Available: {list(df.columns)}"
            )
        return effective

    resolved_keys = [_resolve_one(k) for k in key_names]

    # build key column: string for single key; tuple for composite key
    if len(resolved_keys) == 1:
        k = resolved_keys[0]
        df["__key__"] = df[k].astype("string").str.strip()
    else:
        parts = [df[k].astype("string").str.strip() for k in resolved_keys]
        df["__key__"] = list(zip(*parts))
    # normalize values for all columns
    for c in df.columns:
        df[c] = _normalize_series(c, df[c])
    # drop exact duplicate rows per key keeping the first
    if df["__key__"].duplicated().any():
        df = df.drop_duplicates(subset=["__key__"], keep="first")
    return df


def compare_pair(spec: PairSpec, print_heads: bool = True, sample_mismatches: int = 0) -> Dict[str, int]:
    left = _prep_df(spec.left_path, spec.left_key, spec.left_renames)
    right = _prep_df(spec.right_path, spec.right_key, spec.right_renames)

    # print heads
    if print_heads:
        print(f"\n=== {spec.name}: LEFT head ({os.path.basename(spec.left_path)}) ===")
        print(left.head(5).to_string(index=False))
        print(f"\n=== {spec.name}: RIGHT head ({os.path.basename(spec.right_path)}) ===")
        print(right.head(5).to_string(index=False))

    # decide columns to compare (unified names after renames/lower)
    if spec.compare_columns is not None:
        cols = [c.lower() for c in spec.compare_columns]
    else:
        left_cols = set(left.columns)
        right_cols = set(right.columns)
        ignore = {"__key__"}
        if spec.ignore_columns:
            ignore |= {c.lower() for c in spec.ignore_columns}
        cols = sorted((left_cols & right_cols) - ignore)
    # ensure key included only for merges
    cols = [c for c in cols if c != "__key__"]

    # Build merged view
    lsub = left[["__key__"] + cols].copy()
    rsub = right[["__key__"] + cols].copy()

    merged = lsub.merge(rsub, on="__key__", how="outer", indicator=True, suffixes=("_left", "_right"))

    # missing on either side
    missing_in_right = merged[merged["_merge"] == "left_only"]["__key__"]
    missing_in_left = merged[merged["_merge"] == "right_only"]["__key__"]

    # mismatches on common keys
    common = merged[merged["_merge"] == "both"].copy()
    mismatch_rows = []
    per_col_mismatch = {c: 0 for c in cols}
    for _, row in common.iterrows():
        diffs = {}
        for c in cols:
            lv = row.get(f"{c}_left")
            rv = row.get(f"{c}_right")
            # normalize complex types to comparable strings
            def norm_val(v):
                if isinstance(v, (list, tuple)):
                    # order-insensitive list comparison by stringifying and sorting
                    try:
                        return json.dumps(sorted([str(x) for x in v]))
                    except Exception:
                        return json.dumps(sorted([str(x) for x in list(v)]))
                if isinstance(v, dict):
                    return json.dumps(v, sort_keys=True)
                return v
            if norm_val(lv) != norm_val(rv):
                diffs[c + "__left"] = lv
                diffs[c + "__right"] = rv
                per_col_mismatch[c] += 1
        if diffs:
            mismatch_rows.append({"__key__": row["__key__"], **diffs})

    # Write reports
    def write_csv(fname: str, series_or_df):
        path = os.path.join(LOG_DIR, fname)
        if isinstance(series_or_df, pd.Series):
            series_or_df.to_frame(name="key").to_csv(path, index=False)
        else:
            pd.DataFrame(series_or_df).to_csv(path, index=False)
        return path

    m_left_path = write_csv(f"comparison_{spec.name}_missing_in_left.csv", missing_in_left)
    m_right_path = write_csv(f"comparison_{spec.name}_missing_in_right.csv", missing_in_right)
    mism_path = write_csv(f"comparison_{spec.name}_mismatches.csv", mismatch_rows)
    if sample_mismatches and mismatch_rows:
        sample = mismatch_rows[: sample_mismatches]
        write_csv(f"comparison_{spec.name}_mismatches_sample.csv", sample)

    summary = {
        "left_rows": int(left.shape[0]),
        "right_rows": int(right.shape[0]),
        "missing_in_left": int(missing_in_left.shape[0]),
        "missing_in_right": int(missing_in_right.shape[0]),
        "mismatches": int(len(mismatch_rows)),
    }

    print(
        f"\n[{spec.name}] rows L={summary['left_rows']} R={summary['right_rows']} | "
        f"missing(L)={summary['missing_in_left']} missing(R)={summary['missing_in_right']} "
        f"mismatches={summary['mismatches']}\n"
        f"  -> {os.path.basename(m_left_path)}, {os.path.basename(m_right_path)}, {os.path.basename(mism_path)}"
    )
    # print top per-column mismatch counts
    if summary['mismatches']:
        top_cols = sorted(per_col_mismatch.items(), key=lambda x: x[1], reverse=True)[:5]
        print("  top mismatch columns:")
        for c, n in top_cols:
            if n:
                print(f"    - {c}: {n}")

    return summary


def build_pairs() -> List[PairSpec]:
    return [
        # Complexes: DID is key; compare shared semantic columns
        PairSpec(
            name="complex",
            left_path=os.path.join(DATA_DIR, "ComplexRawData250704.tsv"),
            right_path=os.path.join(PROJECT_ROOT, "tmp", "complex_data.csv"),
            left_key="DID",
            right_key="did",
            left_renames={"DID": "did", "FromPubChem": "from_pubchem"},
            right_renames={},
            compare_columns=["complex_smiles", "metal_info", "from_pubchem", "inactive"],
        ),
        # Ligands: TSV DID corresponds to CSV ligand_did (not csv.did). Compare core fields
        PairSpec(
            name="ligand",
            left_path=os.path.join(DATA_DIR, "LigandRawData250704.tsv"),
            right_path=os.path.join(PROJECT_ROOT, "tmp", "ligand_data.csv"),
            left_key="DID",
            right_key="ligand_did",
            left_renames=None,
            right_renames=None,
            compare_columns=["ligand_smiles"],
        ),
        # Repaired ligands: schema aligns closely
        PairSpec(
            name="repaired_ligand",
            left_path=os.path.join(DATA_DIR, "RepairedLigandFull0704.tsv"),
            right_path=os.path.join(PROJECT_ROOT, "tmp", "repaired_ligand_data.csv"),
            left_key="ligand_new_did",
            right_key="ligand_new_did",
            compare_columns=["source_did", "new_smiles", "old_smiles", "is_repaired"],
        ),
        # PubChem: prefer intersection; ignore auxiliary columns
        PairSpec(
            name="pubchem",
            left_path=os.path.join(DATA_DIR, "PubChemRawData250704.tsv"),
            right_path=os.path.join(PROJECT_ROOT, "tmp", "pubchem_data.csv"),
            left_key="cid",
            right_key="cid",
            ignore_columns=["source_file", "row_index"],
        ),
        # Fragments: unique by (fragment_DID, source_DID); compare core columns
        PairSpec(
            name="fragments",
            left_path=os.path.join(DATA_DIR, "FragmentFull0704.tsv"),
            right_path=os.path.join(PROJECT_ROOT, "tmp", "fragments.csv"),
            left_key=["fragment_DID", "source_DID"],
            right_key=["fragment_DID", "source_DID"],
            left_renames={"fragment_DID": "fragment_did", "source_DID": "source_did"},
            right_renames={"fragment_DID": "fragment_did"},
            compare_columns=["source_did", "fragment_smiles", "fragment_irl_did", "fragment_irl_smiles"],
        ),
    ]

def main(argv: Optional[Sequence[str]] = None):
    parser = argparse.ArgumentParser(description="Compare TSVs under data/mysql with CSVs under tmp")
    parser.add_argument("--pairs", type=str, default="",
                        help="Comma-separated pair names to compare (default: all). Choices: complex,ligand,repaired_ligand,pubchem,fragments")
    parser.add_argument("--no-heads", action="store_true", help="Do not print heads for each pair")
    parser.add_argument("--sample-mismatches", type=int, default=0, help="Also write a small sample of mismatches per pair")
    args = parser.parse_args(argv)

    all_pairs = build_pairs()
    if args.pairs:
        wanted = {p.strip() for p in args.pairs.split(",") if p.strip()}
        pairs = [p for p in all_pairs if p.name in wanted]
    else:
        pairs = all_pairs

    all_summaries: List[Tuple[str, Dict[str, int]]] = []
    for spec in pairs:
        try:
            summary = compare_pair(spec, print_heads=not args.no_heads, sample_mismatches=args.sample_mismatches)
            all_summaries.append((spec.name, summary))
        except Exception as e:
            print(f"\n!! Error comparing {spec.name}: {e}")

    # write overall summary
    summary_lines = [
        "pair,left_rows,right_rows,missing_in_left,missing_in_right,mismatches"
    ]
    for name, s in all_summaries:
        summary_lines.append(
            f"{name},{s['left_rows']},{s['right_rows']},{s['missing_in_left']},{s['missing_in_right']},{s['mismatches']}"
        )
    with open(os.path.join(LOG_DIR, "comparison_summary.csv"), "w", encoding="utf-8") as f:
        f.write("\n".join(summary_lines))
    print(f"\nSummary written to {LOG_DIR}/comparison_summary.csv")


if __name__ == "__main__":
    main()
