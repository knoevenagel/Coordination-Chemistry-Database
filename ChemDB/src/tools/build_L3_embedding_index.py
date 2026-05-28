#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build L3 embedding index files from tmp L3 list and DID→SMILES mapping.
Uses MolCLR (GIN, GCN) and ECFP fingerprints; saves DID–smiles–embedding to data/L3_embedding.

L3 list: tmp/metal_l3_index.csv (unique l3_did).
DID→SMILES: tmp/repaired_ligand_data.csv (ligand_new_did, new_smiles).
Output: data/L3_embedding/L3_embedding_{gin,gcn,ecfp}.npz with keys: dids, smiles, embeddings.
GCN requires torch_sparse / torch_scatter; if missing, GCN index is skipped.
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import sys
import time
from typing import Dict, List, Tuple

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(SRC_DIR)
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)


def setup_logging(log_file: str | None) -> None:
    """Log to file and stdout; flush after each message."""
    handlers: List[logging.Handler] = [logging.StreamHandler(sys.stdout)]
    if log_file:
        os.makedirs(os.path.dirname(log_file) or ".", exist_ok=True)
        handlers.append(logging.FileHandler(log_file, encoding="utf-8"))
    for h in handlers:
        h.setLevel(logging.INFO)
        h.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logging.basicConfig(level=logging.INFO, handlers=handlers, force=True)
    # Ensure stdout is line-buffered so logs show up
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(line_buffering=True)

# Default paths under repo
DEFAULT_TMP_DIR = os.path.join(REPO_ROOT, "tmp")
DEFAULT_OUT_DIR = os.path.join(REPO_ROOT, "data", "L3_embedding")
L3_INDEX_NAME = "metal_l3_index.csv"
DID_SMILES_NAME = "repaired_ligand_data.csv"


def load_l3_did_smiles(tmp_dir: str) -> List[Tuple[str, str]]:
    """Load full L3 list (DID → SMILES) from tmp: metal_l3_index.csv + repaired_ligand_data.csv."""
    l3_path = os.path.join(tmp_dir, L3_INDEX_NAME)
    did_smiles_path = os.path.join(tmp_dir, DID_SMILES_NAME)
    if not os.path.isfile(l3_path):
        raise FileNotFoundError(f"L3 index not found: {l3_path}")
    if not os.path.isfile(did_smiles_path):
        raise FileNotFoundError(f"DID–SMILES CSV not found: {did_smiles_path}")

    # Unique L3 DIDs (keep order)
    l3_dids: List[str] = []
    seen = set()
    with open(l3_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        lid_col = next((c for c in (reader.fieldnames or []) if "l3" in c.lower() and "did" in c.lower()), None)
        if not lid_col:
            lid_col = "l3_did" if "l3_did" in (reader.fieldnames or []) else None
        if not lid_col:
            raise ValueError(f"No L3 DID column in {l3_path}; columns: {reader.fieldnames}")
        for row in reader:
            did = (row.get(lid_col) or "").strip()
            if did and did not in seen:
                seen.add(did)
                l3_dids.append(did)

    # DID → SMILES from repaired_ligand_data (normalize: strip L3_ for lookup)
    did_to_smiles: Dict[str, str] = {}
    with open(did_smiles_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        did_col = next((c for c in (reader.fieldnames or []) if "did" in c.lower() and "new" in c.lower()), None)
        if not did_col:
            did_col = next((c for c in (reader.fieldnames or []) if c in ("ligand_new_did", "did")), "ligand_new_did")
        smiles_col = next((c for c in (reader.fieldnames or []) if "smiles" in c.lower() and "new" in c.lower()), None)
        if not smiles_col:
            smiles_col = next((c for c in (reader.fieldnames or []) if c in ("new_smiles", "smiles")), "new_smiles")
        for row in reader:
            did = (row.get(did_col) or "").strip()
            smiles = (row.get(smiles_col) or "").strip()
            if did and smiles:
                did_to_smiles[did] = smiles

    # Build (L3_DID, SMILES) for L3 that have a SMILES (normalize L3_DID to D... for lookup)
    out: List[Tuple[str, str]] = []
    for l3_did in l3_dids:
        key = l3_did[3:] if l3_did.upper().startswith("L3_") else l3_did  # L3_ -> strip 3 chars
        smiles = did_to_smiles.get(key)
        if smiles:
            out.append((l3_did, smiles))
    return out


def embed_molclr(smiles_list: List[str], model_type: str, device: str = "cpu") -> np.ndarray:
    """Return (N, D) embeddings for valid SMILES; invalid positions get NaN row (caller should drop)."""
    from molclr_api import extract_embeddings, get_model

    model = get_model(model_type, device=device)
    embeddings, valid_indices, _ = extract_embeddings(model, smiles_list, device=device)
    # Map back to full list: valid_indices[i] -> embeddings[i]
    n = len(smiles_list)
    dim = embeddings.shape[1] if len(embeddings) else 512
    out = np.full((n, dim), np.nan, dtype=np.float32)
    for idx, i in enumerate(valid_indices):
        out[i] = embeddings[idx]
    return out


def _log(msg: str) -> None:
    logging.info(msg)
    sys.stdout.flush()
    sys.stderr.flush()


def embed_ecfp(smiles_list: List[str], radius: int = 2, n_bits: int = 2048) -> np.ndarray:
    """Morgan/ECFP fingerprint; invalid SMILES -> zero vector."""
    from rdkit import Chem
    from rdkit.Chem import AllChem, DataStructs
    n = len(smiles_list)
    out = np.zeros((n, n_bits), dtype=np.float32)
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)
            DataStructs.ConvertToNumpyArray(fp, out[i])
    return out


def main():
    parser = argparse.ArgumentParser(description="Build L3 embedding indices (GIN, GCN, ECFP, RDKit FP).")
    parser.add_argument("--tmp-dir", default=DEFAULT_TMP_DIR, help="Tmp dir containing metal_l3_index.csv and repaired_ligand_data.csv")
    parser.add_argument("--out-dir", default=DEFAULT_OUT_DIR, help="Output directory for index files")
    parser.add_argument("--max-samples", type=int, default=None, help="Limit to N L3 (e.g. 100 for testing)")
    parser.add_argument("--estimate", type=int, nargs="?", const=1_000_000, metavar="N", help="Only estimate cost: run small batch and extrapolate to N (default 1e6)")
    parser.add_argument(
        "--backends",
        nargs="+",
        choices=["gin", "gcn", "ecfp"],
        default=["gin", "gcn", "ecfp"],
        help="Which backends to run (default: gin gcn ecfp)",
    )
    parser.add_argument("--device", default="cpu", help="Device for MolCLR")
    parser.add_argument("--radius", type=int, default=2, help="ECFP/RDKit FP radius")
    parser.add_argument("--n-bits", type=int, default=2048, help="ECFP/RDKit FP n_bits")
    parser.add_argument("--chunk-size", type=int, default=50_000, help="Process and write npz per N L3 (default 50k)")
    parser.add_argument("--log-file", default=None, help="Log file path (default: <out-dir>/build.log)")
    args = parser.parse_args()

    log_file = args.log_file
    if log_file is None and args.estimate is None:
        log_file = os.path.join(args.out_dir, "build.log")
    setup_logging(log_file)

    _log("Loading L3 DID → SMILES from tmp...")
    pairs = load_l3_did_smiles(args.tmp_dir)
    _log(f"Total L3 with SMILES: {len(pairs)}")

    if args.estimate is not None:
        target_n = args.estimate
        sample_n = min(200, len(pairs)) if len(pairs) >= 200 else max(50, len(pairs))
        pairs = pairs[:sample_n]
        _log(f"Estimate mode: timing {len(pairs)} L3, extrapolating to {target_n:,} L3")
    elif args.max_samples is not None:
        pairs = pairs[: args.max_samples]
        _log(f"Limited to first {len(pairs)} L3 (--max-samples={args.max_samples})")

    if not pairs:
        _log("No L3 to process. Exiting.")
        return

    dids = [p[0] for p in pairs]
    smiles_list = [p[1] for p in pairs]
    n_total = len(dids)
    os.makedirs(args.out_dir, exist_ok=True)
    chunk_size = max(1, args.chunk_size)
    n_chunks = (n_total + chunk_size - 1) // chunk_size

    def save_npz(dids_out: List[str], smiles_out: List[str], emb: np.ndarray, name: str) -> None:
        path = os.path.join(args.out_dir, name)
        np.savez(path, dids=np.array(dids_out, dtype=object), smiles=np.array(smiles_out, dtype=object), embeddings=emb)
        _log(f"  Saved {emb.shape[0]} rows to {path}")

    def run_backend_chunked(
        name: str,
        embed_fn,  # (List[str]) -> np.ndarray
        drop_invalid: bool = False,
    ) -> None:
        """Process in chunks: compute embeddings per chunk, write part npz immediately, then merge from disk."""
        part_dir = os.path.join(args.out_dir, f"_parts_{name}")
        os.makedirs(part_dir, exist_ok=True)
        part_paths: List[str] = []
        for c in range(n_chunks):
            start = c * chunk_size
            end = min(start + chunk_size, n_total)
            chunk_dids = dids[start:end]
            chunk_smiles = smiles_list[start:end]
            t0 = time.perf_counter()
            emb = embed_fn(chunk_smiles)
            elapsed = time.perf_counter() - t0
            if drop_invalid:
                valid = ~np.any(np.isnan(emb), axis=1)
                emb = emb[valid]
                chunk_dids = [x for x, v in zip(chunk_dids, valid) if v]
                chunk_smiles = [x for x, v in zip(chunk_smiles, valid) if v]
            part_path = os.path.join(part_dir, f"part_{c:04d}.npz")
            np.savez(
                part_path,
                dids=np.array(chunk_dids, dtype=object),
                smiles=np.array(chunk_smiles, dtype=object),
                embeddings=emb,
            )
            part_paths.append(part_path)
            _log(f"  {name} chunk {c+1}/{n_chunks} ({end:,}/{n_total:,} L3) done in {elapsed:.1f}s, wrote {part_path}")
        # Merge: load parts from disk and write final npz
        _log(f"  Merging {len(part_paths)} parts...")
        all_dids = []
        all_smiles = []
        all_emb = []
        for p in part_paths:
            with np.load(p, allow_pickle=True) as f:
                all_dids.extend(f["dids"].tolist())
                all_smiles.extend(f["smiles"].tolist())
                all_emb.append(f["embeddings"])
        final_emb = np.concatenate(all_emb, axis=0)
        all_emb.clear()
        final_path = os.path.join(args.out_dir, f"L3_embedding_{name}.npz")
        np.savez(
            final_path,
            dids=np.array(all_dids, dtype=object),
            smiles=np.array(all_smiles, dtype=object),
            embeddings=final_emb,
        )
        _log(f"  Merged {final_emb.shape[0]} rows -> {final_path}")
        for p in part_paths:
            try:
                os.remove(p)
            except OSError:
                pass
        try:
            os.rmdir(part_dir)
        except OSError:
            pass

    estimate_target = args.estimate
    timings: Dict[str, float] = {}

    if estimate_target is not None:
        # Estimate mode: small batch, no chunking, no part files
        # 1) GIN
        if "gin" in args.backends:
            _log("Computing MolCLR GIN embeddings...")
            try:
                t0 = time.perf_counter()
                embed_molclr(smiles_list, "gin", device=args.device)
                timings["gin"] = time.perf_counter() - t0
                _log(f"  GIN: {timings['gin']:.2f} s for {len(smiles_list)} L3")
            except Exception as e:
                _log(f"  Skipped GIN: {e}")
        # 2) GCN
        if "gcn" in args.backends:
            _log("Computing MolCLR GCN embeddings...")
            try:
                t0 = time.perf_counter()
                embed_molclr(smiles_list, "gcn", device=args.device)
                timings["gcn"] = time.perf_counter() - t0
                _log(f"  GCN: {timings['gcn']:.2f} s for {len(smiles_list)} L3")
            except Exception as e:
                _log(f"  Skipped GCN: {e}")
        # 3) ECFP
        if "ecfp" in args.backends:
            _log("Computing ECFP embeddings...")
            t0 = time.perf_counter()
            embed_ecfp(smiles_list, radius=args.radius, n_bits=args.n_bits)
            timings["ecfp"] = time.perf_counter() - t0
            _log(f"  ECFP: {timings['ecfp']:.2f} s for {len(smiles_list)} L3")
        n = len(smiles_list)
        scale = estimate_target / n
        _log("\n--- 预估 100 万 L3 开销 (Extrapolation to {:,} L3) ---".format(estimate_target))
        total_sec = 0.0
        for name in ("gin", "gcn", "ecfp"):
            if name not in timings:
                continue
            sec = timings[name] * scale
            total_sec += sec
            _log(f"  {name:12s}: {sec/3600:.1f} h ({sec:.0f} s)")
        _log(f"  合计 (串行): {total_sec/3600:.1f} h")
        _log("Done (estimate only; no files written).")
        return

    # Full run: chunked processing, write part npz per chunk then merge
    _log(f"Chunk size: {chunk_size:,} L3, {n_chunks} chunks per backend")

    # 1) MolCLR GIN
    if "gin" in args.backends:
        _log("Computing MolCLR GIN embeddings (chunked)...")
        try:
            t0 = time.perf_counter()
            run_backend_chunked(
                "gin",
                lambda sm: embed_molclr(sm, "gin", device=args.device),
                drop_invalid=True,
            )
            _log(f"  GIN total: {time.perf_counter() - t0:.1f} s")
        except Exception as e:
            _log(f"  Skipped GIN: {e}")

    # 2) MolCLR GCN
    if "gcn" in args.backends:
        _log("Computing MolCLR GCN embeddings (chunked)...")
        try:
            t0 = time.perf_counter()
            run_backend_chunked(
                "gcn",
                lambda sm: embed_molclr(sm, "gcn", device=args.device),
                drop_invalid=True,
            )
            _log(f"  GCN total: {time.perf_counter() - t0:.1f} s")
        except Exception as e:
            _log(f"  Skipped GCN: {e}")

    # 3) ECFP
    if "ecfp" in args.backends:
        _log("Computing ECFP embeddings (chunked)...")
        t0 = time.perf_counter()
        run_backend_chunked(
            "ecfp",
            lambda sm: embed_ecfp(sm, radius=args.radius, n_bits=args.n_bits),
            drop_invalid=False,
        )
        _log(f"  ECFP total: {time.perf_counter() - t0:.1f} s")

    _log("Done.")


if __name__ == "__main__":
    main()
