#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert database identifiers (DIDs) into embeddings using MolCLR (GIN/GCN)
or RDKit/ECFP fingerprints.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from functools import lru_cache
from typing import Dict, Optional, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from molclr_api import get_model, extract_embeddings  # type: ignore

DEFAULT_DID_CSV = os.path.join(PROJECT_ROOT, "tmp", "repaired_ligand_data.csv")
NEO4J_PREFIXES = [
    "L0_", "L1_", "L2_", "L3_", "L4_", "L5_", "M_",
    "METAL_", "COMPLEX_", "LIGAND_", "REPAIREDLIGAND_", "FRAGMENT_", "IRL_",
]


def normalize_did(did: str) -> str:
    text = did.strip()
    upper = text.upper()
    for prefix in NEO4J_PREFIXES:
        if upper.startswith(prefix):
            return text[len(prefix):]
    return text


@lru_cache(maxsize=4)
def load_did_smiles_map(csv_path: str) -> Dict[str, str]:
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"DID-SMILES CSV not found: {csv_path}")

    did_cols = {"ligand_new_did", "did", "source_did"}
    smiles_cols = {"new_smiles", "smiles", "old_smiles"}
    mapping: Dict[str, str] = {}

    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        did_col = next((c for c in reader.fieldnames or [] if c in did_cols), None)
        smiles_col = next((c for c in reader.fieldnames or [] if c in smiles_cols), None)
        if not did_col or not smiles_col:
            raise ValueError(f"CSV {csv_path} lacks DID/SMILES columns {reader.fieldnames}")

        for row in reader:
            did = row.get(did_col, "").strip()
            smiles = row.get(smiles_col, "").strip()
            if did and smiles:
                mapping[did] = smiles

    if not mapping:
        raise ValueError(f"No DID-SMILES entries loaded from {csv_path}")
    return mapping


def get_smiles_from_did(did: str, csv_path: Optional[str]) -> str:
    csv_path = csv_path or DEFAULT_DID_CSV
    did_map = load_did_smiles_map(csv_path)
    did_norm = normalize_did(did)
    smiles = did_map.get(did_norm)
    if not smiles:
        raise KeyError(f"DID {did} not found in {csv_path}")
    return smiles


def embed_with_molclr(smiles: str, model_type: str) -> np.ndarray:
    model = get_model(model_type)
    embeddings, valid_indices, invalid_smiles = extract_embeddings(model, [smiles], device="cpu")
    if not valid_indices:
        raise ValueError(f"MolCLR ({model_type}) failed for SMILES: {smiles}")
    return embeddings[0]


def embed_with_rdkit(smiles: str, radius: int = 2, n_bits: int = 2048) -> np.ndarray:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    bitvect = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)
    array = np.zeros((n_bits,), dtype=np.float32)
    DataStructs.ConvertToNumpyArray(bitvect, array)
    return array


def did_to_embedding(
    did: str,
    backend: str = "molclr_gin",
    did_csv: Optional[str] = None,
    radius: int = 2,
    n_bits: int = 2048,
) -> Tuple[np.ndarray, Dict[str, str]]:
    smiles = get_smiles_from_did(did, did_csv)
    backend_key = backend.lower()

    if backend_key == "molclr_gin":
        vector = embed_with_molclr(smiles, "gin")
    elif backend_key == "molclr_gcn":
        vector = embed_with_molclr(smiles, "gcn")
    elif backend_key in {"rdkit_fp", "ecfp"}:
        vector = embed_with_rdkit(smiles, radius=radius, n_bits=n_bits)
    else:
        raise ValueError(f"Unsupported backend {backend}")

    meta = {
        "did": did,
        "normalized_did": normalize_did(did),
        "smiles": smiles,
        "backend": backend_key,
        "shape": str(vector.shape),
        "dtype": str(vector.dtype),
    }
    return vector, meta


def save_embedding(embedding: np.ndarray, meta: Dict[str, str], output: Optional[str], fmt: str):
    fmt = fmt.lower()
    if fmt == "npy":
        if not output:
            raise ValueError("npy output requires --output path")
        np.save(output, embedding)
        print(f"Saved embedding to {output}")
        return

    payload = {"meta": meta, "embedding": embedding.tolist()}
    text = json.dumps(payload, ensure_ascii=False, indent=2)
    if output:
        with open(output, "w", encoding="utf-8") as f:
            f.write(text + "\n")
        print(f"Saved JSON to {output}")
    else:
        print(text)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert DID to embedding (MolCLR or RDKit/ECFP).")
    parser.add_argument("--did", required=True, help="Target DID (supports Neo4j prefixes)")
    parser.add_argument(
        "--backend",
        default="molclr_gin",
        choices=["molclr_gin", "molclr_gcn", "rdkit_fp", "ecfp"],
        help="Embedding backend",
    )
    parser.add_argument(
        "--did-csv",
        default=None,
        help=f"DID-SMILES mapping CSV path (default: {DEFAULT_DID_CSV})",
    )
    parser.add_argument("--radius", type=int, default=2, help="Radius for RDKit/ECFP backends")
    parser.add_argument("--n-bits", type=int, default=2048, help="Fingerprint size for RDKit/ECFP")
    parser.add_argument("--output", default=None, help="Optional output file path")
    parser.add_argument("--format", default="json", choices=["json", "npy"], help="Output format")
    return parser.parse_args()


def main():
    args = parse_args()
    embedding, meta = did_to_embedding(
        did=args.did,
        backend=args.backend,
        did_csv=args.did_csv,
        radius=args.radius,
        n_bits=args.n_bits,
    )
    save_embedding(embedding, meta, args.output, fmt=args.format)


if __name__ == "__main__":
    main()
