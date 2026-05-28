"""ECFP embedding and RankModel scoring (no MolCLR)."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import torch

from training.data import normalize_did

from .bundle import RunModelBundle, load_run_bundle
from .io import TaskBundle, TaskRow, load_task_bundle, write_json


def compute_ecfp_vector(smiles: str, radius: int = 2, n_bits: int = 2048) -> Optional[np.ndarray]:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator

    smi = (smiles or "").strip()
    if not smi:
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    fp = gen.GetFingerprint(mol)
    out = np.zeros((n_bits,), dtype=np.float32)
    DataStructs.ConvertToNumpyArray(fp, out)
    return out


def embed_row(
    row: TaskRow,
    ligand_lookup: Dict[str, np.ndarray],
    d_l: int,
) -> Tuple[Optional[np.ndarray], Optional[str]]:
    """Prefer smiles ECFP; fallback to did npz lookup."""
    if row.smiles.strip():
        vec = compute_ecfp_vector(row.smiles)
        if vec is not None:
            return vec, None
        return None, "invalid_smiles"
    if row.did.strip():
        key = normalize_did(row.did)
        if key and key in ligand_lookup:
            return ligand_lookup[key].astype(np.float32), None
        return None, "did_not_in_npz"
    return None, "missing_smiles_and_did"


def mean_embeddings(vectors: List[np.ndarray]) -> np.ndarray:
    if not vectors:
        raise ValueError("no embeddings for context")
    return np.stack(vectors, axis=0).mean(axis=0).astype(np.float32)


def score_and_rank_strings(
    model: torch.nn.Module,
    metal_vec: np.ndarray,
    context_vec: np.ndarray,
    candidate_embeddings: np.ndarray,
    candidate_keys: List[str],
    device: torch.device,
) -> List[Tuple[str, int, float]]:
    metal = torch.from_numpy(metal_vec).float().unsqueeze(0).to(device)
    context = torch.from_numpy(context_vec).float().unsqueeze(0).to(device)
    candidates = torch.from_numpy(candidate_embeddings).float().unsqueeze(0).to(device)
    with torch.no_grad():
        logits = model(metal, context, candidates)
    scores = logits[0].cpu().numpy()
    order = np.argsort(-scores)
    return [(candidate_keys[idx], rank, float(scores[idx])) for rank, idx in enumerate(order, start=1)]


def _default_single_output_dir(model_run_root: Path) -> Path:
    return model_run_root / "reports" / "model_after"


def recommend_single_model(
    model_run_root: str | Path,
    task_dir: str | Path,
    output_dir: str | Path | None = None,
    *,
    checkpoint: str | Path | None = None,
    device: str = "cpu",
    model_id: str | None = None,
) -> Dict[str, Any]:
    root = Path(model_run_root).resolve()
    out = Path(output_dir).resolve() if output_dir else _default_single_output_dir(root)
    out.mkdir(parents=True, exist_ok=True)

    bundle = load_run_bundle(root, checkpoint=checkpoint, device=device, model_id=model_id)
    task = load_task_bundle(task_dir)

    if task.metal not in bundle.metal_lookup:
        raise KeyError(f"metal {task.metal!r} not in metal embedding table")

    metal_vec = bundle.metal_lookup[task.metal]
    pos_vecs: List[np.ndarray] = []
    pos_errors: List[str] = []
    for row in task.positives:
        vec, err = embed_row(row, bundle.ligand_lookup, bundle.d_l)
        if vec is not None:
            pos_vecs.append(vec)
        elif err:
            pos_errors.append(f"{row.row_id}:{err}")

    if not pos_vecs:
        raise RuntimeError("no valid positive embeddings for context")
    context_vec = mean_embeddings(pos_vecs)

    valid_rows: List[TaskRow] = []
    valid_vecs: List[np.ndarray] = []
    invalid: List[Dict[str, str]] = []
    for row in task.candidates:
        vec, err = embed_row(row, bundle.ligand_lookup, bundle.d_l)
        if vec is not None:
            valid_rows.append(row)
            valid_vecs.append(vec)
        else:
            invalid.append({"candidate_id": row.row_id, "error": err or "unknown"})

    if not valid_vecs:
        raise RuntimeError("no valid candidate embeddings")

    cand_embs = np.stack(valid_vecs, axis=0)
    keys = [r.row_id for r in valid_rows]
    ranked = score_and_rank_strings(
        bundle.model, metal_vec, context_vec, cand_embs, keys, bundle.device
    )
    rank_map = {k: (rank, score) for k, rank, score in ranked}

    ranking_path = out / "ranking.csv"
    fields = [
        "rank", "score", "candidate_id", "did", "smiles", "name", "source",
        "is_known_positive", "is_known_negative", "valid", "error",
        "model_id", "model_run_root", "checkpoint_path",
    ]
    with open(ranking_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in task.candidates:
            rid = row.row_id
            if rid in rank_map:
                rank, score = rank_map[rid]
                valid = True
                err = ""
            else:
                rank, score = "", ""
                valid = False
                err = next((x["error"] for x in invalid if x["candidate_id"] == rid), "invalid")
            w.writerow({
                "rank": rank,
                "score": score,
                "candidate_id": rid,
                "did": row.did,
                "smiles": row.smiles,
                "name": row.name,
                "source": row.source,
                "is_known_positive": rid in task.positive_ids,
                "is_known_negative": rid in task.negative_ids,
                "valid": valid,
                "error": err,
                "model_id": bundle.model_id,
                "model_run_root": str(bundle.model_run_root),
                "checkpoint_path": str(bundle.ckpt_path),
            })

    report = {
        "status": "success",
        "context_source": "positives_mean",
        "task_id": task.task_id,
        "metal": task.metal,
        "model_id": bundle.model_id,
        "model_run_root": str(bundle.model_run_root),
        "checkpoint_path": str(bundle.ckpt_path),
        "device": str(bundle.device),
        "output_dir": str(out),
        "ranking_path": str(ranking_path),
        "n_candidates": len(task.candidates),
        "n_valid_candidates": len(valid_rows),
        "n_invalid_candidates": len(invalid),
        "positive_embedding_errors": pos_errors,
        "invalid_candidates": invalid,
    }
    write_json(out / "report.json", report)
    return report
