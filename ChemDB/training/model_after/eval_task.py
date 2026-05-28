"""Task-fit evaluation: held-out metrics and score margins."""

from __future__ import annotations

import csv
import itertools
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .bundle import load_run_bundle
from .io import TaskBundle, load_task_bundle, write_json
from .scoring import embed_row, mean_embeddings, recommend_single_model, score_and_rank_strings


def _hit_at(ranked_keys: List[str], heldout: List[str], top_k: int) -> float:
    top = set(ranked_keys[:top_k])
    if not heldout:
        return 0.0
    return sum(1 for h in heldout if h in top) / len(heldout)


def _heldout_metrics(
    bundle,
    task: TaskBundle,
    key_to_emb: Dict[str, np.ndarray],
    metal_vec: np.ndarray,
    max_k: Optional[int] = None,
) -> Tuple[Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], List[dict]]:
    """Leave-one-out on positive row_ids. Returns mrr, hit5,10,20,50, heldout_rows."""
    pos_keys = [r.row_id for r in task.positives if r.row_id in key_to_emb]
    if len(pos_keys) < 2:
        return None, None, None, None, None, []

    pool_ids: List[str] = []
    seen: set[str] = set()
    for r in list(task.candidates) + list(task.positives):
        if r.row_id in key_to_emb and r.row_id not in seen:
            seen.add(r.row_id)
            pool_ids.append(r.row_id)
    all_cand_keys = pool_ids
    upper_k = len(pos_keys) - 1
    if max_k is not None:
        upper_k = min(upper_k, max_k)

    context_rows: List[dict] = []
    heldout_rows: List[dict] = []

    for k in range(1, upper_k + 1):
        for ctx_tuple in itertools.combinations(pos_keys, k):
            ctx_set = set(ctx_tuple)
            heldout = [p for p in pos_keys if p not in ctx_set]
            ctx_vecs = [key_to_emb[i] for i in ctx_tuple]
            context_vec = mean_embeddings(ctx_vecs)
            cand_keys = [c for c in all_cand_keys if c not in ctx_set]
            if not cand_keys:
                continue
            cand_embs = np.stack([key_to_emb[c] for c in cand_keys], axis=0)
            ranked = score_and_rank_strings(
                bundle.model, metal_vec, context_vec, cand_embs, cand_keys, bundle.device
            )
            order = [x[0] for x in sorted(ranked, key=lambda t: t[1])]
            id_to_rank = {cid: rank for cid, rank, _ in ranked}
            ranks_for_mrr = [id_to_rank.get(h, len(order) + 1) for h in heldout]
            mrr = float(np.mean([1.0 / r for r in ranks_for_mrr])) if ranks_for_mrr else 0.0
            context_rows.append({
                "k": k,
                "heldout_count": len(heldout),
                "mrr": mrr,
                "hit_at_5": _hit_at(order, heldout, 5),
                "hit_at_10": _hit_at(order, heldout, 10),
                "hit_at_20": _hit_at(order, heldout, 20),
                "hit_at_50": _hit_at(order, heldout, 50),
            })
            for h in heldout:
                heldout_rows.append({
                    "heldout_positive_id": h,
                    "rank": id_to_rank.get(h, len(order) + 1),
                    "context_k": k,
                })

    if not context_rows:
        return None, None, None, None, None, heldout_rows

    return (
        float(np.mean([r["mrr"] for r in context_rows])),
        float(np.mean([r["hit_at_5"] for r in context_rows])),
        float(np.mean([r["hit_at_10"] for r in context_rows])),
        float(np.mean([r["hit_at_20"] for r in context_rows])),
        float(np.mean([r["hit_at_50"] for r in context_rows])),
        heldout_rows,
    )


def evaluate_task_fit(
    model_run_root: str | Path,
    task_dir: str | Path,
    output_dir: str | Path,
    *,
    checkpoint: str | Path | None = None,
    device: str = "cpu",
    model_id: str | None = None,
) -> Dict[str, Any]:
    out = Path(output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    bundle = load_run_bundle(model_run_root, checkpoint=checkpoint, device=device, model_id=model_id)
    task = load_task_bundle(task_dir)

    rank_report = recommend_single_model(
        model_run_root=bundle.model_run_root,
        task_dir=task.task_dir,
        output_dir=out,
        checkpoint=bundle.ckpt_path,
        device=str(bundle.device),
        model_id=bundle.model_id,
    )

    warnings: List[str] = []
    key_to_emb: Dict[str, np.ndarray] = {}
    for row in list(task.positives) + list(task.candidates):
        if row.row_id in key_to_emb:
            continue
        vec, err = embed_row(row, bundle.ligand_lookup, bundle.d_l)
        if vec is not None:
            key_to_emb[row.row_id] = vec
        elif err and row in task.candidates:
            warnings.append(f"candidate {row.row_id}: {err}")

    metal_vec = bundle.metal_lookup.get(task.metal)
    if metal_vec is None:
        raise KeyError(f"metal {task.metal!r} not in embedding table")

    mrr, h5, h10, h20, h50, heldout_rows = _heldout_metrics(bundle, task, key_to_emb, metal_vec)
    if len(task.positives) < 2:
        warnings.append("positives < 2: held-out mrr/hit@k set to null")
        status = "partial"
    elif mrr is None:
        warnings.append("held-out evaluation produced no cases")
        status = "partial"
    else:
        status = "success"

    id_to_score: Dict[str, float] = {}
    ranking_path = Path(rank_report["ranking_path"])
    if ranking_path.is_file():
        with open(ranking_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row.get("valid") in (True, "True", "true", "1"):
                    try:
                        id_to_score[row["candidate_id"]] = float(row["score"])
                    except (TypeError, ValueError):
                        pass

    pos_scores = [id_to_score[k] for k in task.positive_ids if k in id_to_score]
    neg_scores = [id_to_score[k] for k in task.negative_ids if k in id_to_score]

    positive_score_mean = float(np.mean(pos_scores)) if pos_scores else None
    negative_score_mean = float(np.mean(neg_scores)) if neg_scores else None
    pos_neg_margin = None
    if positive_score_mean is not None and negative_score_mean is not None:
        pos_neg_margin = positive_score_mean - negative_score_mean
    elif not task.negatives:
        warnings.append("no negatives: negative_score_mean and pos_neg_margin are null")

    metrics = {
        "model_id": bundle.model_id,
        "task_id": task.task_id,
        "metal": task.metal,
        "num_candidates": len(task.candidates),
        "num_positives": len(task.positives),
        "num_negatives": len(task.negatives),
        "mrr": mrr,
        "hit_at_5": h5,
        "hit_at_10": h10,
        "hit_at_20": h20,
        "hit_at_50": h50,
        "positive_score_mean": positive_score_mean,
        "negative_score_mean": negative_score_mean,
        "pos_neg_margin": pos_neg_margin,
        "status": status,
        "warnings": warnings,
        "context_source": "positives_mean",
        "model_run_root": str(bundle.model_run_root),
        "checkpoint_path": str(bundle.ckpt_path),
        "ranking_path": str(ranking_path),
    }
    write_json(out / "metrics.json", metrics)
    write_json(out / "evaluation_manifest.json", {
        "task_id": task.task_id,
        "model_id": bundle.model_id,
        "output_dir": str(out),
        "n_heldout_rows": len(heldout_rows or []),
    })

    if heldout_rows:
        held_path = out / "heldout_ranking.csv"
        with open(held_path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(heldout_rows[0].keys()))
            w.writeheader()
            w.writerows(heldout_rows)

    return metrics
