"""Multi model_run_root evaluation from model_runs.csv."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Dict, List, Optional

from .bundle import load_run_bundle
from .eval_task import evaluate_task_fit
from .io import write_json
from .scoring import recommend_single_model


def _selection_score(row: Dict[str, Any]) -> Optional[float]:
    mrr = row.get("mrr")
    margin = row.get("pos_neg_margin")
    if mrr is not None:
        return float(mrr)
    if margin is not None:
        return float(margin)
    return None


def load_model_runs_table(path: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for raw in reader:
            run_root = (raw.get("model_run_root") or "").strip()
            if not run_root:
                continue
            rows.append({
                "model_id": (raw.get("model_id") or Path(run_root).name).strip(),
                "model_run_root": run_root,
                "checkpoint_path": (raw.get("checkpoint_path") or "").strip(),
                "model_name": (raw.get("model_name") or "").strip(),
                "notes": (raw.get("notes") or "").strip(),
            })
    return rows


def evaluate_models(
    model_runs_file: str | Path,
    task_dir: str | Path,
    output_dir: str | Path,
    *,
    device: str = "cpu",
) -> Dict[str, Any]:
    out = Path(output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    models_root = out / "models"
    models_root.mkdir(parents=True, exist_ok=True)

    entries = load_model_runs_table(Path(model_runs_file))
    if not entries:
        raise ValueError(f"no models in {model_runs_file}")

    summary_rows: List[Dict[str, Any]] = []
    errors: List[Dict[str, str]] = []

    for ent in entries:
        mid = ent["model_id"]
        run_root = Path(ent["model_run_root"]).resolve()
        ckpt = ent["checkpoint_path"] or None
        if ckpt:
            ckpt = str(Path(ckpt).resolve())
        model_out = models_root / mid
        row: Dict[str, Any] = {
            "model_id": mid,
            "model_name": ent.get("model_name") or mid,
            "model_run_root": str(run_root),
            "checkpoint_path": ckpt or str(run_root / "training" / "ckpts" / "best.pt"),
            "embedding_backend": "ecfp",
            "status": "failed",
            "mrr": None,
            "hit_at_5": None,
            "hit_at_10": None,
            "hit_at_20": None,
            "hit_at_50": None,
            "positive_score_mean": None,
            "negative_score_mean": None,
            "pos_neg_margin": None,
            "selection_score": None,
            "rank_among_models": None,
            "error": "",
        }
        try:
            if not run_root.is_dir():
                raise FileNotFoundError(f"model_run_root missing: {run_root}")
            metrics = evaluate_task_fit(
                run_root,
                task_dir,
                model_out,
                checkpoint=ckpt,
                device=device,
                model_id=mid,
            )
            row.update({
                "status": metrics.get("status", "success"),
                "mrr": metrics.get("mrr"),
                "hit_at_5": metrics.get("hit_at_5"),
                "hit_at_10": metrics.get("hit_at_10"),
                "hit_at_20": metrics.get("hit_at_20"),
                "hit_at_50": metrics.get("hit_at_50"),
                "positive_score_mean": metrics.get("positive_score_mean"),
                "negative_score_mean": metrics.get("negative_score_mean"),
                "pos_neg_margin": metrics.get("pos_neg_margin"),
                "error": "; ".join(metrics.get("warnings") or []),
            })
            row["selection_score"] = _selection_score(row)
        except Exception as exc:
            row["error"] = str(exc)
            errors.append({"model_id": mid, "error": str(exc)})

        summary_rows.append(row)

    scored = [r for r in summary_rows if r.get("selection_score") is not None]
    scored.sort(key=lambda r: float(r["selection_score"]), reverse=True)
    for i, r in enumerate(scored, start=1):
        r["rank_among_models"] = i
    for r in summary_rows:
        if r.get("rank_among_models") is None:
            r["rank_among_models"] = None

    summary_path = out / "model_selection_summary.csv"
    fields = [
        "model_id", "model_name", "model_run_root", "checkpoint_path", "embedding_backend",
        "status", "mrr", "hit_at_5", "hit_at_10", "hit_at_20", "hit_at_50",
        "positive_score_mean", "negative_score_mean", "pos_neg_margin",
        "selection_score", "rank_among_models", "error",
    ]
    with open(summary_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(summary_rows)

    best: Optional[Dict[str, Any]] = None
    if scored:
        best = scored[0]
        best_payload = {
            "model_id": best["model_id"],
            "model_run_root": best["model_run_root"],
            "checkpoint_path": best["checkpoint_path"],
            "selection_score": best["selection_score"],
            "selection_reason": "highest selection_score (mrr, else pos_neg_margin)",
            "metrics_path": str(models_root / best["model_id"] / "metrics.json"),
            "ranking_path": str(models_root / best["model_id"] / "ranking.csv"),
        }
    else:
        best_payload = {
            "model_id": None,
            "selection_reason": "no model with valid selection_score",
            "errors": errors,
        }
    write_json(out / "best_model.json", best_payload)

    manifest = {
        "task_dir": str(Path(task_dir).resolve()),
        "output_dir": str(out),
        "n_models": len(entries),
        "n_success_with_score": len(scored),
        "summary_path": str(summary_path),
        "errors": errors,
    }
    write_json(out / "model_selection_manifest.json", manifest)
    return manifest
