"""Workspace paths for shared model-after tasks and results."""

from __future__ import annotations

from pathlib import Path

WEB_ROOT = Path(__file__).resolve().parents[2]
WORKSPACE_ROOT = WEB_ROOT / "workspace"
MODEL_AFTER_TASKS_ROOT = WORKSPACE_ROOT / "model_after_tasks"
MODEL_AFTER_RESULTS_ROOT = WORKSPACE_ROOT / "model_after_results"

RESULT_FILENAMES = frozenset({
    "ranking.csv",
    "report.json",
    "metrics.json",
    "evaluation_manifest.json",
    "heldout_ranking.csv",
    "recommend_manifest.json",
    "evaluate_model_manifest.json",
    "evaluate_models_manifest.json",
    "model_selection_summary.csv",
    "best_model.json",
    "model_selection_manifest.json",
})


def task_id_from_dir(task_dir: Path) -> str:
    return task_dir.resolve().name


def batch_id_from_run(model_run_root: Path) -> str:
    """e.g. workspace/projects/p001/runs/run_demo -> p001_run_demo."""
    parts = model_run_root.resolve().parts
    if "projects" in parts and "runs" in parts:
        pi = parts.index("projects")
        ri = parts.index("runs")
        if pi + 1 < len(parts) and ri + 1 < len(parts):
            return f"{parts[pi + 1]}_{parts[ri + 1]}"
    return model_run_root.resolve().name


def default_output_dir(
    task_dir: Path,
    *,
    model_run_root: Path | None = None,
    batch_id: str | None = None,
    command: str = "",
) -> Path:
    """Default: workspace/model_after_results/{task_id}/{batch_id}/."""
    tid = task_id_from_dir(task_dir)
    if batch_id:
        bid = batch_id
    elif model_run_root is not None:
        bid = batch_id_from_run(model_run_root)
        if command == "evaluate-model":
            bid = f"{bid}_eval"
    else:
        bid = "batch_default"
    return MODEL_AFTER_RESULTS_ROOT / tid / bid
