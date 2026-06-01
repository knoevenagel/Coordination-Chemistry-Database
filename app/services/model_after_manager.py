"""Phase 1D CLI: model-after recommend / evaluate / multi-model selection."""

from __future__ import annotations

import argparse
import json
import logging
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from .model_after_paths import default_output_dir
from . import job_manager

CHEMDB_ROOT = Path(__file__).resolve().parents[2] / "ChemDB"


def _ensure_chemdb_on_path() -> None:
    for p in (str(CHEMDB_ROOT), str(CHEMDB_ROOT / "src")):
        if p not in sys.path:
            sys.path.insert(0, p)


def _write_manifest(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)


def _quote(parts: List[str]) -> str:
    return " ".join(shlex.quote(p) for p in parts)


def _batch_id_default(prefix: str = "batch") -> str:
    return f"{prefix}_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"


def create_model_runs_file(
    *,
    workspace_root: str | Path,
    task_id: str,
    models: List[Dict[str, Any]],
    batch_id: str,
) -> Path:
    ws = Path(workspace_root).resolve()
    out = ws / "model_after_results" / task_id / batch_id / "_inputs"
    out.mkdir(parents=True, exist_ok=True)
    path = out / "model_runs.csv"
    lines = ["model_id,model_run_root,checkpoint_path,model_name,notes"]
    for row in models:
        mid = str(row.get("model_id") or "")
        project_id = str(row.get("project_id") or "")
        run_id = str(row.get("run_id") or "")
        run_root = str((ws / "projects" / project_id / "runs" / run_id).resolve())
        ckpt = str(row.get("checkpoint_path") or "")
        model_name = str(row.get("checkpoint_stem") or mid)
        lines.append(f"{mid},{run_root},{ckpt},{model_name},")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def submit_evaluate_models_job(
    *,
    db_path: str | Path,
    workspace_root: str | Path,
    task_id: str,
    models: List[Dict[str, Any]],
    batch_id: Optional[str] = None,
    device: str = "cpu",
) -> Dict[str, Any]:
    if not models:
        raise ValueError("no models selected")
    ws = Path(workspace_root).resolve()
    bid = batch_id or _batch_id_default("batch")
    task_dir = ws / "model_after_tasks" / task_id
    if not task_dir.is_dir():
        raise FileNotFoundError(f"task directory not found: {task_dir}")
    output_dir = ws / "model_after_results" / task_id / bid
    model_runs_file = create_model_runs_file(workspace_root=ws, task_id=task_id, models=models, batch_id=bid)

    eval_cmd = [
        sys.executable,
        "-m",
        "app.services.model_after_manager",
        "evaluate-models",
        "--model-runs-file",
        str(model_runs_file),
        "--task-dir",
        str(task_dir),
        "--output-dir",
        str(output_dir),
        "--device",
        device,
        "--batch-id",
        bid,
    ]
    rebuild_cmd = [
        sys.executable,
        "-m",
        "app.storage.cli",
        "rebuild-index",
        "--workspace-root",
        str(ws),
        "--db-path",
        str(Path(db_path).resolve()),
    ]
    command = ["bash", "-lc", f"{_quote(eval_cmd)} && {_quote(rebuild_cmd)}"]
    job = job_manager.launch_subprocess_job_async(
        db_path,
        job_type="evaluate_models",
        title=f"Evaluate models for task {task_id}/{bid}",
        command=command,
        workspace_root=ws,
        task_id=task_id,
        batch_id=bid,
        cwd=ws,
    )
    return {
        "job": job,
        "task_id": task_id,
        "batch_id": bid,
        "output_dir": str(output_dir.resolve()),
        "model_runs_file": str(model_runs_file.resolve()),
    }


def submit_recommend_job(
    *,
    db_path: str | Path,
    workspace_root: str | Path,
    task_dir: str | Path,
    model_run_root: str | Path,
    model_id: Optional[str] = None,
    batch_id: Optional[str] = None,
    device: str = "cpu",
) -> Dict[str, Any]:
    ws = Path(workspace_root).resolve()
    task_dir_p = Path(task_dir).resolve()
    if not task_dir_p.is_dir():
        raise FileNotFoundError(f"task directory not found: {task_dir_p}")
    model_run_root_p = Path(model_run_root).resolve()
    bid = batch_id or _batch_id_default("recommend")
    output_dir = ws / "recommendation_results" / task_dir_p.name / bid
    rec_cmd = [
        sys.executable,
        "-m",
        "app.services.model_after_manager",
        "recommend",
        "--model-run-root",
        str(model_run_root_p),
        "--task-dir",
        str(task_dir_p),
        "--output-dir",
        str(output_dir),
        "--device",
        device,
        "--batch-id",
        bid,
    ]
    if model_id:
        rec_cmd.extend(["--model-id", model_id])
    command = ["bash", "-lc", _quote(rec_cmd)]
    job = job_manager.launch_subprocess_job_async(
        db_path,
        job_type="recommend",
        title=f"Recommend {task_dir_p.name}/{bid}",
        command=command,
        workspace_root=ws,
        task_id=task_dir_p.name,
        batch_id=bid,
        model_id=model_id,
        cwd=ws,
    )
    return {
        "job": job,
        "task_dir": str(task_dir_p),
        "batch_id": bid,
        "output_dir": str(output_dir.resolve()),
    }


def _resolve_output(
    task_dir: Path,
    *,
    model_run_root: Path | None = None,
    output_dir: str | Path | None = None,
    batch_id: str | None = None,
    command: str = "",
) -> Path:
    if output_dir:
        return Path(output_dir).resolve()
    return default_output_dir(
        task_dir,
        model_run_root=model_run_root,
        batch_id=batch_id,
        command=command,
    )


def cmd_recommend(args: argparse.Namespace) -> int:
    _ensure_chemdb_on_path()
    from training.model_after import recommend_single_model

    run_root = Path(args.model_run_root).resolve()
    task_dir = Path(args.task_dir).resolve()
    out = _resolve_output(
        task_dir,
        model_run_root=run_root,
        output_dir=args.output_dir,
        batch_id=args.batch_id,
        command="recommend",
    )

    logging.info("recommend model_run_root=%s task_dir=%s output_dir=%s", run_root, task_dir, out)
    report = recommend_single_model(
        run_root,
        task_dir,
        out,
        checkpoint=args.checkpoint,
        device=args.device,
        model_id=args.model_id,
    )
    manifest = {
        "command": "recommend",
        "status": report.get("status", "success"),
        "model_run_root": str(run_root),
        "task_dir": str(task_dir),
        "output_dir": str(out),
        "device": args.device,
        "finished_at": datetime.now(timezone.utc).isoformat(),
        "report": report,
    }
    _write_manifest(out / "recommend_manifest.json", manifest)
    logging.info("ranking: %s", report.get("ranking_path"))
    return 0


def cmd_evaluate_model(args: argparse.Namespace) -> int:
    _ensure_chemdb_on_path()
    from training.model_after import evaluate_task_fit

    run_root = Path(args.model_run_root).resolve()
    task_dir = Path(args.task_dir).resolve()
    out = _resolve_output(
        task_dir,
        model_run_root=run_root,
        output_dir=args.output_dir,
        batch_id=args.batch_id,
        command="evaluate-model",
    )

    logging.info("evaluate-model model_run_root=%s task_dir=%s output_dir=%s", run_root, task_dir, out)
    metrics = evaluate_task_fit(
        run_root,
        task_dir,
        out,
        checkpoint=args.checkpoint,
        device=args.device,
        model_id=args.model_id,
    )
    manifest = {
        "command": "evaluate-model",
        "status": metrics.get("status", "success"),
        "model_run_root": str(run_root),
        "task_dir": str(task_dir),
        "output_dir": str(out),
        "device": args.device,
        "finished_at": datetime.now(timezone.utc).isoformat(),
        "metrics": metrics,
    }
    _write_manifest(out / "evaluate_model_manifest.json", manifest)
    logging.info("metrics: %s", out / "metrics.json")
    return 0


def cmd_evaluate_models(args: argparse.Namespace) -> int:
    _ensure_chemdb_on_path()
    from training.model_after import evaluate_models

    task_dir = Path(args.task_dir).resolve()
    runs_file = Path(args.model_runs_file).resolve()
    out = _resolve_output(
        task_dir,
        output_dir=args.output_dir,
        batch_id=args.batch_id or "batch_default",
        command="evaluate-models",
    )

    logging.info("evaluate-models task_dir=%s output_dir=%s", task_dir, out)
    result = evaluate_models(
        runs_file,
        task_dir,
        out,
        device=args.device,
    )
    manifest = {
        "command": "evaluate-models",
        "task_dir": str(task_dir),
        "model_runs_file": str(runs_file),
        "output_dir": str(out),
        "device": args.device,
        "finished_at": datetime.now(timezone.utc).isoformat(),
        "result": result,
    }
    _write_manifest(out / "evaluate_models_manifest.json", manifest)
    logging.info("summary: %s", result.get("summary_path"))
    return 0


def _add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument("--device", default="cpu", help="cpu, cuda, or auto (default: cpu)")
    p.add_argument(
        "--batch-id",
        default=None,
        help="Batch folder under workspace/model_after_results/{task_id}/ (default: derived or batch_default)",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="ChemDB Phase 1D model-after manager")
    parser.add_argument("-v", "--verbose", action="store_true")
    sub = parser.add_subparsers(dest="command", required=True)

    p_rec = sub.add_parser("recommend", help="Rank candidates for one model run")
    p_rec.add_argument("--model-run-root", required=True)
    p_rec.add_argument("--task-dir", required=True)
    p_rec.add_argument("--output-dir", default=None)
    p_rec.add_argument("--checkpoint", default=None)
    p_rec.add_argument("--model-id", default=None)
    _add_common(p_rec)

    p_ev = sub.add_parser("evaluate-model", help="Task-fit metrics for one model run")
    p_ev.add_argument("--model-run-root", required=True)
    p_ev.add_argument("--task-dir", required=True)
    p_ev.add_argument("--output-dir", default=None)
    p_ev.add_argument("--checkpoint", default=None)
    p_ev.add_argument("--model-id", default=None)
    _add_common(p_ev)

    p_ms = sub.add_parser("evaluate-models", help="Compare multiple model_run_root entries")
    p_ms.add_argument("--model-runs-file", required=True)
    p_ms.add_argument("--task-dir", required=True)
    p_ms.add_argument("--output-dir", default=None)
    _add_common(p_ms)

    return parser


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(message)s",
    )
    if args.command == "recommend":
        return cmd_recommend(args)
    if args.command == "evaluate-model":
        return cmd_evaluate_model(args)
    if args.command == "evaluate-models":
        return cmd_evaluate_models(args)
    parser.error(f"unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
