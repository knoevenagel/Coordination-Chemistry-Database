"""Phase 1D CLI: model-after recommend / evaluate / multi-model selection."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional

from .model_after_paths import default_output_dir

CHEMDB_ROOT = Path(__file__).resolve().parents[2] / "ChemDB"


def _ensure_chemdb_on_path() -> None:
    for p in (str(CHEMDB_ROOT), str(CHEMDB_ROOT / "src")):
        if p not in sys.path:
            sys.path.insert(0, p)


def _write_manifest(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)


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
