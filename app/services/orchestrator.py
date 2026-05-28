"""Run-based pipeline orchestration (Phase 1A core + Phase 1C training)."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from .run_context import RunContext
from .step_registry import (
    PIPELINE_ORDER,
    STEP_REGISTRY,
    TRAINING_PIPELINE_ORDER,
    StepSpec,
    get_step,
    pipeline_through,
)
from . import training_manager, ga_registry_manager

WEB_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CHEMDB_REPO = WEB_ROOT / "ChemDB"

SEED_DATA_FILES = ("metal_list.txt", "p_elements_list.txt")
METAL_EMBEDDING_SEED = "data/metal_embedding/element_features.csv"


def _run_subprocess(
    cmd: List[str],
    *,
    cwd: str,
    env: Dict[str, str],
    log_f,
    stream_to_console: bool,
) -> int:
    if stream_to_console:
        proc = subprocess.Popen(
            cmd,
            cwd=cwd,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            bufsize=1,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            log_f.write(line)
            log_f.flush()
            sys.stdout.write(line)
            sys.stdout.flush()
        return proc.wait()
    return subprocess.run(
        cmd,
        cwd=cwd,
        env=env,
        stdout=log_f,
        stderr=subprocess.STDOUT,
    ).returncode


def should_stream_step_output(stream_output: Optional[bool] = None) -> bool:
    if stream_output is not None:
        return stream_output
    return os.environ.get("CHEMDB_STREAM_STEP_OUTPUT", "").lower() in ("1", "true", "yes")


def chemdb_repo_root() -> Path:
    env = os.environ.get("CHEMDB_REPO_ROOT")
    if env:
        return Path(env).resolve()
    return DEFAULT_CHEMDB_REPO.resolve()


def resolve_run_path(ctx: RunContext, rel: str) -> Path:
    return (ctx.run_root / rel).resolve()


def expand_command(spec: StepSpec, ctx: RunContext, extra_argv: Optional[List[str]] = None) -> List[str]:
    mapping = {
        "{data_dir}": str(ctx.data_dir),
        "{tmp_dir}": str(ctx.tmp_dir),
        "{training_dir}": str(ctx.training_dir),
    }
    expanded: List[str] = []
    for token in spec.command_argv:
        out = token
        for k, v in mapping.items():
            out = out.replace(k, v)
        expanded.append(out)
    if extra_argv:
        expanded.extend(extra_argv)
    return expanded


def _pythonpath_env(repo: Path, env: Dict[str, str]) -> Dict[str, str]:
    parts = [str(repo), str(repo / "src")]
    prev = env.get("PYTHONPATH", "")
    if prev:
        parts.append(prev)
    env = dict(env)
    env["PYTHONPATH"] = os.pathsep.join(parts)
    return env


def _resolve_script_path(repo: Path, spec: StepSpec) -> Optional[Path]:
    if spec.script_kind == "python_module":
        return None
    if spec.script_kind == "repo_script":
        return repo / spec.script
    return repo / "src" / spec.script


def build_command(spec: StepSpec, ctx: RunContext, extra_argv: Optional[List[str]] = None) -> List[str]:
    argv = expand_command(spec, ctx, extra_argv)
    repo = chemdb_repo_root()
    if spec.script_kind == "python_module":
        return [sys.executable, "-m", argv[0]] + argv[1:]
    script_path = _resolve_script_path(repo, spec)
    assert script_path is not None
    return [sys.executable, str(script_path)] + argv[1:]


def check_inputs(ctx: RunContext, spec: StepSpec) -> Dict[str, Any]:
    missing: List[str] = []
    checked: List[str] = []
    for rel in spec.required_inputs:
        path = resolve_run_path(ctx, rel)
        checked.append(str(path))
        if rel == "data/pubchem" or rel.endswith("/pubchem"):
            if not path.is_dir():
                missing.append(rel)
            elif not any(path.glob("*.csv")):
                missing.append(f"{rel} (no *.csv)")
        elif path.is_dir():
            if not any(path.iterdir()):
                missing.append(f"{rel} (empty dir)")
        elif not path.is_file():
            missing.append(rel)
    ok = len(missing) == 0
    return {"ok": ok, "missing": missing, "checked": checked}


def check_outputs(
    ctx: RunContext,
    spec: StepSpec,
    *,
    mandatory_only: bool = False,
) -> Dict[str, Any]:
    missing: List[str] = []
    present: List[str] = []
    for rel in spec.expected_outputs:
        path = resolve_run_path(ctx, rel)
        if path.is_file():
            present.append(rel)
        else:
            missing.append(rel)
    optional_missing: List[str] = []
    if not mandatory_only:
        for rel in spec.optional_outputs:
            path = resolve_run_path(ctx, rel)
            if not path.is_file():
                optional_missing.append(rel)
    ok = len(missing) == 0
    return {
        "ok": ok,
        "missing": missing,
        "present": present,
        "optional_missing": optional_missing,
    }


def _paths_under_run(ctx: RunContext, paths: List[str]) -> bool:
    for p in paths:
        if not ctx.resolve_under_run(Path(p)):
            return False
    return True


def _path_entry(path: Path) -> Dict[str, Any]:
    exists = path.exists()
    size = path.stat().st_size if exists and path.is_file() else None
    return {"path": str(path.resolve()), "exists": exists, "size_bytes": size}


def _manifest_io_lists(
    ctx: RunContext,
    spec: StepSpec,
    *,
    include_optional_outputs: bool = False,
) -> tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    inputs = [_path_entry(resolve_run_path(ctx, rel)) for rel in spec.required_inputs]
    out_rels = list(spec.expected_outputs)
    if include_optional_outputs:
        out_rels += list(spec.optional_outputs)
    outputs = [_path_entry(resolve_run_path(ctx, rel)) for rel in out_rels]
    return inputs, outputs


def write_manifest(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)


def _training_integration_mode() -> bool:
    return os.environ.get("CHEMDB_TRAINING_INTEGRATION", "").lower() in ("1", "true", "yes")


def run_in_process_step(ctx: RunContext, spec: StepSpec) -> Dict[str, Any]:
    if spec.step_id == "prepare_training_index":
        result = training_manager.prepare_training_index(ctx)
    elif spec.step_id == "prepare_training_config":
        result = training_manager.prepare_training_config(
            ctx, integration=_training_integration_mode()
        )
    elif spec.step_id == "require_bound_ga":
        result = ga_registry_manager.require_bound_ga(ctx)
    else:
        return {"ok": False, "error": f"unknown in-process step: {spec.step_id}"}
    return result


def init_run(
    ctx: RunContext,
    *,
    pubchem_source: Optional[Path] = None,
) -> Dict[str, Any]:
    ctx.create_dirs()
    repo = chemdb_repo_root()
    copied: List[str] = []
    for name in SEED_DATA_FILES:
        src = repo / "data" / name
        dst = ctx.data_dir / name
        if src.is_file():
            shutil.copy2(src, dst)
            copied.append(name)

    metal_error: Optional[str] = None
    metal_src = repo / "data" / "metal_embedding" / "element_features.csv"
    metal_dst = ctx.data_dir / "metal_embedding" / "element_features.csv"
    if metal_src.is_file():
        metal_dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(metal_src, metal_dst)
        copied.append("metal_embedding/element_features.csv")
    else:
        metal_error = f"Missing metal embedding seed: {metal_src}"

    pubchem_files: List[str] = []
    if pubchem_source is not None:
        src = pubchem_source.resolve()
        pubchem_dir = ctx.data_dir / "pubchem"
        pubchem_dir.mkdir(parents=True, exist_ok=True)
        if src.is_file():
            if src.suffix.lower() != ".csv":
                raise ValueError(f"--pubchem-source must be a CSV file or directory: {src}")
            dst = pubchem_dir / src.name
            shutil.copy2(src, dst)
            pubchem_files.append(dst.name)
        elif src.is_dir():
            for csv in sorted(src.glob("*.csv")):
                shutil.copy2(csv, pubchem_dir / csv.name)
                pubchem_files.append(csv.name)
        else:
            raise FileNotFoundError(f"--pubchem-source not found: {src}")

    payload: Dict[str, Any] = {
        "action": "init-run",
        "run_root": str(ctx.run_root),
        "chemdb_repo_root": str(repo),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "seed_files_copied": copied,
        "pubchem_files": pubchem_files,
    }
    if metal_error:
        payload["metal_embedding_seed_error"] = metal_error
        payload["status"] = "failed"
    else:
        payload["status"] = "success"
    write_manifest(ctx.manifest_dir / "run.json", payload)
    if metal_error:
        raise FileNotFoundError(metal_error)
    return payload


def _finalize_manifest(
    manifest: Dict[str, Any],
    *,
    ctx: RunContext,
    spec: StepSpec,
    started_ts: float,
    exit_code: Optional[int] = None,
    error: Optional[str] = None,
    include_optional_outputs: bool = False,
) -> Dict[str, Any]:
    finished = datetime.now(timezone.utc).isoformat()
    manifest["finished_at"] = finished
    manifest["duration_seconds"] = round(time.perf_counter() - started_ts, 3)
    if exit_code is not None:
        manifest["exit_code"] = exit_code
    out = check_outputs(ctx, spec, mandatory_only=True)
    manifest["output_check"] = out
    inputs, outputs = _manifest_io_lists(
        ctx, spec, include_optional_outputs=include_optional_outputs
    )
    manifest["inputs"] = inputs
    manifest["outputs"] = outputs
    all_paths = [e["path"] for e in inputs + outputs]
    manifest["paths_under_run_root"] = _paths_under_run(ctx, all_paths)
    if error:
        manifest["status"] = "failed"
        manifest["error"] = error
    elif exit_code is not None and exit_code != 0:
        manifest["status"] = "failed"
        manifest["error"] = f"process exit code {exit_code}"
    elif not out["ok"]:
        manifest["status"] = "failed"
        manifest["error"] = "missing expected outputs"
    else:
        manifest["status"] = "success"
        manifest["error"] = None
    write_manifest(ctx.manifest_dir / f"{spec.step_id}.json", manifest)
    return manifest


def run_step(
    ctx: RunContext,
    step_id: str,
    *,
    extra_argv: Optional[List[str]] = None,
    skip_input_check: bool = False,
    stream_output: Optional[bool] = None,
) -> Dict[str, Any]:
    spec = get_step(step_id)
    repo = chemdb_repo_root()
    started = datetime.now(timezone.utc).isoformat()
    started_ts = time.perf_counter()
    inp = check_inputs(ctx, spec) if not skip_input_check else {"ok": True, "missing": [], "skipped": True}
    inputs, _ = _manifest_io_lists(ctx, spec)
    manifest: Dict[str, Any] = {
        "step_id": step_id,
        "status": "running",
        "started_at": started,
        "run_root": str(ctx.run_root.resolve()),
        "cwd": str(repo),
        "input_check": inp,
        "inputs": inputs,
        "outputs": [],
        "error": None,
    }

    if not inp.get("ok", False):
        manifest["status"] = "failed"
        manifest["error"] = "missing required inputs"
        manifest["command"] = []
        return _finalize_manifest(
            manifest, ctx=ctx, spec=spec, started_ts=started_ts, error=manifest["error"]
        )

    log_path = ctx.run_root / spec.log_file
    log_path.parent.mkdir(parents=True, exist_ok=True)

    if spec.runner == "in_process":
        manifest["command"] = [f"in_process:{step_id}"]
        with open(log_path, "w", encoding="utf-8") as log_f:
            log_f.write(f"# in-process: {step_id}\n\n")
            log_f.flush()
            try:
                result = run_in_process_step(ctx, spec)
            except Exception as exc:
                result = {"ok": False, "error": str(exc)}
            log_f.write(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
        manifest["in_process_result"] = result
        if not result.get("ok"):
            return _finalize_manifest(
                manifest,
                ctx=ctx,
                spec=spec,
                started_ts=started_ts,
                exit_code=1,
                error=result.get("error", "in-process step failed"),
            )
        return _finalize_manifest(manifest, ctx=ctx, spec=spec, started_ts=started_ts, exit_code=0)

    full_cmd = build_command(spec, ctx, extra_argv)
    manifest["command"] = full_cmd

    env = os.environ.copy()
    env = _pythonpath_env(repo, env)
    if "CHEMDB_REPO_ROOT" not in env:
        env["CHEMDB_REPO_ROOT"] = str(repo)
    env["CHEMDB_RUN_ROOT"] = str(ctx.run_root)
    stream = should_stream_step_output(stream_output)
    if stream:
        env["PYTHONUNBUFFERED"] = "1"

    with open(log_path, "w", encoding="utf-8") as log_f:
        log_f.write(f"# command: {' '.join(full_cmd)}\n")
        log_f.write(f"# cwd: {repo}\n\n")
        log_f.flush()
        if stream:
            print(f"\n>>> [{step_id}] {' '.join(full_cmd)}", flush=True)
            print(f">>> log: {log_path}\n", flush=True)
        exit_code = _run_subprocess(
            full_cmd,
            cwd=str(repo),
            env=env,
            log_f=log_f,
            stream_to_console=stream,
        )
        if stream:
            print(f"\n>>> [{step_id}] exit code {exit_code}\n", flush=True)

    manifest["log_file"] = str(log_path.resolve())
    return _finalize_manifest(
        manifest, ctx=ctx, spec=spec, started_ts=started_ts, exit_code=exit_code
    )


def run_pipeline(
    ctx: RunContext,
    through_step: str,
    *,
    pipeline: str = "core",
    extra_argv_by_step: Optional[Dict[str, List[str]]] = None,
    stream_output: Optional[bool] = None,
) -> Dict[str, Any]:
    steps = pipeline_through(through_step, pipeline=pipeline)
    results: List[Dict[str, Any]] = []
    overall_ok = True
    stream = should_stream_step_output(stream_output)
    if stream:
        print(
            f"\n=== run-pipeline ({pipeline}) through {through_step} @ {ctx.run_root} ===\n",
            flush=True,
        )
    for spec in steps:
        extra = (extra_argv_by_step or {}).get(spec.step_id)
        m = run_step(ctx, spec.step_id, extra_argv=extra, stream_output=stream_output)
        results.append({
            "step_id": spec.step_id,
            "status": m.get("status"),
            "manifest": str(ctx.manifest_dir / f"{spec.step_id}.json"),
        })
        if m.get("status") != "success":
            overall_ok = False
            break
    summary = {
        "action": "run-pipeline",
        "pipeline": pipeline,
        "through": through_step,
        "run_root": str(ctx.run_root),
        "ok": overall_ok,
        "steps": results,
        "finished_at": datetime.now(timezone.utc).isoformat(),
    }
    write_manifest(ctx.manifest_dir / f"pipeline_{pipeline}.json", summary)
    return summary


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="ChemDB Run-based orchestrator")
    sub = p.add_subparsers(dest="command", required=True)
    run_kw = {"type": str, "required": True, "help": "Path to run workspace root"}

    init_p = sub.add_parser("init-run", help="Create run dirs and copy seed data files")
    init_p.add_argument("--run-root", **run_kw)
    init_p.add_argument(
        "--pubchem-source",
        type=str,
        default=None,
        help="Optional CSV file or directory of CSVs to copy into run/data/pubchem",
    )

    stream_kw = {
        "action": "store_true",
        "help": "Stream step stdout/stderr to console (also set CHEMDB_STREAM_STEP_OUTPUT=1)",
    }
    all_steps = sorted(STEP_REGISTRY.keys())
    step_p = sub.add_parser("run-step", help="Run a single pipeline step")
    step_p.add_argument("--run-root", **run_kw)
    step_p.add_argument("--step", type=str, required=True, choices=all_steps)
    step_p.add_argument("--stream", **stream_kw)

    pipe_p = sub.add_parser("run-pipeline", help="Run steps through a target")
    pipe_p.add_argument("--run-root", **run_kw)
    pipe_p.add_argument(
        "--pipeline",
        type=str,
        default="core",
        choices=("core", "training"),
        help="Pipeline to run (default: core)",
    )
    pipe_p.add_argument(
        "--through",
        type=str,
        required=True,
        help="Last step id to run (inclusive)",
    )
    pipe_p.add_argument("--stream", **stream_kw)

    return p


def _validate_through(through: str, pipeline: str) -> None:
    if pipeline == "core" and through == "step4_5":
        raise SystemExit(
            "core pipeline no longer runs step4_5 (auto GA generation). "
            "Use: python -m app.services.ga_registry_manager generate-ga-from-run | "
            "bind-ga-version-to-run | apply-ga-to-run, then "
            "run-pipeline --through apply_ga_to_run or step6_7."
        )
    order = TRAINING_PIPELINE_ORDER if pipeline == "training" else PIPELINE_ORDER
    if through not in order:
        raise SystemExit(f"--through must be one of {order}, got {through!r}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    ctx = RunContext.from_run_root(args.run_root)
    pipeline = getattr(args, "pipeline", "core")

    if args.command == "init-run":
        pubchem = Path(args.pubchem_source).resolve() if getattr(args, "pubchem_source", None) else None
        try:
            init_run(ctx, pubchem_source=pubchem)
        except FileNotFoundError as exc:
            print(str(exc), file=sys.stderr)
            return 1
        print(f"Initialized run at {ctx.run_root}")
        return 0

    if args.command == "run-step":
        m = run_step(ctx, args.step, stream_output=getattr(args, "stream", False) or None)
        print(json.dumps(m, indent=2, ensure_ascii=False))
        return 0 if m.get("status") == "success" else 1

    if args.command == "run-pipeline":
        _validate_through(args.through, pipeline)
        s = run_pipeline(
            ctx,
            args.through,
            pipeline=pipeline,
            stream_output=getattr(args, "stream", False) or None,
        )
        print(json.dumps(s, indent=2, ensure_ascii=False))
        return 0 if s.get("ok") else 1

    parser.error("Unknown command")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
