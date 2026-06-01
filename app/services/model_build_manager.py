"""Model build orchestration through background jobs (Phase 3B MVP)."""

from __future__ import annotations

import shlex
import sys
from pathlib import Path
from typing import Any, Dict, Optional

from app.services import job_manager


def _quote(parts: list[str]) -> str:
    return " ".join(shlex.quote(p) for p in parts)


def submit_model_build_job(
    *,
    workspace_root: str | Path,
    db_path: str | Path,
    project_id: str,
    run_id: str,
    ga_set_id: str,
    ga_version_id: str,
    pubchem_source: Optional[str] = None,
    hp_set_id: Optional[str] = None,
    hp_version_id: Optional[str] = None,
) -> Dict[str, Any]:
    ws = Path(workspace_root).resolve()
    run_root = ws / "projects" / project_id / "runs" / run_id

    init_cmd = [sys.executable, "-m", "app.services.orchestrator", "init-run", "--run-root", str(run_root)]
    if pubchem_source:
        init_cmd.extend(["--pubchem-source", pubchem_source])
    bind_cmd = [
        sys.executable,
        "-m",
        "app.services.ga_registry_manager",
        "bind-ga-version-to-run",
        "--workspace-root",
        str(ws),
        "--run-root",
        str(run_root),
        "--ga-set-id",
        ga_set_id,
        "--ga-version-id",
        ga_version_id,
    ]
    core_cmd = [
        sys.executable,
        "-m",
        "app.services.orchestrator",
        "run-pipeline",
        "--run-root",
        str(run_root),
        "--pipeline",
        "core",
        "--through",
        "step13",
    ]
    train_steps = [
        "build_l3_embedding_ecfp",
        "zscore_metal_embedding",
        "prepare_training_index",
        "prepare_training_config",
    ]
    train_tail = ["training_data", "training_train"]
    train_cmds = [
        [
            sys.executable,
            "-m",
            "app.services.orchestrator",
            "run-step",
            "--run-root",
            str(run_root),
            "--step",
            step,
        ]
        for step in (train_steps + train_tail)
    ]
    hp_apply_cmd = None
    if hp_set_id and hp_version_id:
        hp_apply_cmd = [
            sys.executable,
            "-m",
            "app.services.hyperparameter_manager",
            "apply-to-run",
            "--workspace-root",
            str(ws),
            "--run-root",
            str(run_root),
            "--hp-set-id",
            hp_set_id,
            "--hp-version-id",
            hp_version_id,
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
    chain_cmds = [init_cmd, bind_cmd, core_cmd] + train_cmds[:4]
    if hp_apply_cmd:
        chain_cmds.append(hp_apply_cmd)
    chain_cmds.extend(train_cmds[4:])
    chain_cmds.append(rebuild_cmd)
    full = " && ".join(_quote(c) for c in chain_cmds)
    command = ["bash", "-lc", full]
    return job_manager.launch_subprocess_job_async(
        db_path,
        job_type="model_build",
        title=f"Model build {project_id}/{run_id}",
        command=command,
        workspace_root=ws,
        project_id=project_id,
        run_id=run_id,
        cwd=ws,
    )
