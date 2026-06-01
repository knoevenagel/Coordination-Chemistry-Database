"""Filesystem → SQLite ingest (Phase 2A read-only registry)."""

from __future__ import annotations

import csv
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None  # type: ignore

STEP_MANIFEST_SKIP = frozenset({
    "run.json",
    "pipeline_core.json",
    "pipeline_training.json",
    "ga_binding.json",
    "ga_stale_report.json",
})


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _read_json(path: Path) -> Optional[Dict[str, Any]]:
    if not path.is_file():
        return None
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        return data if isinstance(data, dict) else None
    except (json.JSONDecodeError, OSError):
        return None


def _read_json_list(path: Path) -> Optional[List[Any]]:
    if not path.is_file():
        return None
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        return data if isinstance(data, list) else None
    except (json.JSONDecodeError, OSError):
        return None


def _count_csv_rows(path: Path) -> int:
    if not path.is_file():
        return 0
    try:
        with open(path, "r", encoding="utf-8") as f:
            reader = csv.reader(f)
            next(reader, None)
            return sum(1 for _ in reader)
    except OSError:
        return 0


def workspace_from_path(path: Path) -> Optional[Path]:
    parts = path.resolve().parts
    for anchor in ("projects", "ga_sets", "model_after_tasks", "model_after_results"):
        if anchor in parts:
            idx = parts.index(anchor)
            return Path(*parts[:idx])
    return None


def parse_project_run_from_run_root(run_root: Path) -> Tuple[Optional[str], Optional[str]]:
    parts = run_root.resolve().parts
    if "projects" in parts and "runs" in parts:
        pi = parts.index("projects")
        ri = parts.index("runs")
        if pi + 1 < len(parts) and ri + 1 < len(parts):
            return parts[pi + 1], parts[ri + 1]
    return None, run_root.name


def parse_project_run_from_model_run_root(model_run_root: Path) -> Tuple[Optional[str], Optional[str]]:
    return parse_project_run_from_run_root(model_run_root)


def artifact_type_from_path(path: str) -> str:
    suffix = Path(path).suffix.lower()
    mapping = {
        ".csv": "csv",
        ".json": "json",
        ".yaml": "yaml",
        ".yml": "yaml",
        ".npz": "npz",
        ".pt": "pt",
        ".pkl": "pkl",
        ".log": "log",
    }
    return mapping.get(suffix, "other")


def synthetic_model_id(project_id: str, run_id: str, checkpoint_path: Path) -> str:
    stem = checkpoint_path.stem
    return f"{project_id}_{run_id}_{stem}"


def _file_size(path: Path) -> Optional[int]:
    try:
        return path.stat().st_size if path.is_file() else None
    except OSError:
        return None


def _rel_to_run(run_root: Path, abs_path: str) -> str:
    try:
        return str(Path(abs_path).resolve().relative_to(run_root.resolve()))
    except ValueError:
        return abs_path


def ingest_project(
    conn: sqlite3.Connection,
    workspace_root: Path,
    project_id: str,
) -> Dict[str, Any]:
    workspace_root = workspace_root.resolve()
    project_root = workspace_root / "projects" / project_id
    result: Dict[str, Any] = {"project_id": project_id, "ingest_status": "ok", "ingest_error": None}

    if not project_root.is_dir():
        result["ingest_status"] = "missing"
        result["ingest_error"] = f"project directory not found: {project_root}"
        conn.execute(
            """
            INSERT INTO projects (project_id, project_root, name, ingest_status, ingest_error, ingested_at)
            VALUES (?, ?, ?, ?, ?, ?)
            ON CONFLICT(project_id) DO UPDATE SET
                project_root=excluded.project_root,
                name=excluded.name,
                ingest_status=excluded.ingest_status,
                ingest_error=excluded.ingest_error,
                ingested_at=excluded.ingested_at
            """,
            (project_id, str(project_root), project_id, result["ingest_status"], result["ingest_error"], _utc_now()),
        )
        conn.commit()
        return result

    name = project_id
    created_at = None
    updated_at = None
    notes = None
    project_json = _read_json(project_root / "project.json")
    if project_json:
        name = str(project_json.get("name") or project_id)
        created_at = project_json.get("created_at")
        updated_at = project_json.get("updated_at")
        notes = project_json.get("notes")

    conn.execute(
        """
        INSERT INTO projects (project_id, project_root, name, created_at, updated_at, notes,
                              ingest_status, ingest_error, ingested_at)
        VALUES (?, ?, ?, ?, ?, ?, 'ok', NULL, ?)
        ON CONFLICT(project_id) DO UPDATE SET
            project_root=excluded.project_root,
            name=excluded.name,
            created_at=excluded.created_at,
            updated_at=excluded.updated_at,
            notes=excluded.notes,
            ingest_status=excluded.ingest_status,
            ingest_error=excluded.ingest_error,
            ingested_at=excluded.ingested_at
        """,
        (project_id, str(project_root), name, created_at, updated_at, notes, _utc_now()),
    )
    conn.commit()
    return result


def ingest_run(conn: sqlite3.Connection, run_root: Path) -> Dict[str, Any]:
    run_root = run_root.resolve()
    project_id, run_id = parse_project_run_from_run_root(run_root)
    result: Dict[str, Any] = {
        "project_id": project_id,
        "run_id": run_id,
        "ingest_status": "ok",
        "ingest_error": None,
    }

    if not project_id or not run_id:
        result["ingest_status"] = "error"
        result["ingest_error"] = f"cannot parse project_id/run_id from {run_root}"
        return result

    ws = workspace_from_path(run_root)
    if ws is not None:
        ingest_project(conn, ws, project_id)

    run_json = _read_json(run_root / "manifests" / "run.json") or {}
    pipeline_core = _read_json(run_root / "manifests" / "pipeline_core.json") or {}
    pipeline_training = _read_json(run_root / "manifests" / "pipeline_training.json") or {}

    pubchem_files = run_json.get("pubchem_files") or []
    pubchem_count = len(pubchem_files) if isinstance(pubchem_files, list) else None

    last_finished = pipeline_training.get("finished_at") or pipeline_core.get("finished_at")

    ga_binding_status = "missing"
    if (run_root / "manifests" / "ga_binding.json").is_file():
        ga_binding_status = "bound"
    elif (run_root / "manifests" / "step4_5.json").is_file():
        ga_binding_status = "legacy_auto_ga"

    conn.execute(
        """
        INSERT INTO runs (
            project_id, run_id, run_root, created_at, chemdb_repo_root, status,
            pubchem_file_count, pipeline_core_ok, pipeline_training_ok, last_finished_at,
            ga_binding_status, ingest_status, ingest_error, ingested_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 'ok', NULL, ?)
        ON CONFLICT(project_id, run_id) DO UPDATE SET
            run_root=excluded.run_root,
            created_at=excluded.created_at,
            chemdb_repo_root=excluded.chemdb_repo_root,
            status=excluded.status,
            pubchem_file_count=excluded.pubchem_file_count,
            pipeline_core_ok=excluded.pipeline_core_ok,
            pipeline_training_ok=excluded.pipeline_training_ok,
            last_finished_at=excluded.last_finished_at,
            ga_binding_status=excluded.ga_binding_status,
            ingest_status=excluded.ingest_status,
            ingest_error=excluded.ingest_error,
            ingested_at=excluded.ingested_at
        """,
        (
            project_id,
            run_id,
            str(run_root),
            run_json.get("created_at"),
            run_json.get("chemdb_repo_root"),
            run_json.get("status"),
            pubchem_count,
            1 if pipeline_core.get("ok") else (0 if "ok" in pipeline_core else None),
            1 if pipeline_training.get("ok") else (0 if "ok" in pipeline_training else None),
            last_finished,
            ga_binding_status,
            _utc_now(),
        ),
    )
    conn.commit()
    return result


def ingest_ga_set(conn: sqlite3.Connection, workspace_root: Path, ga_set_id: str) -> Dict[str, Any]:
    workspace_root = workspace_root.resolve()
    ga_set_root = workspace_root / "ga_sets" / ga_set_id
    result: Dict[str, Any] = {"ga_set_id": ga_set_id, "versions": [], "ingest_status": "ok"}

    meta_path = ga_set_root / "ga_set.json"
    meta = _read_json(meta_path) or {}
    if not ga_set_root.is_dir():
        result["ingest_status"] = "missing"
        result["ingest_error"] = f"ga set not found: {ga_set_root}"
        conn.execute(
            """
            INSERT INTO ga_sets (ga_set_id, ga_set_root, name, ingest_status, ingest_error, ingested_at)
            VALUES (?, ?, ?, ?, ?, ?)
            ON CONFLICT(ga_set_id) DO UPDATE SET
                ga_set_root=excluded.ga_set_root,
                name=excluded.name,
                ingest_status=excluded.ingest_status,
                ingest_error=excluded.ingest_error,
                ingested_at=excluded.ingested_at
            """,
            (ga_set_id, str(ga_set_root), ga_set_id, result["ingest_status"], result["ingest_error"], _utc_now()),
        )
        conn.commit()
        return result

    tags = meta.get("tags")
    source = meta.get("source")
    conn.execute(
        """
        INSERT INTO ga_sets (
            ga_set_id, ga_set_root, name, description, created_at, updated_at,
            tags_json, source_json, ingest_status, ingest_error, ingested_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, 'ok', NULL, ?)
        ON CONFLICT(ga_set_id) DO UPDATE SET
            ga_set_root=excluded.ga_set_root,
            name=excluded.name,
            description=excluded.description,
            created_at=excluded.created_at,
            updated_at=excluded.updated_at,
            tags_json=excluded.tags_json,
            source_json=excluded.source_json,
            ingest_status=excluded.ingest_status,
            ingest_error=excluded.ingest_error,
            ingested_at=excluded.ingested_at
        """,
        (
            ga_set_id,
            str(ga_set_root),
            meta.get("name") or ga_set_id,
            meta.get("description"),
            meta.get("created_at"),
            meta.get("updated_at"),
            json.dumps(tags, ensure_ascii=False) if tags is not None else None,
            json.dumps(source, ensure_ascii=False) if source is not None else None,
            _utc_now(),
        ),
    )

    versions_root = ga_set_root / "versions"
    if versions_root.is_dir():
        for vdir in sorted(p for p in versions_root.iterdir() if p.is_dir()):
            vid = vdir.name
            vmeta = _read_json(vdir / "ga_version.json") or {}
            csv_path = vdir / "GA_with_id.csv"
            conn.execute(
                """
                INSERT INTO ga_versions (
                    ga_set_id, ga_version_id, parent_version_id, created_at, created_by,
                    num_ga, checksum, ga_csv_path, source_json, notes,
                    ingest_status, ingest_error, ingested_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 'ok', NULL, ?)
                ON CONFLICT(ga_set_id, ga_version_id) DO UPDATE SET
                    parent_version_id=excluded.parent_version_id,
                    created_at=excluded.created_at,
                    created_by=excluded.created_by,
                    num_ga=excluded.num_ga,
                    checksum=excluded.checksum,
                    ga_csv_path=excluded.ga_csv_path,
                    source_json=excluded.source_json,
                    notes=excluded.notes,
                    ingest_status=excluded.ingest_status,
                    ingest_error=excluded.ingest_error,
                    ingested_at=excluded.ingested_at
                """,
                (
                    ga_set_id,
                    vid,
                    vmeta.get("parent_version_id"),
                    vmeta.get("created_at"),
                    vmeta.get("created_by"),
                    vmeta.get("num_ga"),
                    vmeta.get("checksum"),
                    str(csv_path.resolve()) if csv_path.is_file() else None,
                    json.dumps(vmeta.get("source"), ensure_ascii=False) if vmeta.get("source") else None,
                    vmeta.get("notes"),
                    _utc_now(),
                ),
            )
            result["versions"].append(vid)

    conn.commit()
    return result


def ingest_run_ga_binding(conn: sqlite3.Connection, run_root: Path) -> Dict[str, Any]:
    run_root = run_root.resolve()
    project_id, run_id = parse_project_run_from_run_root(run_root)
    if not project_id or not run_id:
        return {"ingest_status": "error", "ingest_error": "invalid run_root path"}

    binding_path = run_root / "manifests" / "ga_binding.json"
    binding = _read_json(binding_path)

    if binding is None:
        legacy = (run_root / "manifests" / "step4_5.json").is_file()
        status = "legacy_auto_ga" if legacy else "missing"
        conn.execute(
            """
            INSERT INTO run_ga_bindings (
                project_id, run_id, ga_set_id, ga_version_id, checksum, num_ga,
                source_ga_csv, materialized_path, bound_at, bound_by, status,
                is_current, ingest_status, ingest_error, ingested_at
            ) VALUES (?, ?, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, ?, 1, 'ok', NULL, ?)
            ON CONFLICT(project_id, run_id) DO UPDATE SET
                ga_set_id=excluded.ga_set_id,
                ga_version_id=excluded.ga_version_id,
                checksum=excluded.checksum,
                num_ga=excluded.num_ga,
                source_ga_csv=excluded.source_ga_csv,
                materialized_path=excluded.materialized_path,
                bound_at=excluded.bound_at,
                bound_by=excluded.bound_by,
                status=excluded.status,
                is_current=excluded.is_current,
                ingest_status=excluded.ingest_status,
                ingest_error=excluded.ingest_error,
                ingested_at=excluded.ingested_at
            """,
            (project_id, run_id, status, _utc_now()),
        )
        conn.commit()
        return {"status": status, "ingest_status": "ok"}

    ga_set_id = binding.get("ga_set_id")
    if ga_set_id:
        ws = workspace_from_path(run_root)
        if ws is not None:
            ga_set_root = ws / "ga_sets" / str(ga_set_id)
            if ga_set_root.is_dir():
                ingest_ga_set(conn, ws, str(ga_set_id))
            else:
                conn.execute(
                    """
                    INSERT INTO ga_sets (ga_set_id, ga_set_root, name, ingest_status, ingested_at)
                    VALUES (?, ?, ?, 'stub', ?)
                    ON CONFLICT(ga_set_id) DO NOTHING
                    """,
                    (str(ga_set_id), str(ga_set_root), str(ga_set_id), _utc_now()),
                )

    conn.execute(
        """
        INSERT INTO run_ga_bindings (
            project_id, run_id, ga_set_id, ga_version_id, checksum, num_ga,
            source_ga_csv, materialized_path, bound_at, bound_by, status,
            is_current, ingest_status, ingest_error, ingested_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 1, 'ok', NULL, ?)
        ON CONFLICT(project_id, run_id) DO UPDATE SET
            ga_set_id=excluded.ga_set_id,
            ga_version_id=excluded.ga_version_id,
            checksum=excluded.checksum,
            num_ga=excluded.num_ga,
            source_ga_csv=excluded.source_ga_csv,
            materialized_path=excluded.materialized_path,
            bound_at=excluded.bound_at,
            bound_by=excluded.bound_by,
            status=excluded.status,
            is_current=excluded.is_current,
            ingest_status=excluded.ingest_status,
            ingest_error=excluded.ingest_error,
            ingested_at=excluded.ingested_at
        """,
        (
            project_id,
            run_id,
            binding.get("ga_set_id"),
            binding.get("ga_version_id"),
            binding.get("checksum"),
            binding.get("num_ga"),
            binding.get("source_ga_csv"),
            binding.get("materialized_path"),
            binding.get("bound_at"),
            binding.get("bound_by"),
            binding.get("status") or "bound",
            _utc_now(),
        ),
    )
    conn.commit()
    return {"status": "bound", "ga_set_id": binding.get("ga_set_id"), "ingest_status": "ok"}


def _upsert_artifact(
    conn: sqlite3.Connection,
    *,
    project_id: str,
    run_id: str,
    step_id: str,
    role: str,
    abs_path: str,
    run_root: Path,
) -> None:
    path = Path(abs_path)
    rel = _rel_to_run(run_root, abs_path)
    exists = path.is_file() or path.is_dir()
    size = _file_size(path) if path.is_file() else None
    conn.execute(
        """
        INSERT INTO artifacts (
            project_id, run_id, step_id, rel_path, role, abs_path, size_bytes,
            exists_flag, artifact_type, ingest_status, ingested_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, 'ok', ?)
        ON CONFLICT(project_id, run_id, step_id, rel_path, role) DO UPDATE SET
            abs_path=excluded.abs_path,
            size_bytes=excluded.size_bytes,
            exists_flag=excluded.exists_flag,
            artifact_type=excluded.artifact_type,
            ingested_at=excluded.ingested_at
        """,
        (
            project_id,
            run_id,
            step_id,
            rel,
            role,
            str(path.resolve()) if exists else abs_path,
            size,
            1 if exists else 0,
            artifact_type_from_path(abs_path),
            _utc_now(),
        ),
    )


def ingest_step_executions(conn: sqlite3.Connection, run_root: Path) -> Dict[str, Any]:
    run_root = run_root.resolve()
    project_id, run_id = parse_project_run_from_run_root(run_root)
    if not project_id or not run_id:
        return {"ingest_status": "error", "ingest_error": "invalid run_root", "count": 0}

    manifest_dir = run_root / "manifests"
    count = 0
    if not manifest_dir.is_dir():
        conn.commit()
        return {"count": 0, "ingest_status": "ok"}

    for mf in sorted(manifest_dir.glob("*.json")):
        if mf.name in STEP_MANIFEST_SKIP:
            continue
        data = _read_json(mf)
        if not data or "step_id" not in data:
            continue
        step_id = str(data["step_id"])
        started_at = str(data.get("started_at") or "")
        command = data.get("command")
        in_proc = data.get("in_process_result")
        conn.execute(
            """
            INSERT INTO step_executions (
                project_id, run_id, step_id, started_at, status, command_json, cwd,
                finished_at, duration_seconds, exit_code, error, log_path, manifest_path,
                in_process_result_json, ingest_status, ingested_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 'ok', ?)
            ON CONFLICT(project_id, run_id, step_id, started_at) DO UPDATE SET
                status=excluded.status,
                command_json=excluded.command_json,
                cwd=excluded.cwd,
                finished_at=excluded.finished_at,
                duration_seconds=excluded.duration_seconds,
                exit_code=excluded.exit_code,
                error=excluded.error,
                log_path=excluded.log_path,
                manifest_path=excluded.manifest_path,
                in_process_result_json=excluded.in_process_result_json,
                ingested_at=excluded.ingested_at
            """,
            (
                project_id,
                run_id,
                step_id,
                started_at,
                data.get("status"),
                json.dumps(command, ensure_ascii=False) if command is not None else None,
                data.get("cwd"),
                data.get("finished_at"),
                data.get("duration_seconds"),
                data.get("exit_code"),
                data.get("error"),
                data.get("log_file"),
                str(mf.resolve()),
                json.dumps(in_proc, ensure_ascii=False) if in_proc is not None else None,
                _utc_now(),
            ),
        )
        count += 1

        for role_key, role in (("inputs", "input"), ("outputs", "output")):
            for entry in data.get(role_key) or []:
                if isinstance(entry, dict) and entry.get("path"):
                    _upsert_artifact(
                        conn,
                        project_id=project_id,
                        run_id=run_id,
                        step_id=step_id,
                        role=role,
                        abs_path=str(entry["path"]),
                        run_root=run_root,
                    )

    conn.commit()
    return {"count": count, "ingest_status": "ok"}


def ingest_key_artifacts(conn: sqlite3.Connection, run_root: Path) -> Dict[str, Any]:
    run_root = run_root.resolve()
    project_id, run_id = parse_project_run_from_run_root(run_root)
    if not project_id or not run_id:
        return {"count": 0, "ingest_status": "error"}

    key_paths = [
        "training/index.json",
        "training/config.yaml",
        "training/ckpts/best.pt",
        "training/ckpts/last.pt",
        "training/ckpts/history.json",
        "tmp/GA_with_id.csv",
        "tmp/step13_kl_nl_samples.csv",
        "data/L3_embedding/L3_embedding_ecfp.npz",
    ]
    count = 0
    for rel in key_paths:
        path = run_root / rel
        if path.is_file():
            _upsert_artifact(
                conn,
                project_id=project_id,
                run_id=run_id,
                step_id="",
                role="key",
                abs_path=str(path),
                run_root=run_root,
            )
            count += 1
    conn.commit()
    return {"count": count, "ingest_status": "ok"}


def _load_yaml(path: Path) -> Optional[Dict[str, Any]]:
    if not path.is_file() or yaml is None:
        return None
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
        return data if isinstance(data, dict) else None
    except (OSError, yaml.YAMLError):
        return None


def _best_history_metrics(history_path: Path) -> Tuple[Optional[int], Optional[float], Optional[float], Optional[float]]:
    rows = _read_json_list(history_path)
    if not rows:
        return None, None, None, None
    best = rows[-1]
    if not isinstance(best, dict):
        return None, None, None, None
    return (
        best.get("epoch"),
        best.get("mrr"),
        best.get("recall_at_1"),
        best.get("recall_at_5"),
    )


def ingest_model(conn: sqlite3.Connection, run_root: Path) -> Dict[str, Any]:
    run_root = run_root.resolve()
    project_id, run_id = parse_project_run_from_run_root(run_root)
    if not project_id or not run_id:
        return {"count": 0, "ingest_status": "error"}

    index_path = run_root / "training" / "index.json"
    config_path = run_root / "training" / "config.yaml"
    history_path = run_root / "training" / "ckpts" / "history.json"
    ckpt_dir = run_root / "training" / "ckpts"

    index = _read_json(index_path) or {}
    config = _load_yaml(config_path) or {}
    embedding = None
    if config.get("data") and isinstance(config["data"], dict):
        embedding = config["data"].get("embedding")
    lig = index.get("ligand_embedding") or {}
    if not embedding:
        embedding = lig.get("backend")
    l3_path = lig.get("path")
    metal = index.get("metal_embedding") or {}
    metal_path = metal.get("path")

    binding_row = conn.execute(
        "SELECT ga_set_id, ga_version_id, checksum FROM run_ga_bindings WHERE project_id=? AND run_id=?",
        (project_id, run_id),
    ).fetchone()
    ga_set_id = binding_row["ga_set_id"] if binding_row else None
    ga_version_id = binding_row["ga_version_id"] if binding_row else None
    ga_checksum = binding_row["checksum"] if binding_row else None

    best_epoch, val_mrr, val_r1, val_r5 = _best_history_metrics(history_path)

    count = 0
    if ckpt_dir.is_dir():
        for ckpt in sorted(ckpt_dir.glob("*.pt")):
            mid = synthetic_model_id(project_id, run_id, ckpt)
            conn.execute(
                """
                INSERT INTO models (
                    model_id, project_id, run_id, checkpoint_path, checkpoint_stem,
                    config_path, index_path, embedding_backend, l3_embedding_path,
                    metal_embedding_path, history_path, best_epoch, val_mrr,
                    val_recall_at_1, val_recall_at_5, ga_set_id, ga_version_id,
                    ga_binding_checksum, size_bytes, ingest_status, ingested_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 'ok', ?)
                ON CONFLICT(model_id) DO UPDATE SET
                    checkpoint_path=excluded.checkpoint_path,
                    checkpoint_stem=excluded.checkpoint_stem,
                    config_path=excluded.config_path,
                    index_path=excluded.index_path,
                    embedding_backend=excluded.embedding_backend,
                    l3_embedding_path=excluded.l3_embedding_path,
                    metal_embedding_path=excluded.metal_embedding_path,
                    history_path=excluded.history_path,
                    best_epoch=excluded.best_epoch,
                    val_mrr=excluded.val_mrr,
                    val_recall_at_1=excluded.val_recall_at_1,
                    val_recall_at_5=excluded.val_recall_at_5,
                    ga_set_id=excluded.ga_set_id,
                    ga_version_id=excluded.ga_version_id,
                    ga_binding_checksum=excluded.ga_binding_checksum,
                    size_bytes=excluded.size_bytes,
                    ingested_at=excluded.ingested_at
                """,
                (
                    mid,
                    project_id,
                    run_id,
                    str(ckpt.resolve()),
                    ckpt.stem,
                    str(config_path.resolve()) if config_path.is_file() else None,
                    str(index_path.resolve()) if index_path.is_file() else None,
                    embedding,
                    l3_path,
                    metal_path,
                    str(history_path.resolve()) if history_path.is_file() else None,
                    best_epoch,
                    val_mrr,
                    val_r1,
                    val_r5,
                    ga_set_id,
                    ga_version_id,
                    ga_checksum,
                    _file_size(ckpt),
                    _utc_now(),
                ),
            )
            count += 1

    conn.commit()
    return {"count": count, "ingest_status": "ok"}


def ingest_model_after_task(conn: sqlite3.Connection, task_dir: Path) -> Dict[str, Any]:
    task_dir = task_dir.resolve()
    task_json = _read_json(task_dir / "task.json") or {}
    task_id = str(task_json.get("task_id") or task_dir.name)

    cand_file = task_json.get("candidate_file") or "candidates.csv"
    pos_file = task_json.get("positive_file") or "positives.csv"
    neg_file = task_json.get("negative_file") or "negatives.csv"

    cand_path = task_dir / cand_file
    pos_path = task_dir / pos_file
    neg_path = task_dir / neg_file if neg_file else None

    conn.execute(
        """
        INSERT INTO model_after_tasks (
            task_id, task_dir, metal, embedding, candidate_file, positive_file,
            negative_file, notes, candidate_count, positive_count, negative_count,
            ingest_status, ingested_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 'ok', ?)
        ON CONFLICT(task_id) DO UPDATE SET
            task_dir=excluded.task_dir,
            metal=excluded.metal,
            embedding=excluded.embedding,
            candidate_file=excluded.candidate_file,
            positive_file=excluded.positive_file,
            negative_file=excluded.negative_file,
            notes=excluded.notes,
            candidate_count=excluded.candidate_count,
            positive_count=excluded.positive_count,
            negative_count=excluded.negative_count,
            ingested_at=excluded.ingested_at
        """,
        (
            task_id,
            str(task_dir),
            task_json.get("metal"),
            task_json.get("embedding"),
            cand_file,
            pos_file,
            neg_file or None,
            task_json.get("notes"),
            _count_csv_rows(cand_path),
            _count_csv_rows(pos_path),
            _count_csv_rows(neg_path) if neg_path and neg_path.is_file() else 0,
            _utc_now(),
        ),
    )
    conn.commit()
    return {"task_id": task_id, "ingest_status": "ok"}


def _parse_batch_task_id(batch_dir: Path) -> str:
    parts = batch_dir.resolve().parts
    if "model_after_results" in parts:
        idx = parts.index("model_after_results")
        if idx + 1 < len(parts):
            return parts[idx + 1]
    return batch_dir.parent.name


def _resolve_registry_model_id(
    conn: sqlite3.Connection,
    project_id: Optional[str],
    run_id: Optional[str],
    checkpoint_path: Optional[str],
) -> Optional[str]:
    if not project_id or not run_id or not checkpoint_path:
        return None
    row = conn.execute(
        """
        SELECT model_id FROM models
        WHERE project_id=? AND run_id=? AND checkpoint_path=?
        """,
        (project_id, run_id, checkpoint_path),
    ).fetchone()
    if row:
        return row["model_id"]
    ckpt = Path(checkpoint_path)
    if ckpt.is_file() or ckpt.suffix == ".pt":
        return synthetic_model_id(project_id, run_id, ckpt)
    return None


def _ingest_model_result_row(
    conn: sqlite3.Connection,
    *,
    task_id: str,
    batch_id: str,
    model_id: str,
    metrics: Dict[str, Any],
    summary: Optional[Dict[str, Any]] = None,
    ranking_path: Optional[str] = None,
    metrics_path: Optional[str] = None,
    report_path: Optional[str] = None,
) -> None:
    model_run_root = metrics.get("model_run_root") or (summary or {}).get("model_run_root")
    checkpoint_path = metrics.get("checkpoint_path") or (summary or {}).get("checkpoint_path")
    proj, rid = parse_project_run_from_model_run_root(Path(model_run_root)) if model_run_root else (None, None)
    registry_mid = _resolve_registry_model_id(conn, proj, rid, checkpoint_path)

    conn.execute(
        """
        INSERT INTO model_after_model_results (
            task_id, batch_id, model_id, model_run_root, checkpoint_path,
            project_id, run_id, registry_model_id, model_name,
            mrr, hit_at_5, hit_at_10, hit_at_20, hit_at_50,
            selection_score, rank_among_models, status,
            ranking_path, metrics_path, report_path, ingest_status, ingested_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 'ok', ?)
        ON CONFLICT(task_id, batch_id, model_id) DO UPDATE SET
            model_run_root=excluded.model_run_root,
            checkpoint_path=excluded.checkpoint_path,
            project_id=excluded.project_id,
            run_id=excluded.run_id,
            registry_model_id=excluded.registry_model_id,
            model_name=excluded.model_name,
            mrr=excluded.mrr,
            hit_at_5=excluded.hit_at_5,
            hit_at_10=excluded.hit_at_10,
            hit_at_20=excluded.hit_at_20,
            hit_at_50=excluded.hit_at_50,
            selection_score=excluded.selection_score,
            rank_among_models=excluded.rank_among_models,
            status=excluded.status,
            ranking_path=excluded.ranking_path,
            metrics_path=excluded.metrics_path,
            report_path=excluded.report_path,
            ingested_at=excluded.ingested_at
        """,
        (
            task_id,
            batch_id,
            model_id,
            model_run_root,
            checkpoint_path,
            proj,
            rid,
            registry_mid,
            (summary or {}).get("model_name"),
            metrics.get("mrr") if metrics.get("mrr") is not None else (summary or {}).get("mrr"),
            metrics.get("hit_at_5") if metrics.get("hit_at_5") is not None else (summary or {}).get("hit_at_5"),
            metrics.get("hit_at_10") if metrics.get("hit_at_10") is not None else (summary or {}).get("hit_at_10"),
            metrics.get("hit_at_20") if metrics.get("hit_at_20") is not None else (summary or {}).get("hit_at_20"),
            metrics.get("hit_at_50") if metrics.get("hit_at_50") is not None else (summary or {}).get("hit_at_50"),
            (summary or {}).get("selection_score"),
            (summary or {}).get("rank_among_models"),
            metrics.get("status") or (summary or {}).get("status"),
            ranking_path or metrics.get("ranking_path"),
            metrics_path,
            report_path,
            _utc_now(),
        ),
    )


def ingest_model_after_batch(conn: sqlite3.Connection, batch_dir: Path) -> Dict[str, Any]:
    batch_dir = batch_dir.resolve()
    task_id = _parse_batch_task_id(batch_dir)
    batch_id = batch_dir.name

    ws = workspace_from_path(batch_dir)
    task_root = (ws / "model_after_tasks" / task_id) if ws else batch_dir.parent / task_id
    if task_root.is_dir():
        ingest_model_after_task(conn, task_root)
    else:
        conn.execute(
            """
            INSERT INTO model_after_tasks (task_id, task_dir, ingest_status, ingested_at)
            VALUES (?, ?, 'stub', ?)
            ON CONFLICT(task_id) DO NOTHING
            """,
            (task_id, str(task_root), _utc_now()),
        )

    eval_models = _read_json(batch_dir / "evaluate_models_manifest.json") or {}
    eval_model = _read_json(batch_dir / "evaluate_model_manifest.json") or {}
    selection = _read_json(batch_dir / "model_selection_manifest.json") or {}
    best = _read_json(batch_dir / "best_model.json") or {}

    command = eval_models.get("command") or eval_model.get("command")
    model_runs_file = eval_models.get("model_runs_file")
    finished_at = eval_models.get("finished_at") or eval_model.get("finished_at")
    n_models = selection.get("n_models")
    n_success = selection.get("n_success_with_score")
    best_model_id = best.get("model_id")
    summary_path = selection.get("summary_path") or str(batch_dir / "model_selection_summary.csv")

    conn.execute(
        """
        INSERT INTO model_after_batches (
            task_id, batch_id, output_dir, command, model_runs_file, finished_at,
            n_models, n_success, best_model_id, summary_path, ingest_status, ingested_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 'ok', ?)
        ON CONFLICT(task_id, batch_id) DO UPDATE SET
            output_dir=excluded.output_dir,
            command=excluded.command,
            model_runs_file=excluded.model_runs_file,
            finished_at=excluded.finished_at,
            n_models=excluded.n_models,
            n_success=excluded.n_success,
            best_model_id=excluded.best_model_id,
            summary_path=excluded.summary_path,
            ingested_at=excluded.ingested_at
        """,
        (
            task_id,
            batch_id,
            str(batch_dir),
            command,
            model_runs_file,
            finished_at,
            n_models,
            n_success,
            best_model_id,
            summary_path if Path(summary_path).is_file() else None,
            _utc_now(),
        ),
    )

    summary_by_model: Dict[str, Dict[str, Any]] = {}
    summary_csv = batch_dir / "model_selection_summary.csv"
    if summary_csv.is_file():
        with open(summary_csv, "r", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                mid = row.get("model_id")
                if mid:
                    summary_by_model[mid] = row

    result_count = 0
    models_dir = batch_dir / "models"
    if models_dir.is_dir():
        for model_dir in sorted(p for p in models_dir.iterdir() if p.is_dir()):
            mid = model_dir.name
            metrics = _read_json(model_dir / "metrics.json") or {}
            _ingest_model_result_row(
                conn,
                task_id=task_id,
                batch_id=batch_id,
                model_id=mid,
                metrics=metrics,
                summary=summary_by_model.get(mid),
                ranking_path=str((model_dir / "ranking.csv").resolve()) if (model_dir / "ranking.csv").is_file() else None,
                metrics_path=str((model_dir / "metrics.json").resolve()) if (model_dir / "metrics.json").is_file() else None,
                report_path=str((model_dir / "report.json").resolve()) if (model_dir / "report.json").is_file() else None,
            )
            result_count += 1
    else:
        metrics = _read_json(batch_dir / "metrics.json") or {}
        if metrics:
            mid = str(metrics.get("model_id") or batch_id)
            _ingest_model_result_row(
                conn,
                task_id=task_id,
                batch_id=batch_id,
                model_id=mid,
                metrics=metrics,
                ranking_path=str((batch_dir / "ranking.csv").resolve()) if (batch_dir / "ranking.csv").is_file() else None,
                metrics_path=str((batch_dir / "metrics.json").resolve()) if (batch_dir / "metrics.json").is_file() else None,
                report_path=str((batch_dir / "report.json").resolve()) if (batch_dir / "report.json").is_file() else None,
            )
            result_count += 1

    conn.commit()
    return {"task_id": task_id, "batch_id": batch_id, "result_count": result_count, "ingest_status": "ok"}


def ingest_hyperparameter_set(conn: sqlite3.Connection, workspace_root: Path, hp_set_id: str) -> Dict[str, Any]:
    workspace_root = workspace_root.resolve()
    hp_set_root = workspace_root / "hyperparameter_sets" / hp_set_id
    set_meta = _read_json(hp_set_root / "hp_set.json") or {}
    conn.execute(
        """
        INSERT INTO hyperparameter_sets (
            hp_set_id, hp_set_root, name, description, created_at, updated_at,
            ingest_status, ingest_error, ingested_at
        ) VALUES (?, ?, ?, ?, ?, ?, 'ok', NULL, ?)
        ON CONFLICT(hp_set_id) DO UPDATE SET
            hp_set_root=excluded.hp_set_root,
            name=excluded.name,
            description=excluded.description,
            created_at=excluded.created_at,
            updated_at=excluded.updated_at,
            ingest_status=excluded.ingest_status,
            ingest_error=excluded.ingest_error,
            ingested_at=excluded.ingested_at
        """,
        (
            hp_set_id,
            str(hp_set_root.resolve()),
            set_meta.get("name") or hp_set_id,
            set_meta.get("description"),
            set_meta.get("created_at"),
            set_meta.get("updated_at"),
            _utc_now(),
        ),
    )

    versions_root = hp_set_root / "versions"
    if versions_root.is_dir():
        for vdir in sorted(p for p in versions_root.iterdir() if p.is_dir()):
            version_meta = _read_json(vdir / "hp_version.json") or {}
            yaml_path = vdir / "config.yaml"
            conn.execute(
                """
                INSERT INTO hyperparameter_versions (
                    hp_set_id, hp_version_id, created_at, created_by, notes,
                    yaml_path, source_json, ingest_status, ingest_error, ingested_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, 'ok', NULL, ?)
                ON CONFLICT(hp_set_id, hp_version_id) DO UPDATE SET
                    created_at=excluded.created_at,
                    created_by=excluded.created_by,
                    notes=excluded.notes,
                    yaml_path=excluded.yaml_path,
                    source_json=excluded.source_json,
                    ingest_status=excluded.ingest_status,
                    ingest_error=excluded.ingest_error,
                    ingested_at=excluded.ingested_at
                """,
                (
                    hp_set_id,
                    vdir.name,
                    version_meta.get("created_at"),
                    version_meta.get("created_by"),
                    version_meta.get("notes"),
                    str(yaml_path.resolve()) if yaml_path.is_file() else None,
                    json.dumps(version_meta.get("source"), ensure_ascii=False) if version_meta.get("source") else None,
                    _utc_now(),
                ),
            )
    conn.commit()
    return {"hp_set_id": hp_set_id, "ingest_status": "ok"}


def register_run(conn: sqlite3.Connection, run_root: Path) -> Dict[str, Any]:
    """Full run ingest: run + binding + steps + artifacts + models."""
    out: Dict[str, Any] = {}
    out["run"] = ingest_run(conn, run_root)
    out["ga_binding"] = ingest_run_ga_binding(conn, run_root)
    out["steps"] = ingest_step_executions(conn, run_root)
    out["artifacts"] = ingest_key_artifacts(conn, run_root)
    out["models"] = ingest_model(conn, run_root)
    return out


def rebuild_index(conn: sqlite3.Connection, workspace_root: Path) -> Dict[str, Any]:
    workspace_root = workspace_root.resolve()
    summary: Dict[str, Any] = {
        "projects": 0,
        "runs": 0,
        "ga_sets": 0,
        "tasks": 0,
        "batches": 0,
        "hyperparameter_sets": 0,
        "errors": [],
    }

    projects_root = workspace_root / "projects"
    if projects_root.is_dir():
        for project_dir in sorted(p for p in projects_root.iterdir() if p.is_dir()):
            pid = project_dir.name
            ingest_project(conn, workspace_root, pid)
            summary["projects"] += 1
            runs_root = project_dir / "runs"
            if runs_root.is_dir():
                for run_dir in sorted(p for p in runs_root.iterdir() if p.is_dir()):
                    try:
                        register_run(conn, run_dir)
                        summary["runs"] += 1
                    except sqlite3.Error as exc:
                        summary["errors"].append({"run": str(run_dir), "error": str(exc)})

    ga_root = workspace_root / "ga_sets"
    if ga_root.is_dir():
        for ga_dir in sorted(p for p in ga_root.iterdir() if p.is_dir()):
            try:
                ingest_ga_set(conn, workspace_root, ga_dir.name)
                summary["ga_sets"] += 1
            except sqlite3.Error as exc:
                summary["errors"].append({"ga_set": ga_dir.name, "error": str(exc)})

    tasks_root = workspace_root / "model_after_tasks"
    if tasks_root.is_dir():
        for task_dir in sorted(p for p in tasks_root.iterdir() if p.is_dir() and p.name != ".gitkeep"):
            try:
                ingest_model_after_task(conn, task_dir)
                summary["tasks"] += 1
            except sqlite3.Error as exc:
                summary["errors"].append({"task": str(task_dir), "error": str(exc)})

    results_root = workspace_root / "model_after_results"
    if results_root.is_dir():
        for task_dir in sorted(p for p in results_root.iterdir() if p.is_dir()):
            for batch_dir in sorted(p for p in task_dir.iterdir() if p.is_dir()):
                try:
                    ingest_model_after_batch(conn, batch_dir)
                    summary["batches"] += 1
                except sqlite3.Error as exc:
                    summary["errors"].append({"batch": str(batch_dir), "error": str(exc)})

    hp_root = workspace_root / "hyperparameter_sets"
    if hp_root.is_dir():
        for hp_set_dir in sorted(p for p in hp_root.iterdir() if p.is_dir()):
            try:
                ingest_hyperparameter_set(conn, workspace_root, hp_set_dir.name)
                summary["hyperparameter_sets"] += 1
            except sqlite3.Error as exc:
                summary["errors"].append({"hyperparameter_set": hp_set_dir.name, "error": str(exc)})

    return summary
