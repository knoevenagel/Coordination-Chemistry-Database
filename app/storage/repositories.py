"""Query helpers for Phase 2A SQLite registry."""

from __future__ import annotations

import sqlite3
from typing import Any, Dict, List, Optional


def _rows_to_dicts(rows: List[sqlite3.Row]) -> List[Dict[str, Any]]:
    return [dict(r) for r in rows]


def _table_exists(conn: sqlite3.Connection, table: str) -> bool:
    row = conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
        (table,),
    ).fetchone()
    return row is not None


def list_projects(conn: sqlite3.Connection) -> List[Dict[str, Any]]:
    rows = conn.execute("SELECT * FROM projects ORDER BY project_id").fetchall()
    return _rows_to_dicts(rows)


def list_jobs(conn: sqlite3.Connection) -> List[Dict[str, Any]]:
    rows = conn.execute(
        "SELECT * FROM jobs ORDER BY created_at DESC, job_id DESC"
    ).fetchall()
    return _rows_to_dicts(rows)


def get_job(conn: sqlite3.Connection, job_id: str) -> Optional[Dict[str, Any]]:
    row = conn.execute(
        "SELECT * FROM jobs WHERE job_id=?",
        (job_id,),
    ).fetchone()
    return dict(row) if row else None


def get_recent_jobs(conn: sqlite3.Connection, limit: int = 20) -> List[Dict[str, Any]]:
    rows = conn.execute(
        "SELECT * FROM jobs ORDER BY created_at DESC, job_id DESC LIMIT ?",
        (max(1, int(limit)),),
    ).fetchall()
    return _rows_to_dicts(rows)


def get_registry_summary(conn: sqlite3.Connection) -> Dict[str, int]:
    table_alias = {
        "projects": "projects",
        "runs": "runs",
        "ga_sets": "ga_sets",
        "ga_versions": "ga_versions",
        "models": "models",
        "tasks": "model_after_tasks",
        "batches": "model_after_batches",
        "results": "model_after_model_results",
        "artifacts": "artifacts",
        "steps": "step_executions",
        "hyperparameter_sets": "hyperparameter_sets",
        "hyperparameter_versions": "hyperparameter_versions",
    }
    out: Dict[str, int] = {}
    for alias, table in table_alias.items():
        if _table_exists(conn, table):
            row = conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()
            out[alias] = int(row["n"]) if row else 0
        else:
            out[alias] = 0
    return out


def list_runs(conn: sqlite3.Connection, project_id: Optional[str] = None) -> List[Dict[str, Any]]:
    if project_id:
        rows = conn.execute(
            "SELECT * FROM runs WHERE project_id=? ORDER BY run_id",
            (project_id,),
        ).fetchall()
    else:
        rows = conn.execute("SELECT * FROM runs ORDER BY project_id, run_id").fetchall()
    return _rows_to_dicts(rows)


def get_run(conn: sqlite3.Connection, project_id: str, run_id: str) -> Optional[Dict[str, Any]]:
    row = conn.execute(
        "SELECT * FROM runs WHERE project_id=? AND run_id=?",
        (project_id, run_id),
    ).fetchone()
    return dict(row) if row else None


def list_ga_sets(conn: sqlite3.Connection) -> List[Dict[str, Any]]:
    rows = conn.execute("SELECT * FROM ga_sets ORDER BY ga_set_id").fetchall()
    return _rows_to_dicts(rows)


def get_ga_set(conn: sqlite3.Connection, ga_set_id: str) -> Optional[Dict[str, Any]]:
    row = conn.execute(
        "SELECT * FROM ga_sets WHERE ga_set_id=?",
        (ga_set_id,),
    ).fetchone()
    return dict(row) if row else None


def list_ga_versions(conn: sqlite3.Connection, ga_set_id: str) -> List[Dict[str, Any]]:
    rows = conn.execute(
        "SELECT * FROM ga_versions WHERE ga_set_id=? ORDER BY ga_version_id",
        (ga_set_id,),
    ).fetchall()
    return _rows_to_dicts(rows)


def get_ga_version(
    conn: sqlite3.Connection,
    ga_set_id: str,
    ga_version_id: str,
) -> Optional[Dict[str, Any]]:
    row = conn.execute(
        "SELECT * FROM ga_versions WHERE ga_set_id=? AND ga_version_id=?",
        (ga_set_id, ga_version_id),
    ).fetchone()
    return dict(row) if row else None


def get_run_ga_binding(conn: sqlite3.Connection, project_id: str, run_id: str) -> Optional[Dict[str, Any]]:
    row = conn.execute(
        "SELECT * FROM run_ga_bindings WHERE project_id=? AND run_id=?",
        (project_id, run_id),
    ).fetchone()
    return dict(row) if row else None


def list_models(
    conn: sqlite3.Connection,
    project_id: Optional[str] = None,
    run_id: Optional[str] = None,
) -> List[Dict[str, Any]]:
    if project_id and run_id:
        rows = conn.execute(
            "SELECT * FROM models WHERE project_id=? AND run_id=? ORDER BY model_id",
            (project_id, run_id),
        ).fetchall()
    elif project_id:
        rows = conn.execute(
            "SELECT * FROM models WHERE project_id=? ORDER BY run_id, model_id",
            (project_id,),
        ).fetchall()
    else:
        rows = conn.execute("SELECT * FROM models ORDER BY project_id, run_id, model_id").fetchall()
    return _rows_to_dicts(rows)


def get_model(conn: sqlite3.Connection, model_id: str) -> Optional[Dict[str, Any]]:
    row = conn.execute(
        "SELECT * FROM models WHERE model_id=?",
        (model_id,),
    ).fetchone()
    return dict(row) if row else None


def list_model_after_tasks(conn: sqlite3.Connection) -> List[Dict[str, Any]]:
    rows = conn.execute("SELECT * FROM model_after_tasks ORDER BY task_id").fetchall()
    return _rows_to_dicts(rows)


def get_model_after_task(conn: sqlite3.Connection, task_id: str) -> Optional[Dict[str, Any]]:
    row = conn.execute(
        "SELECT * FROM model_after_tasks WHERE task_id=?",
        (task_id,),
    ).fetchone()
    return dict(row) if row else None


def list_batches(conn: sqlite3.Connection, task_id: str) -> List[Dict[str, Any]]:
    rows = conn.execute(
        "SELECT * FROM model_after_batches WHERE task_id=? ORDER BY batch_id",
        (task_id,),
    ).fetchall()
    return _rows_to_dicts(rows)


def list_all_batches(conn: sqlite3.Connection) -> List[Dict[str, Any]]:
    rows = conn.execute(
        """
        SELECT * FROM model_after_batches
        ORDER BY task_id, finished_at DESC, batch_id
        """
    ).fetchall()
    return _rows_to_dicts(rows)


def get_batch(conn: sqlite3.Connection, task_id: str, batch_id: str) -> Optional[Dict[str, Any]]:
    row = conn.execute(
        """
        SELECT * FROM model_after_batches
        WHERE task_id=? AND batch_id=?
        """,
        (task_id, batch_id),
    ).fetchone()
    return dict(row) if row else None


def list_model_results(conn: sqlite3.Connection, task_id: str, batch_id: str) -> List[Dict[str, Any]]:
    rows = conn.execute(
        """
        SELECT * FROM model_after_model_results
        WHERE task_id=? AND batch_id=?
        ORDER BY rank_among_models, model_id
        """,
        (task_id, batch_id),
    ).fetchall()
    return _rows_to_dicts(rows)


def list_model_results_with_models(
    conn: sqlite3.Connection,
    task_id: str,
    batch_id: str,
) -> List[Dict[str, Any]]:
    rows = conn.execute(
        """
        SELECT
          r.task_id,
          r.batch_id,
          r.model_id AS raw_model_id,
          r.registry_model_id,
          CASE
            WHEN m.model_id IS NULL THEN 'missing_model'
            ELSE 'matched_model'
          END AS model_join_status,
          r.model_name,
          r.project_id,
          r.run_id,
          m.checkpoint_stem,
          m.embedding_backend,
          m.ga_set_id,
          m.ga_version_id,
          m.ga_binding_checksum,
          m.val_mrr,
          r.mrr,
          r.hit_at_5,
          r.hit_at_10,
          r.hit_at_20,
          r.hit_at_50,
          r.selection_score,
          r.rank_among_models,
          r.status,
          r.ranking_path,
          r.metrics_path,
          r.report_path
        FROM model_after_model_results r
        LEFT JOIN models m
          ON r.registry_model_id = m.model_id
        WHERE r.task_id=? AND r.batch_id=?
        ORDER BY r.rank_among_models, r.model_id
        """,
        (task_id, batch_id),
    ).fetchall()
    return _rows_to_dicts(rows)


def list_runs_by_ga_version(
    conn: sqlite3.Connection,
    ga_set_id: str,
    ga_version_id: str,
) -> List[Dict[str, Any]]:
    rows = conn.execute(
        """
        SELECT
          b.project_id,
          b.run_id,
          b.status AS binding_status,
          b.bound_at,
          b.checksum,
          r.run_root,
          r.status AS run_status,
          r.ga_binding_status
        FROM run_ga_bindings b
        LEFT JOIN runs r
          ON b.project_id = r.project_id AND b.run_id = r.run_id
        WHERE b.ga_set_id=? AND b.ga_version_id=?
        ORDER BY b.project_id, b.run_id
        """,
        (ga_set_id, ga_version_id),
    ).fetchall()
    return _rows_to_dicts(rows)


def list_hyperparameter_sets(conn: sqlite3.Connection) -> List[Dict[str, Any]]:
    if not _table_exists(conn, "hyperparameter_sets"):
        return []
    rows = conn.execute(
        "SELECT * FROM hyperparameter_sets ORDER BY hp_set_id"
    ).fetchall()
    return _rows_to_dicts(rows)


def list_hyperparameter_versions(conn: sqlite3.Connection, hp_set_id: str) -> List[Dict[str, Any]]:
    if not _table_exists(conn, "hyperparameter_versions"):
        return []
    rows = conn.execute(
        "SELECT * FROM hyperparameter_versions WHERE hp_set_id=? ORDER BY hp_version_id",
        (hp_set_id,),
    ).fetchall()
    return _rows_to_dicts(rows)


def get_hyperparameter_version(conn: sqlite3.Connection, hp_set_id: str, hp_version_id: str) -> Optional[Dict[str, Any]]:
    if not _table_exists(conn, "hyperparameter_versions"):
        return None
    row = conn.execute(
        "SELECT * FROM hyperparameter_versions WHERE hp_set_id=? AND hp_version_id=?",
        (hp_set_id, hp_version_id),
    ).fetchone()
    return dict(row) if row else None


def get_best_model_for_batch(
    conn: sqlite3.Connection,
    task_id: str,
    batch_id: str,
) -> Optional[Dict[str, Any]]:
    batch = conn.execute(
        "SELECT * FROM model_after_batches WHERE task_id=? AND batch_id=?",
        (task_id, batch_id),
    ).fetchone()
    if not batch:
        return None
    best_id = batch["best_model_id"]
    if best_id:
        row = conn.execute(
            """
            SELECT * FROM model_after_model_results
            WHERE task_id=? AND batch_id=? AND model_id=?
            """,
            (task_id, batch_id, best_id),
        ).fetchone()
        if row:
            return dict(row)
    row = conn.execute(
        """
        SELECT * FROM model_after_model_results
        WHERE task_id=? AND batch_id=?
        ORDER BY selection_score DESC NULLS LAST, mrr DESC NULLS LAST
        LIMIT 1
        """,
        (task_id, batch_id),
    ).fetchone()
    return dict(row) if row else None
