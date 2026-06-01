"""Local subprocess job manager for Phase 3B MVP."""

from __future__ import annotations

import argparse
import json
import os
import sqlite3
import subprocess
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from app.storage.db import init_db

JOB_PENDING = "pending"
JOB_RUNNING = "running"
JOB_SUCCESS = "success"
JOB_FAILED = "failed"
JOB_CANCELLED = "cancelled"


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _ensure_jobs_dir(workspace_root: str | Path) -> Path:
    jobs_dir = Path(workspace_root).resolve() / "jobs"
    jobs_dir.mkdir(parents=True, exist_ok=True)
    return jobs_dir


def create_job(
    conn: sqlite3.Connection,
    *,
    job_type: str,
    title: str,
    command: List[str],
    workspace_root: str | Path,
    db_path: str | Path,
    project_id: Optional[str] = None,
    run_id: Optional[str] = None,
    task_id: Optional[str] = None,
    batch_id: Optional[str] = None,
    model_id: Optional[str] = None,
) -> Dict[str, Any]:
    job_id = f"job_{uuid.uuid4().hex[:12]}"
    log_path = _ensure_jobs_dir(workspace_root) / f"{job_id}.log"
    payload = {
        "job_id": job_id,
        "job_type": job_type,
        "status": JOB_PENDING,
        "title": title,
        "command_json": json.dumps(command, ensure_ascii=False),
        "workspace_root": str(Path(workspace_root).resolve()),
        "db_path": str(Path(db_path).resolve()),
        "project_id": project_id,
        "run_id": run_id,
        "task_id": task_id,
        "batch_id": batch_id,
        "model_id": model_id,
        "log_path": str(log_path.resolve()),
        "result_json": None,
        "error": None,
        "created_at": _utc_now(),
        "started_at": None,
        "finished_at": None,
    }
    conn.execute(
        """
        INSERT INTO jobs (
            job_id, job_type, status, title, command_json, workspace_root, db_path,
            project_id, run_id, task_id, batch_id, model_id, log_path, result_json,
            error, created_at, started_at, finished_at
        ) VALUES (
            :job_id, :job_type, :status, :title, :command_json, :workspace_root, :db_path,
            :project_id, :run_id, :task_id, :batch_id, :model_id, :log_path, :result_json,
            :error, :created_at, :started_at, :finished_at
        )
        """,
        payload,
    )
    conn.commit()
    return payload


def mark_running(conn: sqlite3.Connection, job_id: str) -> None:
    conn.execute(
        "UPDATE jobs SET status=?, started_at=?, error=NULL WHERE job_id=?",
        (JOB_RUNNING, _utc_now(), job_id),
    )
    conn.commit()


def mark_success(conn: sqlite3.Connection, job_id: str, result: Optional[Dict[str, Any]] = None) -> None:
    conn.execute(
        "UPDATE jobs SET status=?, result_json=?, finished_at=? WHERE job_id=?",
        (
            JOB_SUCCESS,
            json.dumps(result, ensure_ascii=False) if result is not None else None,
            _utc_now(),
            job_id,
        ),
    )
    conn.commit()


def mark_failed(conn: sqlite3.Connection, job_id: str, error: str, result: Optional[Dict[str, Any]] = None) -> None:
    conn.execute(
        "UPDATE jobs SET status=?, error=?, result_json=?, finished_at=? WHERE job_id=?",
        (
            JOB_FAILED,
            error,
            json.dumps(result, ensure_ascii=False) if result is not None else None,
            _utc_now(),
            job_id,
        ),
    )
    conn.commit()


def run_subprocess_job(
    db_path: str | Path,
    *,
    job_type: str,
    title: str,
    command: List[str],
    workspace_root: str | Path,
    project_id: Optional[str] = None,
    run_id: Optional[str] = None,
    task_id: Optional[str] = None,
    batch_id: Optional[str] = None,
    model_id: Optional[str] = None,
    cwd: Optional[str | Path] = None,
) -> Dict[str, Any]:
    db_path = Path(db_path).resolve()
    conn = init_db(db_path)
    try:
        job = create_job(
            conn,
            job_type=job_type,
            title=title,
            command=command,
            workspace_root=workspace_root,
            db_path=db_path,
            project_id=project_id,
            run_id=run_id,
            task_id=task_id,
            batch_id=batch_id,
            model_id=model_id,
        )
        mark_running(conn, job["job_id"])
        log_path = Path(job["log_path"])
        working_dir = Path(cwd).resolve() if cwd else Path(workspace_root).resolve()
        with open(log_path, "w", encoding="utf-8") as log_f:
            log_f.write(f"# job_id: {job['job_id']}\n")
            log_f.write(f"# command: {' '.join(command)}\n")
            log_f.write(f"# cwd: {working_dir}\n\n")
            log_f.flush()
            proc = subprocess.run(
                command,
                cwd=str(working_dir),
                stdout=log_f,
                stderr=subprocess.STDOUT,
                text=True,
            )
        result = {"exit_code": proc.returncode}
        if proc.returncode == 0:
            mark_success(conn, job["job_id"], result=result)
        else:
            mark_failed(conn, job["job_id"], error=f"process exit code {proc.returncode}", result=result)
        row = conn.execute("SELECT * FROM jobs WHERE job_id=?", (job["job_id"],)).fetchone()
        return dict(row) if row else job
    finally:
        conn.close()


def launch_subprocess_job_async(
    db_path: str | Path,
    *,
    job_type: str,
    title: str,
    command: List[str],
    workspace_root: str | Path,
    project_id: Optional[str] = None,
    run_id: Optional[str] = None,
    task_id: Optional[str] = None,
    batch_id: Optional[str] = None,
    model_id: Optional[str] = None,
    cwd: Optional[str | Path] = None,
) -> Dict[str, Any]:
    db_path = Path(db_path).resolve()
    conn = init_db(db_path)
    try:
        job = create_job(
            conn,
            job_type=job_type,
            title=title,
            command=command,
            workspace_root=workspace_root,
            db_path=db_path,
            project_id=project_id,
            run_id=run_id,
            task_id=task_id,
            batch_id=batch_id,
            model_id=model_id,
        )
    finally:
        conn.close()

    worker_cmd = [
        sys.executable,
        "-m",
        "app.services.job_manager",
        "worker-run",
        "--db-path",
        str(db_path),
        "--job-id",
        job["job_id"],
    ]
    if cwd:
        worker_cmd += ["--cwd", str(Path(cwd).resolve())]
    env = dict(os.environ)
    env["PYTHONPATH"] = f"{Path(__file__).resolve().parents[2]}{os.pathsep}{env.get('PYTHONPATH', '')}"
    subprocess.Popen(
        worker_cmd,
        cwd=str(Path(__file__).resolve().parents[2]),
        env=env,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        start_new_session=True,
    )
    return job


def worker_run_job(db_path: str | Path, *, job_id: str, cwd: Optional[str | Path] = None) -> int:
    conn = init_db(db_path)
    try:
        row = conn.execute("SELECT * FROM jobs WHERE job_id=?", (job_id,)).fetchone()
        if row is None:
            return 1
        job = dict(row)
        command = json.loads(job.get("command_json") or "[]")
        if not command:
            mark_failed(conn, job_id, "empty command")
            return 1

        mark_running(conn, job_id)
        log_path = Path(job["log_path"])
        working_dir = Path(cwd).resolve() if cwd else Path(job.get("workspace_root") or ".").resolve()
        with open(log_path, "a", encoding="utf-8") as log_f:
            log_f.write(f"\n# worker_start: {_utc_now()}\n")
            log_f.flush()
            proc = subprocess.run(
                command,
                cwd=str(working_dir),
                stdout=log_f,
                stderr=subprocess.STDOUT,
                text=True,
            )
        result = {"exit_code": proc.returncode}
        if proc.returncode == 0:
            mark_success(conn, job_id, result=result)
            return 0
        mark_failed(conn, job_id, error=f"process exit code {proc.returncode}", result=result)
        return proc.returncode
    except Exception as exc:
        mark_failed(conn, job_id, error=str(exc))
        return 1
    finally:
        conn.close()


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Job manager worker entry")
    sub = p.add_subparsers(dest="command", required=True)
    w = sub.add_parser("worker-run", help="Execute a queued job")
    w.add_argument("--db-path", required=True)
    w.add_argument("--job-id", required=True)
    w.add_argument("--cwd", default=None)
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "worker-run":
        return worker_run_job(args.db_path, job_id=args.job_id, cwd=args.cwd)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
