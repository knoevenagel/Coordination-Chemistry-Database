from __future__ import annotations

import sys
from pathlib import Path

from app.services import job_manager
from app.storage.db import init_db
from app.storage import repositories


def test_run_subprocess_job_success(tmp_path: Path) -> None:
    workspace = tmp_path / "workspace"
    db_path = workspace / "chemdb.sqlite"
    workspace.mkdir(parents=True, exist_ok=True)
    init_db(db_path).close()

    row = job_manager.run_subprocess_job(
        db_path,
        job_type="test",
        title="success job",
        command=[sys.executable, "-c", "print('ok')"],
        workspace_root=workspace,
        cwd=workspace,
    )

    assert row["status"] == job_manager.JOB_SUCCESS
    assert row["started_at"] is not None
    assert row["finished_at"] is not None
    log_path = Path(row["log_path"])
    assert log_path.is_file()
    assert "ok" in log_path.read_text(encoding="utf-8")


def test_run_subprocess_job_failed(tmp_path: Path) -> None:
    workspace = tmp_path / "workspace"
    db_path = workspace / "chemdb.sqlite"
    workspace.mkdir(parents=True, exist_ok=True)
    init_db(db_path).close()

    row = job_manager.run_subprocess_job(
        db_path,
        job_type="test",
        title="failed job",
        command=[sys.executable, "-c", "import sys; print('bad'); sys.exit(2)"],
        workspace_root=workspace,
        cwd=workspace,
    )

    assert row["status"] == job_manager.JOB_FAILED
    assert "exit code 2" in (row["error"] or "")
    log_path = Path(row["log_path"])
    assert log_path.is_file()
    assert "bad" in log_path.read_text(encoding="utf-8")


def test_jobs_repositories(tmp_path: Path) -> None:
    workspace = tmp_path / "workspace"
    db_path = workspace / "chemdb.sqlite"
    workspace.mkdir(parents=True, exist_ok=True)
    conn = init_db(db_path)
    created = job_manager.create_job(
        conn,
        job_type="demo",
        title="demo title",
        command=[sys.executable, "-c", "print('x')"],
        workspace_root=workspace,
        db_path=db_path,
    )
    fetched = repositories.get_job(conn, created["job_id"])
    jobs = repositories.list_jobs(conn)
    recent = repositories.get_recent_jobs(conn, limit=1)
    conn.close()

    assert fetched is not None
    assert fetched["job_id"] == created["job_id"]
    assert len(jobs) >= 1
    assert len(recent) == 1


def test_worker_run_job(tmp_path: Path) -> None:
    workspace = tmp_path / "workspace"
    db_path = workspace / "chemdb.sqlite"
    workspace.mkdir(parents=True, exist_ok=True)
    conn = init_db(db_path)
    created = job_manager.create_job(
        conn,
        job_type="worker",
        title="worker title",
        command=[sys.executable, "-c", "print('worker-ok')"],
        workspace_root=workspace,
        db_path=db_path,
    )
    conn.close()

    code = job_manager.worker_run_job(db_path, job_id=created["job_id"], cwd=workspace)
    assert code == 0

    conn = init_db(db_path)
    row = repositories.get_job(conn, created["job_id"])
    conn.close()
    assert row is not None
    assert row["status"] == job_manager.JOB_SUCCESS
