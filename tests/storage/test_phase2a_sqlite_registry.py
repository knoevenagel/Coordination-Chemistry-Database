"""Phase 2A SQLite registry tests (tmp_path only; no real workspace)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from app.storage.db import init_db, table_names
from app.storage import ingest
from app.storage import repositories

WEB_ROOT = Path(__file__).resolve().parents[2]
CHEMDB_TMP = WEB_ROOT / "ChemDB" / "tmp"
SRC_TMP = WEB_ROOT / "ChemDB" / "src" / "tmp"
REAL_DB = WEB_ROOT / "workspace" / "chemdb.sqlite"

EXPECTED_TABLES = {
    "projects",
    "runs",
    "ga_sets",
    "ga_versions",
    "run_ga_bindings",
    "step_executions",
    "artifacts",
    "models",
    "model_after_tasks",
    "model_after_batches",
    "model_after_model_results",
}


@pytest.fixture
def db_path(tmp_path: Path) -> Path:
    return tmp_path / "chemdb.sqlite"


@pytest.fixture
def workspace(tmp_path: Path) -> Path:
    ws = tmp_path / "workspace"
    (ws / "projects" / "p001" / "runs" / "run_demo" / "manifests").mkdir(parents=True)
    (ws / "projects" / "p001" / "runs" / "run_demo" / "training" / "ckpts").mkdir(parents=True)
    (ws / "ga_sets").mkdir(parents=True)
    (ws / "model_after_tasks").mkdir(parents=True)
    (ws / "model_after_results").mkdir(parents=True)
    return ws


def test_init_db_creates_tables(db_path: Path) -> None:
    conn = init_db(db_path)
    names = set(table_names(conn))
    conn.close()
    assert EXPECTED_TABLES <= names


def test_register_project_and_run_from_minimal_workspace(db_path: Path, workspace: Path) -> None:
    run_root = workspace / "projects" / "p001" / "runs" / "run_demo"
    (run_root / "manifests" / "run.json").write_text(
        json.dumps(
            {
                "action": "init-run",
                "run_root": str(run_root),
                "created_at": "2026-05-28T00:00:00+00:00",
                "status": "success",
                "pubchem_files": ["A.csv", "B.csv"],
            }
        ),
        encoding="utf-8",
    )
    conn = init_db(db_path)
    ingest.ingest_project(conn, workspace, "p001")
    ingest.ingest_run(conn, run_root)
    conn.close()

    conn = init_db(db_path)
    projects = repositories.list_projects(conn)
    runs = repositories.list_runs(conn, "p001")
    conn.close()
    assert len(projects) == 1
    assert projects[0]["project_id"] == "p001"
    assert len(runs) == 1
    assert runs[0]["run_id"] == "run_demo"
    assert runs[0]["pubchem_file_count"] == 2


def test_register_ga_set_and_versions(db_path: Path, workspace: Path) -> None:
    ga_root = workspace / "ga_sets" / "test_ga"
    vdir = ga_root / "versions" / "v001"
    vdir.mkdir(parents=True)
    (ga_root / "ga_set.json").write_text(
        json.dumps({"ga_set_id": "test_ga", "name": "Test GA", "tags": []}),
        encoding="utf-8",
    )
    (vdir / "GA_with_id.csv").write_text("GA_SMILES,GA_ID\nc1ccccc1,G1\n", encoding="utf-8")
    (vdir / "ga_version.json").write_text(
        json.dumps(
            {
                "ga_set_id": "test_ga",
                "ga_version_id": "v001",
                "num_ga": 1,
                "checksum": "sha256:abc",
            }
        ),
        encoding="utf-8",
    )

    conn = init_db(db_path)
    ingest.ingest_ga_set(conn, workspace, "test_ga")
    conn.close()

    conn = init_db(db_path)
    versions = repositories.list_ga_versions(conn, "test_ga")
    conn.close()
    assert len(versions) == 1
    assert versions[0]["num_ga"] == 1
    assert versions[0]["checksum"] == "sha256:abc"


def test_register_run_ga_binding(db_path: Path, workspace: Path) -> None:
    run_root = workspace / "projects" / "p001" / "runs" / "run_demo"
    (run_root / "manifests").mkdir(parents=True, exist_ok=True)
    ga_root = workspace / "ga_sets" / "fixture"
    (ga_root / "versions" / "v001").mkdir(parents=True)
    (ga_root / "ga_set.json").write_text(json.dumps({"ga_set_id": "fixture"}), encoding="utf-8")

    conn = init_db(db_path)
    ingest.ingest_run(conn, run_root)
    (run_root / "manifests" / "ga_binding.json").write_text(
        json.dumps(
            {
                "ga_set_id": "fixture",
                "ga_version_id": "v001",
                "checksum": "sha256:deadbeef",
                "num_ga": 10,
                "status": "bound",
            }
        ),
        encoding="utf-8",
    )
    ingest.ingest_run_ga_binding(conn, run_root)
    binding = repositories.get_run_ga_binding(conn, "p001", "run_demo")
    conn.close()
    assert binding is not None
    assert binding["ga_set_id"] == "fixture"
    assert binding["status"] == "bound"


def test_legacy_run_without_ga_binding_does_not_fail(db_path: Path, workspace: Path) -> None:
    run_root = workspace / "projects" / "p001" / "runs" / "legacy_run"
    (run_root / "manifests").mkdir(parents=True)
    (run_root / "manifests" / "run.json").write_text('{"status":"success"}', encoding="utf-8")
    (run_root / "manifests" / "step4_5.json").write_text(
        json.dumps({"step_id": "step4_5", "status": "success", "started_at": "t1"}),
        encoding="utf-8",
    )

    conn = init_db(db_path)
    result = ingest.register_run(conn, run_root)
    binding = repositories.get_run_ga_binding(conn, "p001", "legacy_run")
    run_row = repositories.get_run(conn, "p001", "legacy_run")
    conn.close()

    assert result["ga_binding"]["status"] == "legacy_auto_ga"
    assert binding is not None
    assert binding["status"] == "legacy_auto_ga"
    assert binding["ga_set_id"] is None
    assert run_row["ga_binding_status"] == "legacy_auto_ga"


def test_register_step_executions_and_artifacts(db_path: Path, workspace: Path) -> None:
    run_root = workspace / "projects" / "p001" / "runs" / "run_demo"
    (run_root / "manifests").mkdir(parents=True, exist_ok=True)
    (run_root / "tmp").mkdir(exist_ok=True)
    out_csv = run_root / "tmp" / "ligand_data.csv"
    out_csv.write_text("x\n", encoding="utf-8")
    (run_root / "manifests" / "step1.json").write_text(
        json.dumps(
            {
                "step_id": "step1",
                "status": "success",
                "started_at": "2026-05-28T01:00:00+00:00",
                "finished_at": "2026-05-28T01:01:00+00:00",
                "duration_seconds": 60,
                "exit_code": 0,
                "command": ["python", "step1.py"],
                "inputs": [{"path": str(run_root / "data"), "exists": True}],
                "outputs": [{"path": str(out_csv.resolve()), "exists": True, "size_bytes": 2}],
            }
        ),
        encoding="utf-8",
    )

    conn = init_db(db_path)
    ingest.ingest_run(conn, run_root)
    ingest.ingest_step_executions(conn, run_root)
    steps = conn.execute(
        "SELECT * FROM step_executions WHERE project_id='p001' AND run_id='run_demo'"
    ).fetchall()
    arts = conn.execute(
        "SELECT * FROM artifacts WHERE project_id='p001' AND run_id='run_demo'"
    ).fetchall()
    conn.close()
    assert len(steps) == 1
    assert steps[0]["step_id"] == "step1"
    assert len(arts) >= 1


def test_register_model_from_training_files(db_path: Path, workspace: Path) -> None:
    run_root = workspace / "projects" / "p001" / "runs" / "run_demo"
    train = run_root / "training"
    ckpts = train / "ckpts"
    ckpts.mkdir(parents=True, exist_ok=True)
    (train / "index.json").write_text(
        json.dumps(
            {
                "ligand_embedding": {"backend": "ecfp", "path": "/tmp/npz"},
                "metal_embedding": {"path": "/tmp/metal.csv"},
            }
        ),
        encoding="utf-8",
    )
    (train / "config.yaml").write_text(
        "data:\n  embedding: ecfp\n",
        encoding="utf-8",
    )
    (ckpts / "history.json").write_text(
        json.dumps([{"epoch": 3, "mrr": 0.5, "recall_at_1": 0.1, "recall_at_5": 0.2}]),
        encoding="utf-8",
    )
    (ckpts / "best.pt").write_bytes(b"fake-checkpoint")

    conn = init_db(db_path)
    ingest.ingest_run(conn, run_root)
    result = ingest.ingest_model(conn, run_root)
    models = repositories.list_models(conn, "p001", "run_demo")
    conn.close()
    assert result["count"] == 1
    assert len(models) == 1
    assert models[0]["model_id"] == "p001_run_demo_best"
    assert models[0]["val_mrr"] == 0.5


def test_register_task_counts_rows(db_path: Path, workspace: Path) -> None:
    task_dir = workspace / "model_after_tasks" / "tiny_task"
    task_dir.mkdir(parents=True)
    (task_dir / "task.json").write_text(
        json.dumps(
            {
                "task_id": "tiny_task",
                "metal": "Ni",
                "embedding": "ecfp",
                "candidate_file": "candidates.csv",
                "positive_file": "positives.csv",
                "negative_file": "negatives.csv",
            }
        ),
        encoding="utf-8",
    )
    (task_dir / "candidates.csv").write_text("candidate_id,smiles\n1,C\n2,C\n", encoding="utf-8")
    (task_dir / "positives.csv").write_text("positive_id,smiles\n1,C\n", encoding="utf-8")
    (task_dir / "negatives.csv").write_text("negative_id,smiles\n2,C\n", encoding="utf-8")

    conn = init_db(db_path)
    ingest.ingest_model_after_task(conn, task_dir)
    tasks = repositories.list_model_after_tasks(conn)
    conn.close()
    assert tasks[0]["candidate_count"] == 2
    assert tasks[0]["positive_count"] == 1
    assert tasks[0]["negative_count"] == 1


def test_register_model_after_batch_multimodel(db_path: Path, workspace: Path) -> None:
    batch = workspace / "model_after_results" / "tiny_task" / "batch_001"
    models_dir = batch / "models" / "m_a"
    models_dir.mkdir(parents=True)
    (batch / "model_selection_summary.csv").write_text(
        "model_id,model_name,model_run_root,checkpoint_path,embedding_backend,status,"
        "mrr,hit_at_5,hit_at_10,hit_at_20,hit_at_50,positive_score_mean,negative_score_mean,"
        "pos_neg_margin,selection_score,rank_among_models,error\n"
        "m_a,name,"
        f"{workspace / 'projects/p001/runs/run_demo'},"
        f"{workspace / 'projects/p001/runs/run_demo/training/ckpts/best.pt'},"
        "ecfp,success,0.1,0.2,0.3,0.4,0.5,,,,0.1,1,\n",
        encoding="utf-8",
    )
    (batch / "best_model.json").write_text(
        json.dumps({"model_id": "m_a", "selection_score": 0.1}),
        encoding="utf-8",
    )
    (batch / "model_selection_manifest.json").write_text(
        json.dumps({"n_models": 1, "n_success_with_score": 1}),
        encoding="utf-8",
    )
    (models_dir / "metrics.json").write_text(
        json.dumps(
            {
                "model_id": "m_a",
                "mrr": 0.1,
                "hit_at_5": 0.2,
                "status": "success",
                "model_run_root": str(workspace / "projects/p001/runs/run_demo"),
                "checkpoint_path": str(workspace / "projects/p001/runs/run_demo/training/ckpts/best.pt"),
            }
        ),
        encoding="utf-8",
    )

    conn = init_db(db_path)
    ingest.ingest_model_after_batch(conn, batch)
    results = repositories.list_model_results(conn, "tiny_task", "batch_001")
    conn.close()
    assert len(results) == 1
    assert results[0]["project_id"] == "p001"
    assert results[0]["run_id"] == "run_demo"


def test_rebuild_index_is_idempotent(db_path: Path, workspace: Path) -> None:
    run_root = workspace / "projects" / "p001" / "runs" / "run_demo"
    (run_root / "manifests").mkdir(parents=True, exist_ok=True)
    (run_root / "manifests" / "run.json").write_text('{"status":"success","pubchem_files":[]}', encoding="utf-8")

    conn = init_db(db_path)
    ingest.rebuild_index(conn, workspace)
    c1 = conn.execute("SELECT COUNT(*) AS c FROM runs").fetchone()["c"]
    ingest.rebuild_index(conn, workspace)
    c2 = conn.execute("SELECT COUNT(*) AS c FROM runs").fetchone()["c"]
    conn.close()
    assert c1 == c2 == 1


def test_repository_queries(db_path: Path, workspace: Path) -> None:
    run_root = workspace / "projects" / "p001" / "runs" / "run_demo"
    train = run_root / "training" / "ckpts"
    train.mkdir(parents=True, exist_ok=True)
    (run_root / "manifests").mkdir(parents=True, exist_ok=True)
    (run_root / "manifests" / "run.json").write_text('{"status":"success"}', encoding="utf-8")
    (train / "best.pt").write_bytes(b"x")
    (run_root / "training" / "index.json").write_text("{}", encoding="utf-8")

    batch = workspace / "model_after_results" / "t1" / "b1"
    (batch / "models" / "m1").mkdir(parents=True)
    (batch / "best_model.json").write_text(json.dumps({"model_id": "m1"}), encoding="utf-8")
    (batch / "models" / "m1" / "metrics.json").write_text(
        json.dumps({"model_id": "m1", "mrr": 0.9, "status": "success"}),
        encoding="utf-8",
    )

    conn = init_db(db_path)
    ingest.rebuild_index(conn, workspace)
    models = repositories.list_models(conn)
    tasks = repositories.list_model_after_tasks(conn)
    batches = repositories.list_batches(conn, "t1") if tasks else []
    best = repositories.get_best_model_for_batch(conn, "t1", "b1")
    conn.close()
    assert len(models) >= 1
    assert best is not None
    assert best["model_id"] == "m1"


def test_no_write_to_real_workspace_or_repo(
    db_path: Path,
    workspace: Path,
    repo_tmp_snapshot: set[str],
) -> None:
    real_db_exists_before = REAL_DB.exists()
    run_root = workspace / "projects" / "p001" / "runs" / "run_demo"
    (run_root / "manifests").mkdir(parents=True, exist_ok=True)
    (run_root / "manifests" / "run.json").write_text('{"status":"success"}', encoding="utf-8")

    conn = init_db(db_path)
    ingest.rebuild_index(conn, workspace)
    conn.close()

    assert db_path.is_file()
    assert REAL_DB.exists() == real_db_exists_before
    if CHEMDB_TMP.exists():
        after = {p.name for p in CHEMDB_TMP.iterdir()}
        assert after == repo_tmp_snapshot
    assert not SRC_TMP.exists()


@pytest.fixture(scope="session")
def repo_tmp_snapshot() -> set[str]:
    if not CHEMDB_TMP.exists():
        CHEMDB_TMP.mkdir(parents=True, exist_ok=True)
    return {p.name for p in CHEMDB_TMP.iterdir()}
