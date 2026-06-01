from __future__ import annotations

from pathlib import Path

from app.storage.db import init_db
from app.storage import repositories


def _seed_minimal_registry(conn) -> None:
    conn.execute(
        "INSERT INTO projects (project_id, project_root, ingest_status) VALUES (?, ?, 'ok')",
        ("p001", "/tmp/workspace/projects/p001"),
    )
    conn.execute(
        """
        INSERT INTO runs (project_id, run_id, run_root, status, ga_binding_status, ingest_status)
        VALUES (?, ?, ?, ?, ?, 'ok')
        """,
        ("p001", "run_demo", "/tmp/workspace/projects/p001/runs/run_demo", "success", "bound"),
    )
    conn.execute(
        "INSERT INTO ga_sets (ga_set_id, ga_set_root, name, ingest_status) VALUES (?, ?, ?, 'ok')",
        ("fixture", "/tmp/workspace/ga_sets/fixture", "Fixture"),
    )
    conn.execute(
        """
        INSERT INTO ga_versions (
            ga_set_id, ga_version_id, num_ga, checksum, ga_csv_path, ingest_status
        ) VALUES (?, ?, ?, ?, ?, 'ok')
        """,
        ("fixture", "v001", 2, "sha256:abc", "/tmp/workspace/ga_sets/fixture/versions/v001/GA_with_id.csv"),
    )
    conn.execute(
        """
        INSERT INTO run_ga_bindings (
            project_id, run_id, ga_set_id, ga_version_id, checksum, num_ga, status, ingest_status
        ) VALUES (?, ?, ?, ?, ?, ?, ?, 'ok')
        """,
        ("p001", "run_demo", "fixture", "v001", "sha256:abc", 2, "bound"),
    )
    conn.execute(
        """
        INSERT INTO models (
            model_id, project_id, run_id, checkpoint_path, checkpoint_stem, embedding_backend,
            ga_set_id, ga_version_id, ga_binding_checksum, val_mrr, ingest_status
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 'ok')
        """,
        (
            "p001_run_demo_best",
            "p001",
            "run_demo",
            "/tmp/workspace/projects/p001/runs/run_demo/training/ckpts/best.pt",
            "best",
            "ecfp",
            "fixture",
            "v001",
            "sha256:abc",
            0.12,
        ),
    )
    conn.execute(
        """
        INSERT INTO model_after_tasks (
            task_id, task_dir, metal, embedding, candidate_file, positive_file, negative_file, ingest_status
        ) VALUES (?, ?, ?, ?, ?, ?, ?, 'ok')
        """,
        ("t1", "/tmp/workspace/model_after_tasks/t1", "Ni", "ecfp", "candidates.csv", "positives.csv", "negatives.csv"),
    )
    conn.execute(
        """
        INSERT INTO model_after_batches (
            task_id, batch_id, output_dir, command, finished_at, n_models, n_success, best_model_id, ingest_status
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, 'ok')
        """,
        ("t1", "b2", "/tmp/workspace/model_after_results/t1/b2", "evaluate-models", "2026-01-02T00:00:00Z", 2, 1, "raw_best"),
    )
    conn.execute(
        """
        INSERT INTO model_after_batches (
            task_id, batch_id, output_dir, command, finished_at, n_models, n_success, best_model_id, ingest_status
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, 'ok')
        """,
        ("t1", "b1", "/tmp/workspace/model_after_results/t1/b1", "evaluate-model", "2026-01-01T00:00:00Z", 1, 1, "raw_best"),
    )
    conn.execute(
        """
        INSERT INTO model_after_model_results (
            task_id, batch_id, model_id, project_id, run_id, registry_model_id, model_name,
            mrr, hit_at_5, hit_at_10, hit_at_20, hit_at_50, selection_score, rank_among_models,
            status, ranking_path, metrics_path, report_path, ingest_status
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 'ok')
        """,
        (
            "t1",
            "b2",
            "raw_best",
            "p001",
            "run_demo",
            "p001_run_demo_best",
            "best checkpoint",
            0.3,
            0.4,
            0.5,
            0.6,
            0.7,
            0.3,
            1,
            "success",
            "/tmp/ranking.csv",
            "/tmp/metrics.json",
            "/tmp/report.json",
        ),
    )
    conn.execute(
        """
        INSERT INTO model_after_model_results (
            task_id, batch_id, model_id, project_id, run_id, registry_model_id, model_name,
            mrr, hit_at_5, hit_at_10, hit_at_20, hit_at_50, selection_score, rank_among_models,
            status, ranking_path, metrics_path, report_path, ingest_status
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 'ok')
        """,
        (
            "t1",
            "b2",
            "raw_missing",
            "p001",
            "run_demo",
            "missing_registry_model",
            "missing checkpoint",
            0.1,
            0.1,
            0.1,
            0.1,
            0.1,
            0.1,
            2,
            "success",
            "/tmp/ranking2.csv",
            "/tmp/metrics2.json",
            "/tmp/report2.json",
        ),
    )
    conn.execute(
        """
        INSERT INTO artifacts (project_id, run_id, step_id, rel_path, role, exists_flag, ingest_status)
        VALUES (?, ?, ?, ?, ?, ?, 'ok')
        """,
        ("p001", "run_demo", "step1", "tmp/out.csv", "output", 1),
    )
    conn.execute(
        """
        INSERT INTO step_executions (project_id, run_id, step_id, started_at, status, ingest_status)
        VALUES (?, ?, ?, ?, ?, 'ok')
        """,
        ("p001", "run_demo", "step1", "2026-01-01T00:00:00Z", "success"),
    )
    conn.commit()


def test_get_registry_summary(tmp_path: Path) -> None:
    conn = init_db(tmp_path / "db.sqlite")
    _seed_minimal_registry(conn)
    summary = repositories.get_registry_summary(conn)
    conn.close()
    assert summary["projects"] == 1
    assert summary["runs"] == 1
    assert summary["ga_sets"] == 1
    assert summary["ga_versions"] == 1
    assert summary["models"] == 1
    assert summary["tasks"] == 1
    assert summary["batches"] == 2
    assert summary["results"] == 2
    assert summary["artifacts"] == 1
    assert summary["steps"] == 1


def test_get_helpers(tmp_path: Path) -> None:
    conn = init_db(tmp_path / "db.sqlite")
    _seed_minimal_registry(conn)
    ga_version = repositories.get_ga_version(conn, "fixture", "v001")
    task = repositories.get_model_after_task(conn, "t1")
    batch = repositories.get_batch(conn, "t1", "b2")
    model = repositories.get_model(conn, "p001_run_demo_best")
    ga_set = repositories.get_ga_set(conn, "fixture")
    conn.close()
    assert ga_set is not None and ga_set["ga_set_id"] == "fixture"
    assert ga_version is not None and ga_version["ga_csv_path"].endswith("GA_with_id.csv")
    assert task is not None and task["task_id"] == "t1"
    assert batch is not None and batch["batch_id"] == "b2"
    assert model is not None and model["checkpoint_stem"] == "best"


def test_list_all_batches_order(tmp_path: Path) -> None:
    conn = init_db(tmp_path / "db.sqlite")
    _seed_minimal_registry(conn)
    rows = repositories.list_all_batches(conn)
    conn.close()
    assert [r["batch_id"] for r in rows][:2] == ["b2", "b1"]


def test_list_model_results_with_models_registry_join(tmp_path: Path) -> None:
    conn = init_db(tmp_path / "db.sqlite")
    _seed_minimal_registry(conn)
    rows = repositories.list_model_results_with_models(conn, "t1", "b2")
    conn.close()
    assert len(rows) == 2
    matched = [r for r in rows if r["raw_model_id"] == "raw_best"][0]
    missing = [r for r in rows if r["raw_model_id"] == "raw_missing"][0]

    assert matched["registry_model_id"] == "p001_run_demo_best"
    assert matched["model_join_status"] == "matched_model"
    assert matched["ga_set_id"] == "fixture"
    assert matched["ga_version_id"] == "v001"

    assert missing["registry_model_id"] == "missing_registry_model"
    assert missing["model_join_status"] == "missing_model"
    assert missing["ga_set_id"] is None
