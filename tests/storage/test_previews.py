from __future__ import annotations

import json
from pathlib import Path

from app.storage.db import init_db
from app.storage import previews


def test_read_csv_preview_limit(tmp_path: Path) -> None:
    csv_path = tmp_path / "big.csv"
    rows = ["a,b"] + [f"{i},{i}" for i in range(200)]
    csv_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    out = previews.read_csv_preview(csv_path, limit=5)
    assert out["ok"] is True
    assert len(out["data"]) == 5


def test_read_csv_preview_missing_file(tmp_path: Path) -> None:
    out = previews.read_csv_preview(tmp_path / "missing.csv", limit=10)
    assert out["ok"] is False
    assert out["data"].empty


def test_read_json_preview_missing_file(tmp_path: Path) -> None:
    out = previews.read_json_preview(tmp_path / "missing.json")
    assert out["ok"] is False
    assert out["data"] == {}


def test_resolve_task_file(tmp_path: Path) -> None:
    task_row = {
        "task_dir": str(tmp_path / "task1"),
        "candidate_file": "candidates.csv",
        "positive_file": "positives.csv",
        "negative_file": "negatives.csv",
    }
    p = previews.resolve_task_file(task_row, "candidates")
    assert p is not None
    assert str(p).endswith("task1/candidates.csv")


def test_preview_ga_version(tmp_path: Path) -> None:
    conn = init_db(tmp_path / "db.sqlite")
    ga_csv = tmp_path / "GA_with_id.csv"
    ga_csv.write_text("GA_SMILES,GA_ID\nc1ccccc1,G1\nc1ccncc1,G2\n", encoding="utf-8")

    conn.execute(
        "INSERT INTO ga_sets (ga_set_id, ga_set_root, name, ingest_status) VALUES (?, ?, ?, 'ok')",
        ("fixture", str(tmp_path / "ga_sets/fixture"), "Fixture"),
    )
    conn.execute(
        """
        INSERT INTO ga_versions (ga_set_id, ga_version_id, num_ga, checksum, ga_csv_path, ingest_status)
        VALUES (?, ?, ?, ?, ?, 'ok')
        """,
        ("fixture", "v001", 2, "sha256:abc", str(ga_csv)),
    )
    conn.commit()

    out = previews.preview_ga_version(conn, "fixture", "v001", limit=1)
    conn.close()
    assert out["ok"] is True
    assert len(out["data"]) == 1


def test_preview_task_files_handles_missing_negative(tmp_path: Path) -> None:
    task_dir = tmp_path / "task1"
    task_dir.mkdir(parents=True)
    (task_dir / "candidates.csv").write_text("id,smiles\n1,C\n", encoding="utf-8")
    (task_dir / "positives.csv").write_text("id,smiles\n1,C\n", encoding="utf-8")

    task_row = {
        "task_dir": str(task_dir),
        "candidate_file": "candidates.csv",
        "positive_file": "positives.csv",
        "negative_file": "",
    }
    out = previews.preview_task_files(task_row, limit=100)
    assert out["candidates"]["ok"] is True
    assert out["positives"]["ok"] is True
    assert out["negatives"]["ok"] is False


def test_preview_model_result(tmp_path: Path) -> None:
    ranking = tmp_path / "ranking.csv"
    metrics = tmp_path / "metrics.json"
    report = tmp_path / "report.json"
    ranking.write_text("rank,score\n1,0.1\n2,0.2\n", encoding="utf-8")
    metrics.write_text(json.dumps({"mrr": 0.1}), encoding="utf-8")
    report.write_text(json.dumps({"status": "success"}), encoding="utf-8")

    row = {
        "ranking_path": str(ranking),
        "metrics_path": str(metrics),
        "report_path": str(report),
    }
    out = previews.preview_model_result(row, limit=1)
    assert out["ranking"]["ok"] is True
    assert len(out["ranking"]["data"]) == 1
    assert out["metrics"]["ok"] is True
    assert out["report"]["ok"] is True
