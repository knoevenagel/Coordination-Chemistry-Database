from __future__ import annotations

from pathlib import Path

from app.services import model_after_manager as mam


def test_submit_evaluate_models_job(monkeypatch, tmp_path: Path) -> None:
    ws = tmp_path / "workspace"
    task_dir = ws / "model_after_tasks" / "task_a"
    task_dir.mkdir(parents=True, exist_ok=True)
    (task_dir / "task.json").write_text("{}", encoding="utf-8")
    captured = {}

    def _fake_launch(db_path, **kwargs):
        captured["db_path"] = str(db_path)
        captured["kwargs"] = kwargs
        return {"job_id": "job_eval_1"}

    monkeypatch.setattr("app.services.model_after_manager.job_manager.launch_subprocess_job_async", _fake_launch)
    out = mam.submit_evaluate_models_job(
        db_path=ws / "chemdb.sqlite",
        workspace_root=ws,
        task_id="task_a",
        models=[
            {
                "model_id": "m1",
                "project_id": "p001",
                "run_id": "r001",
                "checkpoint_path": "/tmp/best.pt",
                "checkpoint_stem": "best",
            }
        ],
        batch_id="batch_x",
        device="cpu",
    )
    assert out["job"]["job_id"] == "job_eval_1"
    assert "evaluate-models" in captured["kwargs"]["command"][2]


def test_submit_recommend_job(monkeypatch, tmp_path: Path) -> None:
    ws = tmp_path / "workspace"
    task_dir = ws / "model_after_tasks" / "task_a"
    run_root = ws / "projects" / "p001" / "runs" / "r001"
    task_dir.mkdir(parents=True, exist_ok=True)
    run_root.mkdir(parents=True, exist_ok=True)
    captured = {}

    def _fake_launch(db_path, **kwargs):
        captured["kwargs"] = kwargs
        return {"job_id": "job_rec_1"}

    monkeypatch.setattr("app.services.model_after_manager.job_manager.launch_subprocess_job_async", _fake_launch)
    out = mam.submit_recommend_job(
        db_path=ws / "chemdb.sqlite",
        workspace_root=ws,
        task_dir=task_dir,
        model_run_root=run_root,
        model_id="m1",
        batch_id="recommend_x",
    )
    assert out["job"]["job_id"] == "job_rec_1"
    assert "recommend" in captured["kwargs"]["command"][2]
