from __future__ import annotations

from pathlib import Path

from app.services import model_build_manager


def test_submit_model_build_job_composes_command(monkeypatch, tmp_path: Path) -> None:
    captured = {}

    def _fake_launch(db_path, **kwargs):
        captured["db_path"] = str(db_path)
        captured["kwargs"] = kwargs
        return {"job_id": "job_test_123"}

    monkeypatch.setattr(
        "app.services.model_build_manager.job_manager.launch_subprocess_job_async",
        _fake_launch,
    )

    ws = tmp_path / "workspace"
    db = ws / "chemdb.sqlite"
    out = model_build_manager.submit_model_build_job(
        workspace_root=ws,
        db_path=db,
        project_id="p001",
        run_id="run_new",
        ga_set_id="fixture_run_demo",
        ga_version_id="v001",
        pubchem_source=str(ws / "pubchem"),
    )

    assert out["job_id"] == "job_test_123"
    cmd = captured["kwargs"]["command"]
    assert cmd[0:2] == ["bash", "-lc"]
    script = cmd[2]
    assert "app.services.orchestrator init-run" in script
    assert "app.services.ga_registry_manager bind-ga-version-to-run" in script
    assert "app.storage.cli rebuild-index" in script
