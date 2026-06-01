from __future__ import annotations

import json
from pathlib import Path

import pytest

from app.services import hyperparameter_manager as hpm


def test_create_hyperparameter_version_from_upload(tmp_path: Path) -> None:
    pytest.importorskip("yaml")
    ws = tmp_path / "workspace"
    ws.mkdir(parents=True, exist_ok=True)
    out = hpm.create_hyperparameter_version_from_upload(
        workspace_root=ws,
        hp_set_id="hp_demo",
        yaml_content="batch_size: 32\noptimizer:\n  lr: 0.001\n",
        notes="demo",
    )
    assert out["hp_set_id"] == "hp_demo"
    assert out["hp_version_id"] == "v001"
    assert Path(out["yaml_path"]).is_file()


def test_apply_hyperparameters_to_run(tmp_path: Path) -> None:
    pytest.importorskip("yaml")
    ws = tmp_path / "workspace"
    run_root = ws / "projects" / "p001" / "runs" / "r001"
    (run_root / "training").mkdir(parents=True, exist_ok=True)
    (run_root / "manifests").mkdir(parents=True, exist_ok=True)
    (run_root / "training" / "config.yaml").write_text(
        "batch_size: 16\noptimizer:\n  lr: 0.01\nmodel:\n  hidden_dim: 64\n",
        encoding="utf-8",
    )
    hpm.create_hyperparameter_version_from_upload(
        workspace_root=ws,
        hp_set_id="hp_demo",
        yaml_content="batch_size: 64\noptimizer:\n  lr: 0.002\nnew_field: 1\n",
    )
    out = hpm.apply_hyperparameters_to_run(
        workspace_root=ws,
        run_root=run_root,
        hp_set_id="hp_demo",
        hp_version_id="v001",
    )
    assert "batch_size" in out["applied"]
    assert "optimizer.lr" in out["applied"]
    assert "new_field" in out["ignored"]
    binding = json.loads((run_root / "manifests" / "hyperparameter_binding.json").read_text(encoding="utf-8"))
    assert binding["hp_set_id"] == "hp_demo"
