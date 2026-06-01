from __future__ import annotations

import json
from pathlib import Path

import pytest

from app.services import task_registry_manager as trm


def test_create_task_from_uploads_success(tmp_path: Path) -> None:
    ws = tmp_path / "workspace"
    ws.mkdir(parents=True, exist_ok=True)
    result = trm.create_task_from_uploads(
        task_id="task_demo_1",
        workspace_root=ws,
        candidates_csv="id,smiles\n1,C\n2,CC\n",
        positives_csv="id,smiles\n1,C\n",
        negatives_csv="id,smiles\n3,CCC\n",
        metal="Fe",
        embedding="ecfp",
        notes="demo",
    )
    assert result["task_id"] == "task_demo_1"
    task_dir = Path(result["task_dir"])
    assert (task_dir / "task.json").is_file()
    payload = json.loads((task_dir / "task.json").read_text(encoding="utf-8"))
    assert payload["candidate_count"] == 2
    assert payload["positive_count"] == 1
    assert payload["negative_count"] == 1


def test_create_task_from_uploads_rollback_on_invalid_csv(tmp_path: Path) -> None:
    ws = tmp_path / "workspace"
    ws.mkdir(parents=True, exist_ok=True)
    with pytest.raises(ValueError):
        trm.create_task_from_uploads(
            task_id="task_demo_2",
            workspace_root=ws,
            candidates_csv="id,smiles\n",
            positives_csv="id,smiles\n1,C\n",
            negatives_csv=None,
        )
    assert not (ws / "model_after_tasks" / "task_demo_2").exists()
