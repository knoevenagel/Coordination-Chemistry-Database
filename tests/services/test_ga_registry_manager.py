from __future__ import annotations

import json
from pathlib import Path

import pytest

from app.services import ga_registry_manager as grm


def test_create_ga_set_or_version_from_upload_strict(tmp_path: Path) -> None:
    pytest.importorskip("rdkit")
    ws = tmp_path / "workspace"
    ws.mkdir(parents=True, exist_ok=True)
    csv_text = "GA_SMILES,GA_ID\nc1ccccc1,G000001\nCCO,G000002\n"

    result = grm.create_ga_set_or_version_from_upload(
        "upload_demo",
        csv_text,
        workspace_root=ws,
    )

    assert result["ga_set_id"] == "upload_demo"
    assert result["ga_version_id"] == "v001"
    csv_path = Path(result["ga_csv"])
    assert csv_path.is_file()
    assert "GA_SMILES,GA_ID" in csv_path.read_text(encoding="utf-8")
    meta = json.loads((csv_path.parent / "ga_version.json").read_text(encoding="utf-8"))
    assert meta["source"]["type"] == "upload"


def test_create_ga_set_or_version_from_upload_rejects_missing_ga_id(tmp_path: Path) -> None:
    pytest.importorskip("rdkit")
    ws = tmp_path / "workspace"
    ws.mkdir(parents=True, exist_ok=True)
    csv_text = "GA_SMILES,GA_ID\nc1ccccc1,\n"

    with pytest.raises(ValueError):
        grm.create_ga_set_or_version_from_upload(
            "upload_demo",
            csv_text,
            workspace_root=ws,
        )

    assert not (ws / "ga_sets" / "upload_demo" / "versions").exists()
