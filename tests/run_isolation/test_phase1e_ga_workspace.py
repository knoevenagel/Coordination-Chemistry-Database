"""Phase 1E: workspace GA sets, run binding, core pipeline without auto GA."""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from unittest.mock import patch

import pytest

from app.services import ga_registry_manager as grm
from app.services.orchestrator import run_pipeline, run_step
from app.services.run_context import RunContext
from app.services.step_registry import PIPELINE_ORDER
from tests.run_isolation.conftest import FIXTURES, REPO_TMP, SRC_TMP, mock_run_subprocess, seed_run_from_fixtures
from tests.run_isolation.integration_helpers import bind_run_ga_fixture, install_ga_fixture_to_workspace

GA_FIXTURES = FIXTURES / "ga_sets"
CHEMDB_ROOT = Path(__file__).resolve().parents[2] / "ChemDB"


@pytest.fixture
def ga_workspace(tmp_path, monkeypatch):
    ws = tmp_path / "workspace"
    ga_root = ws / "ga_sets"
    for name in ("fixture_run_demo", "fixture_restructured"):
        shutil.copytree(GA_FIXTURES / name, ga_root / name)
    monkeypatch.setattr(grm, "WORKSPACE_ROOT", ws)
    monkeypatch.setattr(grm, "GA_SETS_ROOT", ga_root)
    return ga_root


def _ctx_with_repaired(tmp_path: Path) -> RunContext:
    ctx = seed_run_from_fixtures(tmp_path / "run")
    (ctx.tmp_dir / "repaired_ligand_data.csv").write_text(
        "ligand_new_did,source_did,new_smiles,old_smiles\n"
        "LIG_D1,D1,c1ccccc1,c1ccccc1\n"
        "LIG_D2,D2,c1ccncc1,c1ccncc1\n",
        encoding="utf-8",
    )
    return ctx


def test_pipeline_order_excludes_step4_5():
    assert "step4_5" not in PIPELINE_ORDER
    assert "require_bound_ga" in PIPELINE_ORDER
    assert "apply_ga_to_run" in PIPELINE_ORDER


def test_create_ga_version_from_csv_immutable(ga_workspace, tmp_path):
    pytest.importorskip("rdkit")
    src = ga_workspace / "fixture_run_demo" / "versions" / "v001" / "GA_with_id.csv"
    out = grm.create_ga_version_from_csv("fixture_run_demo", src)
    v1 = ga_workspace / "fixture_run_demo" / "versions" / "v001" / "GA_with_id.csv"
    v2 = ga_workspace / "fixture_run_demo" / "versions" / out["ga_version_id"] / "GA_with_id.csv"
    assert v2.is_file()
    assert out["ga_version_id"] == "v002"
    assert v1.read_text(encoding="utf-8") == (GA_FIXTURES / "fixture_run_demo" / "versions" / "v001" / "GA_with_id.csv").read_text(
        encoding="utf-8"
    )


def test_invalid_ga_smiles_rejected(ga_workspace, tmp_path):
    pytest.importorskip("rdkit")
    bad = tmp_path / "bad_ga.csv"
    bad.write_text("GA_SMILES,GA_ID\nnot_a_smiles,G000001\n", encoding="utf-8")
    with pytest.raises(ValueError, match="invalid GA_SMILES"):
        grm.create_ga_version_from_csv("fixture_run_demo", bad)
    assert not (ga_workspace / "fixture_run_demo" / "versions" / "v003").exists()


def test_missing_ga_id_generated(ga_workspace, tmp_path):
    pytest.importorskip("rdkit")
    raw = tmp_path / "input.csv"
    raw.write_text("GA_SMILES,GA_ID\nc1ccccc1,\n", encoding="utf-8")
    out = grm.create_ga_version_from_csv("fixture_run_demo", raw)
    csv_path = ga_workspace / "fixture_run_demo" / "versions" / out["ga_version_id"] / "GA_with_id.csv"
    text = csv_path.read_text(encoding="utf-8")
    assert "G" in text.splitlines()[1].split(",")[1]


def test_bind_materializes_snapshot(ga_workspace, tmp_path):
    ctx = _ctx_with_repaired(tmp_path)
    binding = grm.bind_ga_version_to_run(ctx.run_root, "fixture_run_demo", "v001")
    mat = ctx.tmp_dir / "GA_with_id.csv"
    assert mat.is_file()
    src = ga_workspace / "fixture_run_demo" / "versions" / "v001" / "GA_with_id.csv"
    assert grm.ga_csv_checksum(mat)[0] == grm.ga_csv_checksum(src)[0]
    assert binding["checksum"] == grm.ga_csv_checksum(mat)[0]
    assert (ctx.manifest_dir / "ga_binding.json").is_file()


def test_two_runs_same_version_independent_snapshots(ga_workspace, tmp_path):
    ctx_a = _ctx_with_repaired(tmp_path / "a")
    ctx_b = _ctx_with_repaired(tmp_path / "b")
    grm.bind_ga_version_to_run(ctx_a.run_root, "fixture_run_demo", "v001")
    grm.bind_ga_version_to_run(ctx_b.run_root, "fixture_run_demo", "v001")
    ga_a = ctx_a.tmp_dir / "GA_with_id.csv"
    ga_b = ctx_b.tmp_dir / "GA_with_id.csv"
    assert ga_a.is_file() and ga_b.is_file()
    assert ga_a.resolve() != ga_b.resolve()
    ga_b.write_text("GA_SMILES,GA_ID\nc1ccccc1,G999999\n", encoding="utf-8")
    assert "G874135" in ga_a.read_text(encoding="utf-8")


def test_two_runs_different_versions(ga_workspace, tmp_path):
    ctx_a = _ctx_with_repaired(tmp_path / "a")
    ctx_b = _ctx_with_repaired(tmp_path / "b")
    grm.bind_ga_version_to_run(ctx_a.run_root, "fixture_run_demo", "v001")
    grm.bind_ga_version_to_run(ctx_b.run_root, "fixture_restructured", "v001")
    n_a = len((ctx_a.tmp_dir / "GA_with_id.csv").read_text(encoding="utf-8").splitlines()) - 1
    n_b = len((ctx_b.tmp_dir / "GA_with_id.csv").read_text(encoding="utf-8").splitlines()) - 1
    assert n_a < n_b


def test_checksum_mismatch_fails_require_bound(ga_workspace, tmp_path):
    ctx = _ctx_with_repaired(tmp_path)
    grm.bind_ga_version_to_run(ctx.run_root, "fixture_run_demo", "v001")
    (ctx.tmp_dir / "GA_with_id.csv").write_text("GA_SMILES,GA_ID\nc1ccccc1,G000001\n", encoding="utf-8")
    result = grm.require_bound_ga(ctx)
    assert result["ok"] is False
    assert "checksum" in result["error"].lower()


def test_pipeline_without_binding_fails(tmp_path):
    ctx = _ctx_with_repaired(tmp_path)
    summary = run_pipeline(ctx, "require_bound_ga")
    assert summary["ok"] is False
    assert not (ctx.tmp_dir / "ligand_with_gac.csv").is_file()


def test_core_pipeline_does_not_auto_generate_ga(tmp_path):
    ctx = _ctx_with_repaired(tmp_path)

    def handler(cmd, log_f):
        log_f.write("mock\n")
        if "step1" in cmd[1] or cmd[1].endswith("step1.py"):
            (ctx.tmp_dir / "complex_data.csv").write_text("x\n", encoding="utf-8")
            (ctx.tmp_dir / "ligand_data.csv").write_text(
                "ligand_did,ligand_smiles,source_complex_did\nD1,C,CX\n", encoding="utf-8"
            )
        elif "step2" in cmd[1]:
            (ctx.tmp_dir / "repaired_ligand_data.csv").write_text(
                "ligand_new_did,source_did,new_smiles,old_smiles\nD1,D1,C,C\n", encoding="utf-8"
            )

    with mock_run_subprocess(handler):
        summary = run_pipeline(ctx, "step2")
    assert summary["ok"] is True
    assert not (ctx.tmp_dir / "GA_with_id.csv").is_file()
    assert not (ctx.tmp_dir / "ligand_with_gac.csv").is_file()


def test_rebind_writes_stale_report(ga_workspace, tmp_path):
    ctx = _ctx_with_repaired(tmp_path)
    grm.bind_ga_version_to_run(ctx.run_root, "fixture_run_demo", "v001")
    (ctx.tmp_dir / "ligand_with_gac.csv").write_text("DID,SMILES,GAC\n", encoding="utf-8")
    (ctx.tmp_dir / "IRL_filtered.csv").write_text("DID,SMILES,GAC\n", encoding="utf-8")
    grm.bind_ga_version_to_run(ctx.run_root, "fixture_restructured", "v001")
    stale = ctx.manifest_dir / "ga_stale_report.json"
    assert stale.is_file()
    data = json.loads(stale.read_text(encoding="utf-8"))
    paths = {e["path"] for e in data["stale_files"]}
    assert "tmp/ligand_with_gac.csv" in paths
    assert (ctx.tmp_dir / "ligand_with_gac.csv").is_file()


def test_apply_ga_uses_run_tmp_ga(ga_workspace, tmp_path):
    pytest.importorskip("rdkit")
    ctx = _ctx_with_repaired(tmp_path)
    grm.bind_ga_version_to_run(ctx.run_root, "fixture_run_demo", "v001")
    calls = []

    def fake_run(mode, input_dir, output_dir, **kw):
        calls.append((mode, Path(output_dir) / "GA_with_id.csv"))
        if mode != "apply-ga":
            raise AssertionError(mode)
        if not (Path(output_dir) / "GA_with_id.csv").is_file():
            raise FileNotFoundError("GA_with_id.csv")
        (Path(output_dir) / "ligand_with_gac.csv").write_text(
            "DID,SMILES,GAC\nLIG_D1,c1ccccc1,5\n", encoding="utf-8"
        )
        (Path(output_dir) / "IRL_filtered.csv").write_text(
            "DID,SMILES,GAC\nLIG_D1,c1ccccc1,5\n", encoding="utf-8"
        )
        (Path(output_dir) / "IRL_filtered_cleaned.csv").write_text(
            "DID,SMILES,GAC\n", encoding="utf-8"
        )

    with patch.object(grm, "_run_step4_5", side_effect=fake_run):
        grm.apply_ga_to_run(ctx.run_root)
    assert calls[0][0] == "apply-ga"
    assert calls[0][1].parent == ctx.tmp_dir.resolve()


def test_through_step4_5_rejected():
    from app.services.orchestrator import _validate_through

    with pytest.raises(SystemExit, match="no longer runs step4_5"):
        _validate_through("step4_5", "core")


def test_no_write_repo_tmp_on_bind(ga_workspace, tmp_path):
    ctx = _ctx_with_repaired(tmp_path)
    grm.bind_ga_version_to_run(ctx.run_root, "fixture_run_demo", "v001")
    assert not SRC_TMP.exists()
    if REPO_TMP.exists():
        assert not any(REPO_TMP.glob("GA_with_id.csv"))


def test_require_bound_ga_step_success(ga_workspace, tmp_path):
    ctx = _ctx_with_repaired(tmp_path)
    grm.bind_ga_version_to_run(ctx.run_root, "fixture_run_demo", "v001")
    m = run_step(ctx, "require_bound_ga")
    assert m["status"] == "success"
