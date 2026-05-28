"""1% PubChem integration tests (require CHEMDB_INTEGRATION_PUBCHEM_SOURCE)."""

from __future__ import annotations

from pathlib import Path

import pytest

from app.services.ga_registry_manager import apply_ga_to_run
from app.services.orchestrator import run_pipeline, run_step
from tests.run_isolation.integration_helpers import (
    assert_step_success,
    bind_run_ga_fixture,
    init_run_with_pubchem,
    integration_pubchem_source_or_skip,
)

# --- Tier A (Phase 1B minimum) ---


@pytest.mark.integration_1pct
@pytest.mark.integration_1pct_tier_a
def test_1pct_tier_a_init_run_step1_step2(tmp_path: Path) -> None:
    pytest.importorskip("rdkit")
    pubchem_source = integration_pubchem_source_or_skip()
    run_root = tmp_path / "run_1pct_tier_a"
    ctx = init_run_with_pubchem(run_root, pubchem_source)

    m1 = run_step(ctx, "step1")
    assert m1["status"] == "success", m1
    m2 = run_step(ctx, "step2")
    assert m2["status"] == "success", m2

    assert (ctx.tmp_dir / "complex_data.csv").is_file()
    assert (ctx.tmp_dir / "repaired_ligand_data.csv").is_file()
    assert_step_success(ctx, "step1")
    assert_step_success(ctx, "step2")


@pytest.mark.integration_1pct
@pytest.mark.integration_1pct_tier_a
def test_1pct_two_runs_do_not_share_outputs(tmp_path: Path) -> None:
    pytest.importorskip("rdkit")
    pubchem_source = integration_pubchem_source_or_skip()
    ctx1 = init_run_with_pubchem(tmp_path / "run_1pct_a", pubchem_source)
    ctx2 = init_run_with_pubchem(tmp_path / "run_1pct_b", pubchem_source)

    (ctx1.tmp_dir / "sentinel_a.txt").write_text("a", encoding="utf-8")
    (ctx2.tmp_dir / "sentinel_b.txt").write_text("b", encoding="utf-8")

    assert not (ctx1.tmp_dir / "sentinel_b.txt").exists()
    assert not (ctx2.tmp_dir / "sentinel_a.txt").exists()
    assert ctx1.data_dir.resolve() != ctx2.data_dir.resolve()


# --- Tier B (Phase 1E: bind GA + apply-ga, no auto step4_5) ---


@pytest.mark.integration_1pct
@pytest.mark.integration_1pct_tier_b
def test_1pct_tier_b_preprocess_bind_apply_ga(tmp_path: Path) -> None:
    pytest.importorskip("rdkit")
    pubchem_source = integration_pubchem_source_or_skip()
    run_root = tmp_path / "run_1pct_tier_b"
    ctx = init_run_with_pubchem(run_root, pubchem_source)

    summary = run_pipeline(ctx, "step2")
    assert summary["ok"] is True, summary
    for step_id in ("step1", "step2"):
        assert_step_success(ctx, step_id)

    bind_run_ga_fixture(ctx, "fixture_run_demo", "v001")
    assert (ctx.manifest_dir / "ga_binding.json").is_file()
    assert (ctx.tmp_dir / "GA_with_id.csv").is_file()

    apply_ga_to_run(ctx.run_root)
    for name in ("ligand_with_gac.csv", "IRL_filtered.csv"):
        assert (ctx.tmp_dir / name).is_file(), f"missing {name}"


# --- Tier C (optional) ---


@pytest.mark.integration_1pct
@pytest.mark.integration_1pct_tier_c
def test_1pct_tier_c_pipeline_through_step6_7(tmp_path: Path) -> None:
    pytest.importorskip("rdkit")
    pubchem_source = integration_pubchem_source_or_skip()
    ctx = init_run_with_pubchem(tmp_path / "run_1pct_tier_c", pubchem_source)

    pre = run_pipeline(ctx, "step2")
    assert pre["ok"] is True, pre
    bind_run_ga_fixture(ctx, "fixture_run_demo", "v001")
    apply_ga_to_run(ctx.run_root)

    summary = run_pipeline(ctx, "step6_7")
    assert summary["ok"] is True, summary
    assert (ctx.tmp_dir / "fragments.csv").is_file()
    assert_step_success(ctx, "step6_7")


# --- Enhanced (NOT Phase 1B minimum) ---


@pytest.mark.integration_1pct_enhanced
def test_1pct_enhanced_step8(tmp_path: Path) -> None:
    pytest.importorskip("rdkit")
    pubchem_source = integration_pubchem_source_or_skip()
    ctx = init_run_with_pubchem(tmp_path / "run_1pct_enhanced_8", pubchem_source)
    pre = run_pipeline(ctx, "step2")
    assert pre["ok"] is True, pre
    bind_run_ga_fixture(ctx, "fixture_run_demo", "v001")
    apply_ga_to_run(ctx.run_root)
    pre2 = run_pipeline(ctx, "step6_7")
    assert pre2["ok"] is True, pre2

    m = run_step(ctx, "step8")
    assert m["status"] == "success", m
    assert (ctx.tmp_dir / "neo4j_metals.csv").is_file()
