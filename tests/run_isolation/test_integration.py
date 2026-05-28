"""Real subprocess integration tests using minimal_pubchem fixture only."""

from __future__ import annotations

from pathlib import Path

import pytest

from app.services.orchestrator import run_step
from tests.run_isolation.integration_helpers import (
    assert_step_success,
    init_run_minimal,
)

pytestmark = pytest.mark.integration


def test_integration_minimal_step1_step2(tmp_path: Path) -> None:
    pytest.importorskip("rdkit")
    run_root = tmp_path / "run_integration_minimal"
    ctx = init_run_minimal(run_root)

    m1 = run_step(ctx, "step1")
    assert m1["status"] == "success", m1
    assert (ctx.tmp_dir / "ligand_data.csv").is_file()
    assert (ctx.data_dir / "pubchem" / "minimal_pubchem.csv").is_file()

    m2 = run_step(ctx, "step2")
    assert m2["status"] == "success", m2
    assert (ctx.tmp_dir / "repaired_ligand_data.csv").is_file()

    assert_step_success(ctx, "step1")
    assert_step_success(ctx, "step2")
