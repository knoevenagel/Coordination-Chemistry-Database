"""Shared fixtures for run isolation tests."""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Set

import pytest

from app.services.run_context import RunContext

WEB_ROOT = Path(__file__).resolve().parents[2]
CHEMDB_ROOT = WEB_ROOT / "ChemDB"
REPO_TMP = CHEMDB_ROOT / "tmp"
SRC_TMP = CHEMDB_ROOT / "src" / "tmp"
FIXTURES = Path(__file__).parent / "fixtures"

_STREAM_MARKERS = {
    "integration",
    "integration_1pct",
    "integration_1pct_tier_a",
    "integration_1pct_tier_b",
    "integration_1pct_tier_c",
    "integration_1pct_enhanced",
    "integration_model_after_rank",
    "integration_model_after_eval",
    "integration_model_after_models",
}


@pytest.fixture(autouse=True)
def _stream_step_output_for_integration_tests(request) -> None:
    """Mirror step logs to stdout during real subprocess integration tests."""
    names = {m.name for m in request.node.iter_markers()}
    if not names & _STREAM_MARKERS:
        yield
        return
    prev_stream = os.environ.get("CHEMDB_STREAM_STEP_OUTPUT")
    prev_sort = os.environ.get("CHEMDB_STEP1_SORT_BY_LENGTH")
    prev_row_to = os.environ.get("CHEMDB_STEP1_ROW_TIMEOUT_SEC")
    os.environ["CHEMDB_STREAM_STEP_OUTPUT"] = "1"
    # Worker count: use step defaults (min(cpu, 200)); team confirmed some rows need long RDKit time.
    # Long SMILES first: reduces tqdm appearing stuck at 19436/19437 while a straggler finishes.
    os.environ["CHEMDB_STEP1_SORT_BY_LENGTH"] = "1"
    # Per-row cap (seconds); 0 = disabled. Prevents infinite hang on pathological structures.
    os.environ["CHEMDB_STEP1_ROW_TIMEOUT_SEC"] = "600"
    yield
    for key, prev in (
        ("CHEMDB_STREAM_STEP_OUTPUT", prev_stream),
        ("CHEMDB_STEP1_SORT_BY_LENGTH", prev_sort),
        ("CHEMDB_STEP1_ROW_TIMEOUT_SEC", prev_row_to),
    ):
        if prev is None:
            os.environ.pop(key, None)
        else:
            os.environ[key] = prev


@pytest.fixture(scope="session")
def repo_tmp_snapshot() -> Set[str]:
    if not REPO_TMP.exists():
        REPO_TMP.mkdir(parents=True, exist_ok=True)
    return {p.name for p in REPO_TMP.iterdir()}


@pytest.fixture(autouse=True)
def assert_repo_tmp_unchanged(repo_tmp_snapshot: Set[str]) -> None:
    yield
    if REPO_TMP.exists():
        after = {p.name for p in REPO_TMP.iterdir()}
        assert after == repo_tmp_snapshot, f"repo tmp/ changed: before={repo_tmp_snapshot}, after={after}"
    assert not SRC_TMP.exists(), f"unexpected {SRC_TMP}"


def seed_run_from_fixtures(run_root: Path) -> RunContext:
    ctx = RunContext.from_run_root(run_root)
    ctx.create_dirs()
    shutil.copy2(FIXTURES / "metal_list.txt", ctx.data_dir / "metal_list.txt")
    shutil.copy2(FIXTURES / "p_elements_list.txt", ctx.data_dir / "p_elements_list.txt")
    shutil.copy2(FIXTURES / "minimal_pubchem.csv", ctx.data_dir / "pubchem" / "minimal_pubchem.csv")
    return ctx


@pytest.fixture
def run_root(tmp_path: Path) -> Path:
    root = tmp_path / "run_001"
    seed_run_from_fixtures(root)
    return root


@pytest.fixture
def run_ctx(run_root: Path) -> RunContext:
    return RunContext.from_run_root(run_root)


def mock_run_subprocess(handler):
    """Patch orchestrator subprocess runner (works with stream mode)."""
    from unittest.mock import patch

    def fake_subprocess(cmd, *, cwd, env, log_f, stream_to_console):
        handler(cmd, log_f)
        return 0

    return patch("app.services.orchestrator._run_subprocess", side_effect=fake_subprocess)
