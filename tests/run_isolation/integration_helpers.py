"""Helpers for integration and 1% integration tests."""

from __future__ import annotations

import json
import os
import shutil
from pathlib import Path
from typing import Any, Dict, List

import pytest

from app.services.run_context import RunContext
from app.services.orchestrator import init_run

WEB_ROOT = Path(__file__).resolve().parents[2]
FIXTURES = Path(__file__).parent / "fixtures"
GA_FIXTURES = FIXTURES / "ga_sets"
WORKSPACE_GA_SETS = WEB_ROOT / "workspace" / "ga_sets"


def _path_under_base(child: Path, base: Path) -> bool:
    """True if child is under base (Python 3.8 compatible; no Path.is_relative_to)."""
    try:
        child.resolve().relative_to(base.resolve())
        return True
    except ValueError:
        return False


def integration_pubchem_source_or_skip() -> Path:
    raw = os.environ.get("CHEMDB_INTEGRATION_PUBCHEM_SOURCE")
    if not raw:
        pytest.skip("CHEMDB_INTEGRATION_PUBCHEM_SOURCE is not set")
    path = Path(raw).resolve()
    if not path.exists():
        pytest.skip(f"CHEMDB_INTEGRATION_PUBCHEM_SOURCE not found: {path}")
    if path.is_file():
        if path.suffix.lower() != ".csv":
            pytest.skip(f"CHEMDB_INTEGRATION_PUBCHEM_SOURCE must be .csv or directory: {path}")
        return path
    if not path.is_dir():
        pytest.skip(f"CHEMDB_INTEGRATION_PUBCHEM_SOURCE invalid: {path}")
    if not any(path.glob("*.csv")):
        pytest.skip(f"No *.csv under CHEMDB_INTEGRATION_PUBCHEM_SOURCE: {path}")
    manifest = path / "MANIFEST.json"
    if not manifest.is_file():
        pytest.skip(f"MANIFEST.json missing under {path}")
    return path


def init_run_with_pubchem(run_root: Path, pubchem_source: Path) -> RunContext:
    ctx = RunContext.from_run_root(run_root)
    init_run(ctx, pubchem_source=pubchem_source)
    pubchem_dir = ctx.data_dir / "pubchem"
    assert pubchem_dir.is_dir(), "run data/pubchem must exist after init-run"
    assert any(pubchem_dir.glob("*.csv")), "run data/pubchem must contain CSVs"
    # step1 must read run-local pubchem, not external integration_data
    for csv in pubchem_dir.glob("*.csv"):
        assert _path_under_base(csv, ctx.run_root), f"pubchem file must be under run_root: {csv}"
    if pubchem_source.is_dir():
        ext = pubchem_source.resolve()
        for csv in pubchem_dir.glob("*.csv"):
            assert not csv.resolve().samefile(ext / csv.name), (
                "pubchem must be copied into run_root, not same inode as integration_data"
            )
    return ctx


def init_run_minimal(run_root: Path) -> RunContext:
    ctx = RunContext.from_run_root(run_root)
    init_run(ctx, pubchem_source=FIXTURES / "minimal_pubchem.csv")
    return ctx


def install_ga_fixture_to_workspace(ga_set_id: str = "fixture_run_demo") -> Path:
    """Copy tests/fixtures/ga_sets/{id} into workspace/ga_sets/{id} for bind CLI."""
    src = GA_FIXTURES / ga_set_id
    if not src.is_dir():
        raise FileNotFoundError(f"missing GA fixture: {src}")
    dest = WORKSPACE_GA_SETS / ga_set_id
    dest.mkdir(parents=True, exist_ok=True)
    for name in ("ga_set.json", "versions"):
        s = src / name
        d = dest / name
        if s.is_dir():
            if d.exists():
                shutil.rmtree(d)
            shutil.copytree(s, d)
        elif s.is_file():
            shutil.copy2(s, d)
    return dest


def bind_run_ga_fixture(
    ctx: RunContext,
    ga_set_id: str = "fixture_run_demo",
    ga_version_id: str = "v001",
) -> Dict[str, Any]:
    install_ga_fixture_to_workspace(ga_set_id)
    from app.services.ga_registry_manager import bind_ga_version_to_run

    return bind_ga_version_to_run(ctx.run_root, ga_set_id, ga_version_id)


def load_manifest(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def assert_manifest_under_run(manifest: Dict[str, Any], ctx: RunContext) -> None:
    run_s = str(ctx.run_root.resolve())
    assert manifest.get("run_root") == run_s or run_s in str(manifest.get("run_root", ""))
    log_file = manifest.get("log_file")
    if log_file:
        assert _path_under_base(Path(log_file), ctx.run_root)
    cmd = manifest.get("command")
    if isinstance(cmd, list):
        cmd_str = " ".join(str(x) for x in cmd)
        assert run_s in cmd_str or str(ctx.data_dir) in cmd_str or str(ctx.tmp_dir) in cmd_str
    if manifest.get("paths_under_run_root") is not None:
        assert manifest["paths_under_run_root"] is True
    inp = manifest.get("input_check") or {}
    for p in inp.get("checked") or []:
        assert str(ctx.run_root.resolve()) in str(p) or _path_under_base(Path(p), ctx.run_root)


def assert_step_success(ctx: RunContext, step_id: str) -> Dict[str, Any]:
    path = ctx.manifest_dir / f"{step_id}.json"
    assert path.is_file(), f"missing manifest {path}"
    m = load_manifest(path)
    assert m.get("status") == "success", m
    assert_manifest_under_run(m, ctx)
    assert (ctx.log_dir / f"{step_id}.log").is_file()
    return m


def collect_paths_from_manifest(manifest: Dict[str, Any]) -> List[str]:
    paths: List[str] = []
    for key in ("log_file", "run_root"):
        v = manifest.get(key)
        if v:
            paths.append(str(v))
    out = manifest.get("output_check") or {}
    for rel in out.get("present") or []:
        paths.append(rel)
    return paths
