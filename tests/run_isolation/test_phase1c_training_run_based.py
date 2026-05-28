"""Phase 1C: per-run ECFP embedding + training pipeline isolation tests."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Set
import numpy as np
import pytest
import yaml

from app.services.orchestrator import build_command, init_run, run_step
from app.services.run_context import RunContext
from app.services.step_registry import STEP_REGISTRY, TRAINING_PIPELINE_ORDER
from app.services import training_manager
from tests.run_isolation.conftest import mock_run_subprocess

WEB_ROOT = Path(__file__).resolve().parents[2]
CHEMDB_ROOT = WEB_ROOT / "ChemDB"
BUILD_SCRIPT = CHEMDB_ROOT / "src" / "tools" / "build_L3_embedding_index.py"
ZSCORE_SCRIPT = CHEMDB_ROOT / "data" / "metal_embedding" / "zscore_element_features.py"
REPO_L3_EMB = CHEMDB_ROOT / "data" / "L3_embedding"
REPO_METAL_EMB = CHEMDB_ROOT / "data" / "metal_embedding"
REPO_TRAINING = CHEMDB_ROOT / "training"

REPO_ARTIFACT_PATTERNS = [
    REPO_L3_EMB / "*.npz",
    REPO_METAL_EMB / "element_features_zscore.csv",
    REPO_TRAINING / "index.json",
    REPO_TRAINING / "config.yaml",
    REPO_TRAINING / "split_index.json",
    REPO_TRAINING / "*_records.pkl",
    REPO_TRAINING / "ckpts",
    REPO_TRAINING / "evaluation" / "results",
]


def _snapshot_repo_artifacts() -> Dict[str, Set[str]]:
    snap: Dict[str, Set[str]] = {}
    for pattern in REPO_ARTIFACT_PATTERNS:
        key = str(pattern)
        if pattern.suffix == ".npz" or pattern.name.endswith(".csv") or pattern.suffix in (".json", ".yaml", ".pkl"):
            parent = pattern.parent
            glob_pat = pattern.name
            snap[key] = {p.name for p in parent.glob(glob_pat)} if parent.exists() else set()
        elif pattern.name == "ckpts":
            snap[key] = (
                {str(p.relative_to(pattern)) for p in pattern.rglob("*") if p.is_file()}
                if pattern.exists()
                else set()
            )
        elif "evaluation" in str(pattern):
            snap[key] = (
                {str(p.relative_to(pattern)) for p in pattern.rglob("*") if p.is_file()}
                if pattern.exists()
                else set()
            )
    return snap


def _assert_repo_artifacts_unchanged(before: Dict[str, Set[str]], after: Dict[str, Set[str]]) -> None:
    assert before == after, f"repo training/embedding artifacts changed: before={before}, after={after}"


def _seed_training_inputs(ctx: RunContext) -> None:
    ctx.create_dirs()
    (ctx.tmp_dir / "step13_kl_nl_samples.csv").write_text(
        "target_did,metal_symbol,kl_did,nl_did\nL3_T1,Fe,D1,D2\n",
        encoding="utf-8",
    )
    (ctx.tmp_dir / "m_l3_pairs.csv").write_text(
        "ligand_did,metal_symbol\nL3_T1,Fe\n",
        encoding="utf-8",
    )
    (ctx.tmp_dir / "l3_gac.json").write_text("{}", encoding="utf-8")
    # build_L3: metal_l3_index uses L3_*; repaired keys strip L3_ prefix for lookup
    (ctx.tmp_dir / "repaired_ligand_data.csv").write_text(
        "ligand_new_did,source_did,new_smiles,old_smiles\nD1,D1,C,C\n",
        encoding="utf-8",
    )
    (ctx.tmp_dir / "metal_l3_index.csv").write_text("l3_did\nL3_D1\n", encoding="utf-8")
    emb_dir = ctx.data_dir / "L3_embedding"
    emb_dir.mkdir(parents=True, exist_ok=True)
    np.savez(
        emb_dir / "L3_embedding_ecfp.npz",
        dids=np.array(["D1"], dtype=object),
        smiles=np.array(["C"], dtype=object),
        embeddings=np.zeros((1, 2048), dtype=np.float32),
    )
    metal_dir = ctx.data_dir / "metal_embedding"
    metal_dir.mkdir(parents=True, exist_ok=True)
    (metal_dir / "element_features.csv").write_text(
        "element,symbol,val\nFe,Fe,1.0\n",
        encoding="utf-8",
    )
    (metal_dir / "element_features_zscore.csv").write_text(
        "element,symbol,val\nFe,Fe,0.0\n",
        encoding="utf-8",
    )


def test_init_run_copies_metal_embedding_seed(tmp_path: Path) -> None:
    ctx = RunContext.from_run_root(tmp_path / "run_metal_seed")
    init_run(ctx)
    dst = ctx.data_dir / "metal_embedding" / "element_features.csv"
    assert dst.is_file()
    assert dst.stat().st_size > 0


def test_zscore_cli_writes_run_output(tmp_path: Path) -> None:
    before = _snapshot_repo_artifacts()
    ctx = RunContext.from_run_root(tmp_path / "run_zscore")
    ctx.create_dirs()
    inp = ctx.data_dir / "metal_embedding" / "element_features.csv"
    inp.parent.mkdir(parents=True, exist_ok=True)
    inp.write_text("element,symbol,x\nFe,Fe,10\nCu,Cu,20\n", encoding="utf-8")
    out = ctx.data_dir / "metal_embedding" / "element_features_zscore.csv"
    subprocess.run(
        [sys.executable, str(ZSCORE_SCRIPT), "--input", str(inp), "--output", str(out)],
        check=True,
        cwd=str(CHEMDB_ROOT),
    )
    assert out.is_file()
    repo_zscore = REPO_METAL_EMB / "element_features_zscore.csv"
    assert repo_zscore.is_file()
    after = _snapshot_repo_artifacts()
    _assert_repo_artifacts_unchanged(before, after)


def test_build_l3_embedding_ecfp_outputs_under_run(tmp_path: Path) -> None:
    pytest.importorskip("rdkit")
    before = _snapshot_repo_artifacts()
    ctx = RunContext.from_run_root(tmp_path / "run_ecfp")
    ctx.create_dirs()
    # build_L3: metal_l3_index uses L3_*; repaired keys strip L3_ prefix for lookup
    (ctx.tmp_dir / "repaired_ligand_data.csv").write_text(
        "ligand_new_did,source_did,new_smiles,old_smiles\nD1,D1,C,C\n",
        encoding="utf-8",
    )
    (ctx.tmp_dir / "metal_l3_index.csv").write_text("l3_did\nL3_D1\n", encoding="utf-8")
    out_dir = ctx.data_dir / "L3_embedding"
    subprocess.run(
        [
            sys.executable,
            str(BUILD_SCRIPT),
            "--tmp-dir",
            str(ctx.tmp_dir),
            "--out-dir",
            str(out_dir),
            "--backends",
            "ecfp",
            "--max-samples",
            "1",
            "--chunk-size",
            "1",
        ],
        check=True,
        cwd=str(CHEMDB_ROOT),
        env={**__import__("os").environ, "PYTHONPATH": str(CHEMDB_ROOT / "src")},
    )
    out_npz = out_dir / "L3_embedding_ecfp.npz"
    assert out_npz.is_file()
    after = _snapshot_repo_artifacts()
    _assert_repo_artifacts_unchanged(before, after)


def test_no_molclr_required_for_phase1c() -> None:
    text = BUILD_SCRIPT.read_text(encoding="utf-8")
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("def ") or stripped.startswith("class "):
            break
        if stripped.startswith(("import ", "from ")):
            assert "molclr" not in stripped.lower()


def test_prepare_training_index_uses_ecfp(tmp_path: Path) -> None:
    ctx = RunContext.from_run_root(tmp_path / "run_index")
    _seed_training_inputs(ctx)
    result = training_manager.prepare_training_index(ctx)
    assert result["ok"] is True
    index = json.loads((ctx.training_dir / "index.json").read_text(encoding="utf-8"))
    le = index["ligand_embedding"]
    assert le["backend"] == "ecfp"
    assert le["mode"] is None
    assert le["path"] == str((ctx.data_dir / "L3_embedding" / "L3_embedding_ecfp.npz").resolve())
    for key in ("samples_csv", "m_l3_pairs", "l3_gac", "m_l3_pairs_path"):
        assert str(ctx.run_root.resolve()) in index[key]
    assert str(ctx.run_root.resolve()) in index["kl_nl_ul_index"]["path"]


def test_prepare_training_config_paths_under_run(tmp_path: Path) -> None:
    ctx = RunContext.from_run_root(tmp_path / "run_cfg")
    _seed_training_inputs(ctx)
    training_manager.prepare_training_index(ctx)
    result = training_manager.prepare_training_config(ctx, integration=True)
    assert result["ok"] is True
    cfg = yaml.safe_load((ctx.training_dir / "config.yaml").read_text(encoding="utf-8"))
    assert str(ctx.training_dir.resolve()) == cfg["data"]["training_dir"]
    assert str((ctx.training_dir / "index.json").resolve()) == cfg["data"]["index_path"]
    assert str((ctx.training_dir / "ckpts").resolve()) == cfg["train"]["ckpt_dir"]
    assert cfg["train"]["epochs"] == 3
    assert cfg["train"]["early_stop_patience"] == 0
    assert cfg["train"]["scheduler"] is None
    assert str(cfg["train"]["device"]).startswith("cuda")
    assert "ChemDB/training" not in yaml.dump(cfg)


def test_training_registry_steps() -> None:
    for step_id in TRAINING_PIPELINE_ORDER:
        assert step_id in STEP_REGISTRY
    assert len(TRAINING_PIPELINE_ORDER) == 6


def test_training_manifest_paths_under_run_mock(tmp_path: Path) -> None:
    ctx = RunContext.from_run_root(tmp_path / "run_manifest_train")
    init_run(ctx)
    _seed_training_inputs(ctx)

    def handler(cmd, log_f):
        log_f.write("ok\n")
        cmd_s = " ".join(cmd)
        if "zscore" in cmd_s:
            out = ctx.data_dir / "metal_embedding" / "element_features_zscore.csv"
            out.write_text("element,symbol,x\nFe,Fe,0\n", encoding="utf-8")
        elif "build_L3" in cmd_s:
            out_dir = ctx.data_dir / "L3_embedding"
            out_dir.mkdir(parents=True, exist_ok=True)
            np.savez(
                out_dir / "L3_embedding_ecfp.npz",
                dids=np.array(["D1"], dtype=object),
                smiles=np.array(["C"], dtype=object),
                embeddings=np.zeros((1, 8), dtype=np.float32),
            )
        elif "training.data" in cmd_s:
            for name in ("split_index.json", "train_records.pkl", "val_records.pkl", "test_records.pkl"):
                (ctx.training_dir / name).write_bytes(b"stub")
        elif "training.train" in cmd_s:
            ckpt = ctx.training_dir / "ckpts"
            ckpt.mkdir(parents=True, exist_ok=True)
            (ckpt / "best.pt").write_bytes(b"pt")
            (ckpt / "history.json").write_text("{}", encoding="utf-8")

    with mock_run_subprocess(handler):
        for step_id in TRAINING_PIPELINE_ORDER:
            m = run_step(ctx, step_id)
            assert m["status"] == "success", step_id
            man = json.loads((ctx.manifest_dir / f"{step_id}.json").read_text(encoding="utf-8"))
            assert (ctx.log_dir / f"{step_id}.log").is_file()
            for entry in man.get("inputs", []) + man.get("outputs", []):
                p = Path(entry["path"])
                assert ctx.resolve_under_run(p), f"{step_id}: {p} not under run_root"


def test_no_write_to_repo_embedding_training(tmp_path: Path) -> None:
    before = _snapshot_repo_artifacts()
    ctx = RunContext.from_run_root(tmp_path / "run_no_repo_write")
    init_run(ctx)
    _seed_training_inputs(ctx)
    training_manager.prepare_training_index(ctx)
    training_manager.prepare_training_config(ctx, integration=True)
    (ctx.training_dir / "split_index.json").write_text("{}", encoding="utf-8")
    after = _snapshot_repo_artifacts()
    _assert_repo_artifacts_unchanged(before, after)


def test_two_runs_ecfp_training_isolated(tmp_path: Path) -> None:
    run1 = RunContext.from_run_root(tmp_path / "run_001")
    run2 = RunContext.from_run_root(tmp_path / "run_002")
    for ctx in (run1, run2):
        _seed_training_inputs(ctx)
        training_manager.prepare_training_index(ctx)
        training_manager.prepare_training_config(ctx, integration=True)
        (ctx.training_dir / "train_records.pkl").write_bytes(b"r1" if ctx is run1 else b"r2")
        ckpt = ctx.training_dir / "ckpts" / "best.pt"
        ckpt.parent.mkdir(parents=True, exist_ok=True)
        ckpt.write_bytes(b"ckpt1" if ctx is run1 else b"ckpt2")
    assert (run1.training_dir / "train_records.pkl") != (run2.training_dir / "train_records.pkl")
    assert (run1.training_dir / "train_records.pkl").read_bytes() == b"r1"
    assert (run2.training_dir / "train_records.pkl").read_bytes() == b"r2"
    idx1 = json.loads((run1.training_dir / "index.json").read_text(encoding="utf-8"))
    idx2 = json.loads((run2.training_dir / "index.json").read_text(encoding="utf-8"))
    assert run1.run_root.as_posix() in idx1["ligand_embedding"]["path"]
    assert run2.run_root.as_posix() in idx2["ligand_embedding"]["path"]
    assert idx1["ligand_embedding"]["path"] != idx2["ligand_embedding"]["path"]


def test_build_l3_command_uses_run_paths(tmp_path: Path) -> None:
    ctx = RunContext.from_run_root(tmp_path / "run_cmd")
    ctx.create_dirs()
    spec = STEP_REGISTRY["build_l3_embedding_ecfp"]
    cmd = build_command(spec, ctx)
    assert str(ctx.tmp_dir) in " ".join(cmd)
    assert str(ctx.data_dir / "L3_embedding") in " ".join(cmd)
    assert "--backends" in cmd and "ecfp" in cmd


@pytest.mark.integration_embedding
def test_integration_embedding_real(tmp_path: Path) -> None:
    pytest.importorskip("rdkit")
    ctx = RunContext.from_run_root(tmp_path / "run_int_emb")
    init_run(ctx)
    (ctx.tmp_dir / "repaired_ligand_data.csv").write_text(
        "ligand_new_did,source_did,new_smiles,old_smiles\nL3_D1,D1,CCO,CCO\n",
        encoding="utf-8",
    )
    (ctx.tmp_dir / "metal_l3_index.csv").write_text("l3_did\nL3_D1\n", encoding="utf-8")
    m1 = run_step(ctx, "zscore_metal_embedding")
    assert m1["status"] == "success"
    m2 = run_step(ctx, "build_l3_embedding_ecfp")
    assert m2["status"] == "success"
    assert (ctx.data_dir / "metal_embedding" / "element_features_zscore.csv").is_file()
    assert (ctx.data_dir / "L3_embedding" / "L3_embedding_ecfp.npz").is_file()


@pytest.mark.integration_training_data
def test_integration_training_data_real(tmp_path: Path) -> None:
    pytest.importorskip("torch")
    ctx = RunContext.from_run_root(tmp_path / "run_int_data")
    init_run(ctx)
    _seed_training_inputs(ctx)
    training_manager.prepare_training_index(ctx)
    m = run_step(ctx, "training_data")
    assert m["status"] == "success"
    for name in ("split_index.json", "train_records.pkl", "val_records.pkl", "test_records.pkl"):
        p = ctx.training_dir / name
        assert p.is_file(), name
        assert ctx.resolve_under_run(p)


@pytest.mark.integration_training_train
def test_integration_training_train_real(tmp_path: Path) -> None:
    pytest.importorskip("torch")
    ctx = RunContext.from_run_root(tmp_path / "run_int_train")
    init_run(ctx)
    _seed_training_inputs(ctx)
    import os

    os.environ["CHEMDB_TRAINING_INTEGRATION"] = "1"
    try:
        training_manager.prepare_training_index(ctx)
        training_manager.prepare_training_config(ctx, integration=True)
        run_step(ctx, "training_data")
        m = run_step(ctx, "training_train")
    finally:
        os.environ.pop("CHEMDB_TRAINING_INTEGRATION", None)
    assert m["status"] == "success"
    assert (ctx.training_dir / "ckpts" / "best.pt").is_file()
    assert (ctx.training_dir / "ckpts" / "history.json").is_file()
