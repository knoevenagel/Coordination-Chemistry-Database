"""Phase 1D: run-based model-after (recommend, evaluate, multi-model selection)."""

from __future__ import annotations

import csv
import json
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Set
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
import yaml

WEB_ROOT = Path(__file__).resolve().parents[2]
CHEMDB_ROOT = WEB_ROOT / "ChemDB"
WORKSPACE = WEB_ROOT / "workspace"
SHARED_TASK_NI = WORKSPACE / "model_after_tasks" / "ni_lkb_p_ni"
RUN_DEMO = WORKSPACE / "projects" / "p001" / "runs" / "run_demo"

TASK_RESULT_NAMES = frozenset({
    "ranking.csv",
    "report.json",
    "metrics.json",
    "evaluation_manifest.json",
    "heldout_ranking.csv",
    "recommend_manifest.json",
    "evaluate_model_manifest.json",
    "evaluate_models_manifest.json",
    "model_selection_summary.csv",
    "best_model.json",
    "model_selection_manifest.json",
})


def _task_dir_files(task_dir: Path) -> Set[str]:
    if not task_dir.is_dir():
        return set()
    return {p.name for p in task_dir.iterdir() if p.is_file()}

REPO_GUARD_PATTERNS = [
    CHEMDB_ROOT / "training" / "evaluation" / "results",
    CHEMDB_ROOT / "training" / "models",
    CHEMDB_ROOT / "result_analysis",
    CHEMDB_ROOT / "tmp",
    CHEMDB_ROOT / "src" / "tmp",
    WEB_ROOT / "ChemDB_restructured" / "ChemDB" / "training" / "evaluation" / "results",
]


def _snapshot_guard_dirs() -> Dict[str, Set[str]]:
    snap: Dict[str, Set[str]] = {}
    for base in REPO_GUARD_PATTERNS:
        key = str(base)
        if base.is_file():
            snap[key] = {base.name}
        elif base.exists():
            snap[key] = {str(p.relative_to(base)) for p in base.rglob("*") if p.is_file()}
        else:
            snap[key] = set()
    return snap


def _assert_guard_unchanged(before: Dict[str, Set[str]], after: Dict[str, Set[str]]) -> None:
    assert before == after, f"repo guard dirs changed: before={before}, after={after}"


def _write_tiny_task(task_dir: Path, *, with_negatives: bool = True, n_pos: int = 2) -> None:
    task_dir.mkdir(parents=True, exist_ok=True)
    spec = {
        "task_id": "tiny_task",
        "metal": "Fe",
        "embedding": "ecfp",
        "candidate_file": "candidates.csv",
        "positive_file": "positives.csv",
        "negative_file": "negatives.csv",
        "notes": "",
    }
    (task_dir / "task.json").write_text(json.dumps(spec, indent=2), encoding="utf-8")
    with open(task_dir / "candidates.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["candidate_id", "smiles", "did", "name", "source"])
        w.writerow(["c1", "C", "", "methane", "test"])
        w.writerow(["c2", "CC", "", "ethane", "test"])
        w.writerow(["c3", "CCC", "", "propane", "test"])
    with open(task_dir / "positives.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["positive_id", "smiles", "did", "name", "source"])
        for i in range(n_pos):
            w.writerow([f"p{i + 1}", "C" if i == 0 else "CC", "", f"pos{i}", "test"])
    if with_negatives:
        with open(task_dir / "negatives.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["negative_id", "smiles", "did", "name", "source"])
            w.writerow(["n1", "CCC", "", "neg", "test"])


def _seed_mock_run(run_root: Path) -> None:
    run_root.mkdir(parents=True, exist_ok=True)
    lig = run_root / "data" / "L3_embedding" / "L3_embedding_ecfp.npz"
    lig.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        lig,
        dids=np.array(["D1"], dtype=object),
        smiles=np.array(["C"], dtype=object),
        embeddings=np.zeros((1, 2048), dtype=np.float32),
    )
    metal = run_root / "data" / "metal_embedding" / "element_features_zscore.csv"
    metal.parent.mkdir(parents=True, exist_ok=True)
    metal.write_text("element,symbol,val\nFe,Fe,0.0\nPd,Pd,0.1\n", encoding="utf-8")
    train = run_root / "training"
    train.mkdir(parents=True, exist_ok=True)
    index = {
        "ligand_embedding": {"backend": "ecfp", "path": str(lig.resolve())},
        "metal_embedding": {"path": str(metal.resolve())},
    }
    (train / "index.json").write_text(json.dumps(index), encoding="utf-8")
    (train / "config.yaml").write_text(
        yaml.dump({"data": {"embedding": "ecfp"}, "model": {"activation": "gelu", "dropout": 0.1}}),
        encoding="utf-8",
    )
    ckpt = train / "ckpts" / "best.pt"
    ckpt.parent.mkdir(parents=True, exist_ok=True)
    ckpt.write_bytes(b"mock-ckpt")


def test_workspace_level_task_schema_validation() -> None:
    from training.model_after.io import load_task_bundle

    assert SHARED_TASK_NI.is_dir(), f"missing shared task fixture: {SHARED_TASK_NI}"
    bundle = load_task_bundle(SHARED_TASK_NI)
    assert bundle.task_id == "ni_lkb_p_ni"
    assert bundle.metal == "Ni"
    assert len(bundle.candidates) >= 2
    assert len(bundle.positives) >= 2
    assert len(bundle.negatives) == 0


def test_task_bundle_schema_validation(tmp_path: Path) -> None:
    from training.model_after.io import load_task_bundle

    task_dir = tmp_path / "task_ok"
    _write_tiny_task(task_dir)
    bundle = load_task_bundle(task_dir)
    assert bundle.task_id == "tiny_task"
    assert bundle.metal == "Fe"
    assert len(bundle.candidates) == 3
    assert len(bundle.positives) == 2

    bad = tmp_path / "task_bad"
    bad.mkdir()
    (bad / "task.json").write_text(json.dumps({"task_id": "x", "metal": ""}), encoding="utf-8")
    with pytest.raises(ValueError, match="metal"):
        load_task_bundle(bad)


def test_model_bundle_paths_under_run(tmp_path: Path) -> None:
    pytest.importorskip("torch")
    from training.model_after.bundle import load_run_bundle

    run_root = tmp_path / "mock_run"
    _seed_mock_run(run_root)
    fake_model = MagicMock()
    fake_model.eval = MagicMock()
    with patch("training.model_after.bundle.torch.load", return_value={"model_state_dict": {}, "config": {}}):
        with patch("training.model_after.bundle.model_module.RankModel", return_value=fake_model):
            with patch("training.model_after.bundle.load_ligand_embedding", return_value=({"D1": np.zeros(2048)}, 2048)):
                with patch("training.model_after.bundle.load_metal_embedding", return_value=({"Fe": np.zeros(1)}, 1)):
                    bundle = load_run_bundle(run_root, device="cpu")
    assert bundle.index_path == run_root / "training" / "index.json"
    assert bundle.ckpt_path == run_root / "training" / "ckpts" / "best.pt"
    assert str(run_root) in str(bundle.index["ligand_embedding"]["path"])


def test_recommend_single_model_mock(tmp_path: Path) -> None:
    from training.model_after import scoring

    task_dir = tmp_path / "task"
    out_dir = tmp_path / "out_rec"
    _write_tiny_task(task_dir)
    run_root = tmp_path / "run"
    _seed_mock_run(run_root)

    fake_report = {
        "status": "success",
        "context_source": "positives_mean",
        "ranking_path": str(out_dir / "ranking.csv"),
        "output_dir": str(out_dir),
    }

    with patch("training.model_after.scoring.load_run_bundle") as mock_load:
        mock_bundle = MagicMock()
        mock_bundle.model_run_root = run_root
        mock_bundle.model_id = "mock_run"
        mock_bundle.ckpt_path = run_root / "training" / "ckpts" / "best.pt"
        mock_bundle.device = "cpu"
        mock_load.return_value = mock_bundle

        def _fake_recommend(*args, **kwargs):
            out_dir.mkdir(parents=True, exist_ok=True)
            ranking = out_dir / "ranking.csv"
            with open(ranking, "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(
                    f,
                    fieldnames=[
                        "rank", "score", "candidate_id", "valid",
                        "is_known_positive", "is_known_negative",
                        "model_id", "model_run_root", "checkpoint_path",
                    ],
                )
                w.writeheader()
                w.writerow({
                    "rank": 1, "score": 0.9, "candidate_id": "c1", "valid": True,
                    "is_known_positive": True, "is_known_negative": False,
                    "model_id": "mock_run",
                    "model_run_root": str(run_root),
                    "checkpoint_path": str(mock_bundle.ckpt_path),
                })
            report = {**fake_report, "ranking_path": str(ranking)}
            (out_dir / "report.json").write_text(json.dumps(report), encoding="utf-8")
            return report

        with patch("training.model_after.scoring.recommend_single_model", side_effect=_fake_recommend):
            report = scoring.recommend_single_model(run_root, task_dir, out_dir)

    assert (out_dir / "ranking.csv").is_file()
    assert (out_dir / "report.json").is_file()
    assert report["context_source"] == "positives_mean"


def test_evaluate_single_model_mock(tmp_path: Path) -> None:
    from training.model_after import eval_task

    task_dir = tmp_path / "task"
    out_dir = tmp_path / "out_eval"
    _write_tiny_task(task_dir, n_pos=1)
    run_root = tmp_path / "run"
    _seed_mock_run(run_root)

    with patch("training.model_after.eval_task.recommend_single_model") as mock_rec:
        mock_rec.return_value = {
            "ranking_path": str(out_dir / "ranking.csv"),
        }
        out_dir.mkdir(parents=True, exist_ok=True)
        with open(out_dir / "ranking.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["candidate_id", "score", "valid"])
            w.writeheader()
            w.writerow({"candidate_id": "c1", "score": "0.5", "valid": "True"})

        with patch("training.model_after.eval_task.load_run_bundle") as mock_load:
            mock_bundle = MagicMock()
            mock_bundle.model_id = "mock"
            mock_bundle.model_run_root = run_root
            mock_bundle.ckpt_path = run_root / "training" / "ckpts" / "best.pt"
            mock_bundle.device = "cpu"
            mock_bundle.metal_lookup = {"Fe": np.zeros(4, dtype=np.float32)}
            mock_bundle.ligand_lookup = {}
            mock_bundle.d_l = 2048
            mock_load.return_value = mock_bundle
            with patch(
                "training.model_after.eval_task._heldout_metrics",
                return_value=(None, None, None, None, None, []),
            ):
                metrics = eval_task.evaluate_task_fit(run_root, task_dir, out_dir)

    assert metrics["status"] == "partial"
    assert metrics["mrr"] is None
    assert any("positives < 2" in w for w in metrics["warnings"])
    assert (out_dir / "metrics.json").is_file()


def test_evaluate_models_mock(tmp_path: Path) -> None:
    from training.model_after import batch

    task_dir = tmp_path / "task"
    _write_tiny_task(task_dir)
    out_dir = tmp_path / "multi_out"
    runs_csv = tmp_path / "model_runs.csv"
    run_a = tmp_path / "run_a"
    run_b = tmp_path / "run_b"
    _seed_mock_run(run_a)
    _seed_mock_run(run_b)
    with open(runs_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["model_id", "model_run_root", "checkpoint_path", "model_name", "notes"])
        w.writerow(["a", str(run_a), "", "A", ""])
        w.writerow(["b", str(run_b), "", "B", ""])

    def fake_eval(run_root, task, model_out, **kw):
        mid = kw.get("model_id") or Path(run_root).name
        model_out = Path(model_out)
        model_out.mkdir(parents=True, exist_ok=True)
        mrr = 0.8 if mid == "a" else 0.3
        metrics = {
            "status": "success",
            "mrr": mrr,
            "hit_at_5": mrr,
            "hit_at_10": mrr,
            "hit_at_20": mrr,
            "hit_at_50": mrr,
            "positive_score_mean": 0.1,
            "negative_score_mean": -0.1,
            "pos_neg_margin": 0.2,
            "warnings": [],
        }
        (model_out / "metrics.json").write_text(json.dumps(metrics), encoding="utf-8")
        (model_out / "ranking.csv").write_text("rank,score\n1,0.9\n", encoding="utf-8")
        return metrics

    with patch("training.model_after.batch.evaluate_task_fit", side_effect=fake_eval):
        manifest = batch.evaluate_models(runs_csv, task_dir, out_dir)

    assert (out_dir / "model_selection_summary.csv").is_file()
    assert (out_dir / "best_model.json").is_file()
    best = json.loads((out_dir / "best_model.json").read_text(encoding="utf-8"))
    assert best["model_id"] == "a"
    assert (out_dir / "models" / "a" / "metrics.json").is_file()
    assert manifest["n_success_with_score"] == 2


def test_no_write_repo_model_after(tmp_path: Path) -> None:
    from training.model_after.io import load_task_bundle

    before = _snapshot_guard_dirs()
    task_dir = tmp_path / "task"
    _write_tiny_task(task_dir)
    load_task_bundle(task_dir)
    after = _snapshot_guard_dirs()
    _assert_guard_unchanged(before, after)


def test_recommend_with_shared_task_cross_project(tmp_path: Path, monkeypatch) -> None:
    from app.services import model_after_manager as mam

    ws = tmp_path / "workspace"
    task_dir = ws / "model_after_tasks" / "ni_lkb_p_ni"
    shutil.copytree(SHARED_TASK_NI, task_dir)
    run_root = ws / "projects" / "p001" / "runs" / "run_demo"
    _seed_mock_run(run_root)
    out_dir = ws / "model_after_results" / "ni_lkb_p_ni" / "p001_run_demo"
    import app.services.model_after_paths as paths

    monkeypatch.setattr(paths, "WORKSPACE_ROOT", ws)
    monkeypatch.setattr(paths, "MODEL_AFTER_RESULTS_ROOT", ws / "model_after_results")

    before = _task_dir_files(task_dir)
    with patch("training.model_after.recommend_single_model") as mock_rec:
        def _fake(*args, **kwargs):
            out = Path(kwargs.get("output_dir") or args[2])
            out.mkdir(parents=True, exist_ok=True)
            (out / "ranking.csv").write_text("rank,score\n1,0.9\n", encoding="utf-8")
            return {"status": "success", "ranking_path": str(out / "ranking.csv")}

        mock_rec.side_effect = _fake
        rc = mam.main([
            "recommend",
            "--model-run-root", str(run_root),
            "--task-dir", str(task_dir),
            "--output-dir", str(out_dir),
            "--device", "cpu",
        ])
    assert rc == 0
    assert (out_dir / "ranking.csv").is_file()
    assert not (task_dir / "ranking.csv").exists()
    assert not (run_root / "reports").exists() or not list((run_root / "reports").rglob("ranking.csv"))
    assert _task_dir_files(task_dir) == before


def test_evaluate_models_cross_project(tmp_path: Path, monkeypatch) -> None:
    from training.model_after import batch
    import app.services.model_after_paths as paths

    ws = tmp_path / "workspace"
    task_dir = ws / "model_after_tasks" / "ni_lkb_p_ni"
    shutil.copytree(SHARED_TASK_NI, task_dir)
    run_a = ws / "projects" / "p001" / "runs" / "run_demo"
    run_b = ws / "projects" / "p002" / "runs" / "run_demo_2"
    _seed_mock_run(run_a)
    _seed_mock_run(run_b)
    out_dir = ws / "model_after_results" / "ni_lkb_p_ni" / "batch_001"
    runs_csv = tmp_path / "model_runs_cross_project.csv"
    with open(runs_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["model_id", "model_run_root", "checkpoint_path", "model_name", "notes"])
        w.writerow(["p001_demo", str(run_a), "", "p001", ""])
        w.writerow(["p002_demo2", str(run_b), "", "p002", ""])

    def fake_eval(run_root, task, model_out, **kw):
        model_out = Path(model_out)
        model_out.mkdir(parents=True, exist_ok=True)
        mrr = 0.9 if kw.get("model_id") == "p001_demo" else 0.2
        metrics = {
            "status": "success",
            "mrr": mrr,
            "hit_at_5": mrr,
            "hit_at_10": mrr,
            "hit_at_20": mrr,
            "hit_at_50": mrr,
            "positive_score_mean": 0.1,
            "negative_score_mean": -0.1,
            "pos_neg_margin": 0.2,
            "warnings": [],
        }
        (model_out / "metrics.json").write_text(json.dumps(metrics), encoding="utf-8")
        (model_out / "ranking.csv").write_text("rank,score\n1,0.9\n", encoding="utf-8")
        return metrics

    before = _task_dir_files(task_dir)
    with patch("training.model_after.batch.evaluate_task_fit", side_effect=fake_eval):
        manifest = batch.evaluate_models(runs_csv, task_dir, out_dir)

    assert (out_dir / "model_selection_summary.csv").is_file()
    best = json.loads((out_dir / "best_model.json").read_text(encoding="utf-8"))
    assert best["model_id"] == "p001_demo"
    assert (out_dir / "models" / "p001_demo" / "metrics.json").is_file()
    assert (out_dir / "models" / "p002_demo2" / "metrics.json").is_file()
    assert manifest["n_models"] == 2
    assert _task_dir_files(task_dir) == before


def test_task_dir_is_not_modified_by_default(tmp_path: Path, monkeypatch) -> None:
    from app.services import model_after_manager as mam
    import app.services.model_after_paths as paths

    ws = tmp_path / "workspace"
    task_dir = ws / "model_after_tasks" / "ni_lkb_p_ni"
    shutil.copytree(SHARED_TASK_NI, task_dir)
    run_root = ws / "projects" / "p001" / "runs" / "run_demo"
    _seed_mock_run(run_root)
    monkeypatch.setattr(paths, "WORKSPACE_ROOT", ws)
    monkeypatch.setattr(paths, "MODEL_AFTER_RESULTS_ROOT", ws / "model_after_results")

    before = _task_dir_files(task_dir)

    with patch("training.model_after.recommend_single_model") as mock_rec:
        mock_rec.return_value = {
            "status": "success",
            "ranking_path": str(ws / "model_after_results" / "ni_lkb_p_ni" / "p001_run_demo" / "ranking.csv"),
        }
        assert mam.main([
            "recommend",
            "--model-run-root", str(run_root),
            "--task-dir", str(task_dir),
            "--device", "cpu",
        ]) == 0

    with patch("training.model_after.evaluate_task_fit") as mock_ev:
        mock_ev.return_value = {"status": "partial", "mrr": None, "warnings": []}
        assert mam.main([
            "evaluate-model",
            "--model-run-root", str(run_root),
            "--task-dir", str(task_dir),
            "--device", "cpu",
        ]) == 0

    runs_csv = tmp_path / "model_runs_cross.csv"
    runs_csv.write_text(
        "model_id,model_run_root,checkpoint_path,model_name,notes\n"
        f"p001,{run_root},,,\n",
        encoding="utf-8",
    )
    with patch("training.model_after.evaluate_models") as mock_ms:
        mock_ms.return_value = {"summary_path": str(ws / "model_after_results" / "ni_lkb_p_ni" / "batch_default" / "model_selection_summary.csv")}
        assert mam.main([
            "evaluate-models",
            "--model-runs-file", str(runs_csv),
            "--task-dir", str(task_dir),
            "--device", "cpu",
        ]) == 0

    after = _task_dir_files(task_dir)
    assert after == before
    assert not any(name in after for name in TASK_RESULT_NAMES)
    assert (ws / "model_after_results" / "ni_lkb_p_ni").exists()


def test_outputs_under_workspace_model_after_results(tmp_path: Path) -> None:
    from training.model_after import batch

    ws = tmp_path / "workspace"
    task_dir = ws / "model_after_tasks" / "tiny"
    out_dir = ws / "model_after_results" / "tiny" / "batch_x"
    _write_tiny_task(task_dir)
    runs_csv = tmp_path / "runs.csv"
    run_root = tmp_path / "run_only"
    _seed_mock_run(run_root)
    with open(runs_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["model_id", "model_run_root", "checkpoint_path", "model_name", "notes"])
        w.writerow(["m1", str(run_root), "", "", ""])

    def _fake_eval(run_root, task, model_out, **kw):
        model_out = Path(model_out)
        model_out.mkdir(parents=True, exist_ok=True)
        (model_out / "metrics.json").write_text("{}", encoding="utf-8")
        return {"status": "partial", "mrr": None, "warnings": []}

    with patch("training.model_after.batch.evaluate_task_fit", side_effect=_fake_eval):
        batch.evaluate_models(runs_csv, task_dir, out_dir)

    assert str(out_dir) in str((out_dir / "model_selection_summary.csv").resolve())
    assert "model_after_results" in out_dir.as_posix()
    assert (out_dir / "models" / "m1" / "metrics.json").is_file()


def test_outputs_under_output_dir(tmp_path: Path) -> None:
    from training.model_after import batch

    task_dir = tmp_path / "task"
    _write_tiny_task(task_dir)
    out_dir = tmp_path / "isolated_out"
    runs_csv = tmp_path / "runs.csv"
    run_root = tmp_path / "run_only"
    _seed_mock_run(run_root)
    with open(runs_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["model_id", "model_run_root", "checkpoint_path", "model_name", "notes"])
        w.writerow(["m1", str(run_root), "", "", ""])

    def _fake_eval(run_root, task, model_out, **kw):
        model_out = Path(model_out)
        model_out.mkdir(parents=True, exist_ok=True)
        (model_out / "metrics.json").write_text("{}", encoding="utf-8")
        return {
            "status": "partial",
            "mrr": None,
            "hit_at_5": None,
            "hit_at_10": None,
            "hit_at_20": None,
            "hit_at_50": None,
            "positive_score_mean": None,
            "negative_score_mean": None,
            "pos_neg_margin": None,
            "warnings": ["positives < 2"],
        }

    with patch("training.model_after.batch.evaluate_task_fit", side_effect=_fake_eval):
        batch.evaluate_models(runs_csv, task_dir, out_dir)

    summary = out_dir / "model_selection_summary.csv"
    assert summary.is_file()
    assert str(out_dir) in str(summary.resolve())
    assert (out_dir / "models" / "m1" / "metrics.json").is_file()


def test_no_molclr_required_for_model_after() -> None:
    import importlib
    import training.model_after as ma

    importlib.import_module("training.model_after.bundle")
    importlib.import_module("training.model_after.scoring")
    assert "molclr_api" not in sys.modules

    pkg = Path(ma.__file__).resolve().parent
    for py in pkg.glob("*.py"):
        text = py.read_text(encoding="utf-8")
        assert "import molclr" not in text
        assert "from molclr" not in text
        assert "molclr_api" not in text


def test_two_task_outputs_do_not_overwrite(tmp_path: Path) -> None:
    from training.model_after import scoring

    def _write_out(out: Path, tag: str) -> None:
        out.mkdir(parents=True, exist_ok=True)
        (out / "ranking.csv").write_text(f"tag,{tag}\n", encoding="utf-8")
        (out / "report.json").write_text(json.dumps({"tag": tag}), encoding="utf-8")

    task1 = tmp_path / "task1"
    task2 = tmp_path / "task2"
    _write_tiny_task(task1)
    _write_tiny_task(task2)
    out1 = tmp_path / "out1"
    out2 = tmp_path / "out2"
    run_root = tmp_path / "run"

    with patch("training.model_after.scoring.recommend_single_model") as mock_fn:
        def side_effect(mrr, task_dir, output_dir, **kw):
            tag = Path(task_dir).name
            _write_out(Path(output_dir), tag)
            return {"ranking_path": str(Path(output_dir) / "ranking.csv"), "context_source": "positives_mean"}

        mock_fn.side_effect = side_effect
        scoring.recommend_single_model(run_root, task1, out1)
        scoring.recommend_single_model(run_root, task2, out2)

    assert "task1" in (out1 / "report.json").read_text(encoding="utf-8")
    assert "task2" in (out2 / "report.json").read_text(encoding="utf-8")
    assert out1 / "ranking.csv" != out2 / "ranking.csv"


@pytest.mark.integration_model_after_rank
def test_integration_model_after_rank(tmp_path: Path) -> None:
    pytest.importorskip("torch")
    pytest.importorskip("rdkit")
    if not (RUN_DEMO / "training" / "ckpts" / "best.pt").is_file():
        pytest.skip("run_demo best.pt not available")

    from training.model_after import recommend_single_model

    task_dir = tmp_path / "task"
    _write_tiny_task(task_dir)
    out_dir = tmp_path / "rank_out"
    report = recommend_single_model(RUN_DEMO, task_dir, out_dir, device="cpu")
    assert (out_dir / "ranking.csv").is_file()
    assert report["n_valid_candidates"] >= 1


@pytest.mark.integration_model_after_eval
def test_integration_model_after_eval(tmp_path: Path) -> None:
    pytest.importorskip("torch")
    pytest.importorskip("rdkit")
    if not (RUN_DEMO / "training" / "ckpts" / "best.pt").is_file():
        pytest.skip("run_demo best.pt not available")

    from training.model_after import evaluate_task_fit

    task_dir = tmp_path / "task"
    _write_tiny_task(task_dir, n_pos=2)
    out_dir = tmp_path / "eval_out"
    metrics = evaluate_task_fit(RUN_DEMO, task_dir, out_dir, device="cpu")
    assert (out_dir / "metrics.json").is_file()
    assert metrics["status"] in ("success", "partial")


@pytest.mark.integration_model_after_models
def test_integration_model_after_models(tmp_path: Path) -> None:
    pytest.importorskip("torch")
    pytest.importorskip("rdkit")
    if not (RUN_DEMO / "training" / "ckpts" / "best.pt").is_file():
        pytest.skip("run_demo best.pt not available")

    from training.model_after import evaluate_models

    run_b = tmp_path / "run_demo_copy"
    if run_b.exists():
        shutil.rmtree(run_b)
    shutil.copytree(RUN_DEMO, run_b)

    task_dir = tmp_path / "task"
    _write_tiny_task(task_dir, n_pos=2)
    out_dir = tmp_path / "selection"
    runs_csv = tmp_path / "model_runs.csv"
    with open(runs_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["model_id", "model_run_root", "checkpoint_path", "model_name", "notes"])
        w.writerow(["demo", str(RUN_DEMO), "", "demo", ""])
        w.writerow(["copy", str(run_b), "", "copy", ""])

    manifest = evaluate_models(runs_csv, task_dir, out_dir, device="cpu")
    assert (out_dir / "model_selection_summary.csv").is_file()
    assert manifest["n_models"] == 2


def test_model_after_manager_cli_recommend_mock(tmp_path: Path, capsys) -> None:
    from app.services import model_after_manager

    task_dir = tmp_path / "task"
    _write_tiny_task(task_dir)
    out_dir = tmp_path / "cli_out"
    run_root = tmp_path / "run"
    _seed_mock_run(run_root)

    with patch("training.model_after.recommend_single_model") as mock_rec:
        mock_rec.return_value = {
            "status": "success",
            "ranking_path": str(out_dir / "ranking.csv"),
            "context_source": "positives_mean",
        }
        rc = model_after_manager.main([
            "recommend",
            "--model-run-root", str(run_root),
            "--task-dir", str(task_dir),
            "--output-dir", str(out_dir),
            "--device", "cpu",
        ])
    assert rc == 0
    assert (out_dir / "recommend_manifest.json").is_file()
