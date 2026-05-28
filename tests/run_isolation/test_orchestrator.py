import json
from pathlib import Path
from app.services.orchestrator import init_run, run_step
from app.services.run_context import RunContext
from tests.run_isolation.conftest import FIXTURES, mock_run_subprocess, seed_run_from_fixtures


def test_init_run_creates_manifest_and_seed_files(tmp_path: Path) -> None:
    ctx = RunContext.from_run_root(tmp_path / "run_init")
    init_run(ctx)
    assert (ctx.data_dir / "metal_list.txt").is_file()
    assert (ctx.manifest_dir / "run.json").is_file()
    manifest = json.loads((ctx.manifest_dir / "run.json").read_text(encoding="utf-8"))
    assert manifest["run_root"] == str(ctx.run_root.resolve())
    assert "metal_list.txt" in manifest["seed_files_copied"]
    assert (ctx.data_dir / "metal_embedding" / "element_features.csv").is_file()
    assert "metal_embedding/element_features.csv" in manifest["seed_files_copied"]


def test_init_run_pubchem_source_file(tmp_path: Path) -> None:
    ctx = RunContext.from_run_root(tmp_path / "run_pubchem")
    src = FIXTURES / "minimal_pubchem.csv"
    init_run(ctx, pubchem_source=src)
    assert (ctx.data_dir / "pubchem" / "minimal_pubchem.csv").is_file()


def test_orchestrator_manifest_paths_under_run_root(tmp_path: Path) -> None:
    ctx = seed_run_from_fixtures(tmp_path / "run_manifest")
    ctx.tmp_dir.mkdir(parents=True, exist_ok=True)
    (ctx.tmp_dir / "ligand_data.csv").write_text(
        "ligand_did,ligand_smiles,source_complex_did\nD1,C,CX\n",
        encoding="utf-8",
    )

    def handler(cmd, log_f):
        log_f.write("ok\n")
        (ctx.tmp_dir / "repaired_ligand_data.csv").write_text(
            "ligand_new_did,source_did,new_smiles,old_smiles\nD1,D1,C,C\n",
            encoding="utf-8",
        )

    with mock_run_subprocess(handler):
        m = run_step(ctx, "step2")

    assert m["status"] == "success"
    assert m.get("paths_under_run_root") is True
    step_manifest = json.loads((ctx.manifest_dir / "step2.json").read_text(encoding="utf-8"))
    for key in ("command", "log_file", "run_root"):
        assert key in step_manifest
    assert str(ctx.run_root) in step_manifest["log_file"]
    assert (ctx.log_dir / "step2.log").read_text(encoding="utf-8").startswith("# command:")
