"""Repo tmp guard is enforced via autouse fixture in conftest.py."""

from pathlib import Path

from app.services.orchestrator import chemdb_repo_root, run_step
from tests.run_isolation.conftest import seed_run_from_fixtures


def test_chemdb_repo_root_points_at_chem_db() -> None:
    root = chemdb_repo_root()
    assert root.name == "ChemDB"
    assert (root / "src" / "step1.py").is_file()


def test_run_step_does_not_touch_repo_tmp(tmp_path: Path) -> None:
    ctx = seed_run_from_fixtures(tmp_path / "run_step_mock")
    ctx.tmp_dir.mkdir(parents=True, exist_ok=True)
    (ctx.tmp_dir / "ligand_data.csv").write_text("ligand_new_did,ligand_smiles\n", encoding="utf-8")

    from tests.run_isolation.conftest import mock_run_subprocess

    def handler(cmd, log_f):
        log_f.write("mock ok\n")
        (ctx.tmp_dir / "repaired_ligand_data.csv").write_text(
            "ligand_new_did,source_did,new_smiles,old_smiles\n",
            encoding="utf-8",
        )

    with mock_run_subprocess(handler):
        manifest = run_step(ctx, "step2", skip_input_check=False)

    assert manifest["status"] == "success"
    assert (ctx.manifest_dir / "step2.json").is_file()
    assert (ctx.log_dir / "step2.log").is_file()
