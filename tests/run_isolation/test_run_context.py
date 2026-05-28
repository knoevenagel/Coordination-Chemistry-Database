from pathlib import Path

from app.services.run_context import RunContext


def test_run_context_paths(tmp_path: Path) -> None:
    root = tmp_path / "runs" / "run_001"
    ctx = RunContext.from_run_root(root)
    ctx.create_dirs()

    assert ctx.run_root == root.resolve()
    assert ctx.data_dir == root.resolve() / "data"
    assert ctx.tmp_dir == root.resolve() / "tmp"
    assert ctx.training_dir == root.resolve() / "training"
    assert ctx.log_dir == root.resolve() / "logs"
    assert ctx.report_dir == root.resolve() / "reports"
    assert ctx.manifest_dir == root.resolve() / "manifests"
    assert ctx.data_dir.is_dir()
    assert (ctx.data_dir / "pubchem").is_dir()
