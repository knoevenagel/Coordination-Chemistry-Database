from pathlib import Path

from app.services.run_context import RunContext
from app.services.step_registry import STEP_REGISTRY, all_registered_paths_under_run


def test_step_registry_required_outputs_under_run_root(tmp_path: Path) -> None:
    ctx = RunContext.from_run_root(tmp_path / "run_registry")
    for spec in STEP_REGISTRY.values():
        for rel in all_registered_paths_under_run(spec):
            assert rel.startswith(("tmp/", "data/", "training/", "manifests/")), (
                f"{spec.step_id}: path {rel} must be under tmp/, data/, training/, or manifests/"
            )
            full = ctx.run_root / rel
            assert ctx.resolve_under_run(full), f"{spec.step_id}: {rel} escapes run_root"
