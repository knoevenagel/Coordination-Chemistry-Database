"""Run workspace paths for ChemDB pipeline isolation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class RunContext:
    run_root: Path
    data_dir: Path
    tmp_dir: Path
    training_dir: Path
    log_dir: Path
    report_dir: Path
    manifest_dir: Path

    @classmethod
    def from_run_root(cls, run_root: str | Path) -> RunContext:
        root = Path(run_root).resolve()
        return cls(
            run_root=root,
            data_dir=root / "data",
            tmp_dir=root / "tmp",
            training_dir=root / "training",
            log_dir=root / "logs",
            report_dir=root / "reports",
            manifest_dir=root / "manifests",
        )

    def create_dirs(self) -> None:
        for d in (
            self.data_dir,
            self.data_dir / "pubchem",
            self.data_dir / "metal_embedding",
            self.data_dir / "L3_embedding",
            self.tmp_dir,
            self.training_dir,
            self.training_dir / "ckpts",
            self.log_dir,
            self.report_dir,
            self.manifest_dir,
        ):
            d.mkdir(parents=True, exist_ok=True)

    def resolve_under_run(self, path: Path) -> bool:
        """True if path is inside run_root."""
        try:
            path.resolve().relative_to(self.run_root)
            return True
        except ValueError:
            return False
