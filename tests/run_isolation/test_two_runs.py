from pathlib import Path

from tests.run_isolation.conftest import seed_run_from_fixtures


def test_two_runs_do_not_overwrite(tmp_path: Path) -> None:
    run1 = tmp_path / "run_001"
    run2 = tmp_path / "run_002"
    ctx1 = seed_run_from_fixtures(run1)
    ctx2 = seed_run_from_fixtures(run2)

    f1 = ctx1.tmp_dir / "sentinel_run1.txt"
    f2 = ctx2.tmp_dir / "sentinel_run2.txt"
    f1.write_text("run1", encoding="utf-8")
    f2.write_text("run2", encoding="utf-8")

    assert f1.read_text(encoding="utf-8") == "run1"
    assert f2.read_text(encoding="utf-8") == "run2"
    assert not (ctx1.tmp_dir / "sentinel_run2.txt").exists()
    assert not (ctx2.tmp_dir / "sentinel_run1.txt").exists()
