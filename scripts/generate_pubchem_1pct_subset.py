#!/usr/bin/env python3
"""Generate a stratified 1% PubChem subset from ChemDB/data/pubchem (per-element CSV)."""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Sequence


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            block = f.read(chunk_size)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def derived_seed(global_seed: int, filename: str) -> int:
    return hash((global_seed, filename)) % (2**32)


def sample_indices(n_data: int, rate: float, min_rows: int, rng: random.Random) -> List[int]:
    if n_data <= 0:
        return []
    k = max(min_rows, round(n_data * rate))
    k = min(k, n_data)
    return sorted(rng.sample(range(n_data), k))


def process_csv(
    src: Path,
    dst: Path,
    rate: float,
    global_seed: int,
    min_rows: int,
    skip_source_checksum: bool,
) -> Dict[str, object]:
    with open(src, "r", encoding="utf-8-sig", newline="") as f:
        header_line = f.readline()
        data_lines = f.readlines()

    # Normalize header (BOM/leading spaces break step1 usecols=['cid', ...])
    header_cols = [c.strip() for c in header_line.strip().split(",")]
    header = ",".join(header_cols) + "\n"

    n_data = len(data_lines)
    rng = random.Random(derived_seed(global_seed, src.name))
    picked = set(sample_indices(n_data, rate, min_rows, rng))

    dst.parent.mkdir(parents=True, exist_ok=True)
    with open(dst, "w", encoding="utf-8-sig", newline="") as out:
        out.write(header)
        for i, line in enumerate(data_lines):
            if i in picked:
                out.write(line)

    entry: Dict[str, object] = {
        "name": src.name,
        "source_rows": n_data,
        "output_rows": len(picked),
        "source_bytes": src.stat().st_size,
        "output_bytes": dst.stat().st_size,
        "output_sha256": sha256_file(dst),
    }
    if not skip_source_checksum:
        entry["source_sha256"] = sha256_file(src)
    return entry


def write_readme(output_dir: Path, manifest: Dict[str, object]) -> None:
    agg = manifest["aggregate"]
    text = f"""# PubChem 1% integration subset

Generated from stratified per-element CSV sampling.

## Provenance

| Field | Value |
|-------|-------|
| Source | `{manifest["source_dir"]}` |
| Sample rate | {manifest["sample_rate"]} |
| Random seed | {manifest["random_seed"]} |
| Source files | {manifest["source_file_count"]} |
| Source data rows | {manifest["source_total_data_rows"]} |
| Output data rows | {manifest["output_total_data_rows"]} |
| Generated at | {manifest["generated_at"]} |

See `MANIFEST.json` for per-file row counts and SHA-256 checksums.

## Regenerate

```bash
cd ChemDBWebVersion
python scripts/generate_pubchem_1pct_subset.py \\
  --source ChemDB/data/pubchem \\
  --output integration_data/pubchem_1pct
```

## Run integration tests

```bash
export CHEMDB_INTEGRATION_PUBCHEM_SOURCE="$(pwd)/integration_data/pubchem_1pct"
pytest tests/run_isolation -m "integration_1pct_tier_a or integration_1pct_tier_b" -v
```

Optional Tier C: `-m integration_1pct_tier_c`

Enhanced (not Phase 1B minimum): `-m integration_1pct_enhanced`

Exclude enhanced when using parent marker on other suites:

```bash
pytest -m "integration_1pct and not integration_1pct_enhanced"
```

## Git

CSV files under this directory are **not** committed. Only this README and `MANIFEST.json` may be tracked.
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    web_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description="Generate stratified 1% PubChem subset")
    parser.add_argument(
        "--source",
        type=Path,
        default=web_root / "ChemDB" / "data" / "pubchem",
        help="Full PubChem directory (per-element *.csv)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=web_root / "integration_data" / "pubchem_1pct",
        help="Output directory for subset CSVs",
    )
    parser.add_argument("--rate", type=float, default=0.01, help="Sampling rate per file")
    parser.add_argument("--seed", type=int, default=42, help="Global random seed")
    parser.add_argument(
        "--min-rows-per-file",
        type=int,
        default=1,
        help="Minimum sampled rows per non-empty file",
    )
    parser.add_argument(
        "--skip-source-checksum",
        action="store_true",
        help="Skip SHA-256 of source files (faster; output_sha256 still computed)",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print stats only, do not write CSVs")
    args = parser.parse_args(argv)

    source_dir = args.source.resolve()
    output_dir = args.output.resolve()
    if not source_dir.is_dir():
        raise SystemExit(f"Source directory not found: {source_dir}")

    csv_files = sorted(source_dir.glob("*.csv"))
    if not csv_files:
        raise SystemExit(f"No *.csv in {source_dir}")

    file_entries: List[Dict[str, object]] = []
    source_total_rows = 0
    output_total_rows = 0
    source_total_bytes = 0
    output_total_bytes = 0

    for src in csv_files:
        with open(src, "r", encoding="utf-8-sig", newline="") as f:
            n_data = sum(1 for _ in f) - 1
        source_total_rows += max(n_data, 0)
        source_total_bytes += src.stat().st_size
        if args.dry_run:
            rng = random.Random(derived_seed(args.seed, src.name))
            k = len(sample_indices(n_data, args.rate, args.min_rows_per_file, rng))
            output_total_rows += k
            continue
        entry = process_csv(
            src,
            output_dir / src.name,
            args.rate,
            args.seed,
            args.min_rows_per_file,
            args.skip_source_checksum,
        )
        file_entries.append(entry)
        output_total_rows += int(entry["output_rows"])
        output_total_bytes += int(entry["output_bytes"])

    if args.dry_run:
        print(f"files={len(csv_files)} source_rows={source_total_rows} approx_output_rows={output_total_rows}")
        return 0

    manifest: Dict[str, object] = {
        "schema_version": "1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source_dir": str(source_dir),
        "output_dir": str(output_dir),
        "sample_rate": args.rate,
        "random_seed": args.seed,
        "sampling_method": "stratified_per_file",
        "min_rows_per_file": args.min_rows_per_file,
        "skip_source_checksum": args.skip_source_checksum,
        "source_file_count": len(csv_files),
        "source_total_data_rows": source_total_rows,
        "output_total_data_rows": output_total_rows,
        "files": file_entries,
        "aggregate": {
            "source_total_bytes": source_total_bytes,
            "output_total_bytes": output_total_bytes,
        },
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    write_readme(output_dir, manifest)
    print(f"Wrote {len(file_entries)} files to {output_dir}")
    print(f"Rows: {source_total_rows} -> {output_total_rows}")
    print(f"MANIFEST: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
