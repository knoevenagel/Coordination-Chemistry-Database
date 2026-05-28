# PubChem 1% integration subset

Generated from stratified per-element CSV sampling.

## Provenance

| Field | Value |
|-------|-------|
| Source | `/home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB/data/pubchem` |
| Sample rate | 0.01 |
| Random seed | 42 |
| Source files | 79 |
| Source data rows | 1946214 |
| Output data rows | 19460 |
| Generated at | 2026-05-20T13:08:21.481768+00:00 |

See `MANIFEST.json` for per-file row counts and SHA-256 checksums.

## Regenerate

```bash
cd ChemDBWebVersion
python scripts/generate_pubchem_1pct_subset.py \
  --source ChemDB/data/pubchem \
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
