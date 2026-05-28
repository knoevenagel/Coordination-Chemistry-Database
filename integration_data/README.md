# Integration data (Phase 1B)

This directory holds **generated** integration datasets. Large CSV files are **not** committed to git.

## PubChem 1% subset

Generate once from full `ChemDB/data/pubchem`:

```bash
cd ChemDBWebVersion
python scripts/generate_pubchem_1pct_subset.py \
  --source ChemDB/data/pubchem \
  --output integration_data/pubchem_1pct
```

Faster regeneration (skip source checksums):

```bash
python scripts/generate_pubchem_1pct_subset.py --skip-source-checksum
```

After generation, set:

```bash
export CHEMDB_INTEGRATION_PUBCHEM_SOURCE="$(pwd)/integration_data/pubchem_1pct"
```

Integration tests call `init-run --pubchem-source "$CHEMDB_INTEGRATION_PUBCHEM_SOURCE"`, which **copies** CSVs into `run_root/data/pubchem/`. Step scripts never read this directory directly.

## Pytest

**Use conda env `scidata`** (RDKit + torch). Canonical runner:

```bash
./scripts/run_pytest.sh    # from ChemDBWebVersion root
```

Direct `pytest` without `scidata` will exit with an environment error (`tests/conftest.py`).

| Command | What runs |
|---------|-----------|
| `./scripts/run_pytest.sh` | Default lightweight tests only (integration excluded via `tests/conftest.py`) |
| `pytest -m integration` | Minimal fixture, real step1+step2 |
| `pytest -m "integration_1pct_tier_a or integration_1pct_tier_b"` | Phase 1B minimum (1% data) |

Integration tests (via `tests/run_isolation/conftest.py`) set:

| Variable | Value | Purpose |
|----------|-------|---------|
| `CHEMDB_STEP1_SORT_BY_LENGTH` | `1` | Process long SMILES early (reduces 99% tqdm stall) |
| `CHEMDB_STEP1_ROW_TIMEOUT_SEC` | `600` | Skip a single row after 10 minutes |

Worker counts are **not** overridden in tests; step1/step2 use the same default as production (`min(cpu_count(), 200)`). Some structures legitimately need a long time on a worker.

**Observed runtime (scidata + RDKit, 1% data, default workers):** full Tier A+B run about **20–25 minutes** (`3 passed` in ~21 min). Step1 tqdm may sit at `19436/19437` for many minutes while the last hard rows finish; that is expected, not a pytest hang.

Use the **scidata** conda env (or any env with RDKit):

```bash
conda activate scidata
export CHEMDB_INTEGRATION_PUBCHEM_SOURCE="$(pwd)/integration_data/pubchem_1pct"
```

**See live step output during pytest** (logs also written under run `logs/`):

```bash
export CHEMDB_INTEGRATION_PUBCHEM_SOURCE="$(pwd)/integration_data/pubchem_1pct"
pytest -s tests/run_isolation/test_integration_1pct.py -v \
  -m "integration_1pct_tier_a or integration_1pct_tier_b"
```

`-s` disables pytest capture; `CHEMDB_STREAM_STEP_OUTPUT=1` is set automatically for integration markers. Manual runs: add `--stream` to `run-step` / `run-pipeline`.
| `pytest -m integration_1pct_tier_c` | Optional through step6_7 |
| `pytest -m integration_1pct_enhanced` | step8/12/13 (not Phase 1B minimum) |

**Do not** use bare `pytest -m integration_1pct` if any test is marked with both `integration_1pct` and `integration_1pct_enhanced`. Enhanced tests use **only** `integration_1pct_enhanced`. To run all non-enhanced 1% tiers:

```bash
pytest -m "integration_1pct and not integration_1pct_enhanced"
```

Or explicitly:

```bash
pytest -m "integration_1pct_tier_a or integration_1pct_tier_b or integration_1pct_tier_c"
```
