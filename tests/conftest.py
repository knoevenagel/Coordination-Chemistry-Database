"""Default pytest collection: scidata env required; integration excluded unless -m set."""

from __future__ import annotations

import os
import sys

import pytest

_SCIDATA_ENV = os.environ.get("CHEMDB_TEST_CONDA_ENV", "scidata")

_DEFAULT_MARKEXPR = (
    "not integration and not integration_1pct and not integration_1pct_tier_a "
    "and not integration_1pct_tier_b and not integration_1pct_tier_c "
    "and not integration_1pct_enhanced and not integration_embedding "
    "and not integration_training_data and not integration_training_train "
    "and not integration_model_after_rank and not integration_model_after_eval "
    "and not integration_model_after_models"
)


def _in_scidata_env() -> bool:
    if os.environ.get("CHEMDB_ALLOW_NON_SCIDATA_TEST", "").lower() in ("1", "true", "yes"):
        return True
    if os.environ.get("CONDA_DEFAULT_ENV") == _SCIDATA_ENV:
        return True
    prefix = os.environ.get("CONDA_PREFIX", "").rstrip("/\\")
    return prefix.endswith(f"/{_SCIDATA_ENV}") or prefix.endswith(f"\\{_SCIDATA_ENV}")


def pytest_configure(config) -> None:
    if not _in_scidata_env():
        pytest.exit(
            f"Tests must run in conda env '{_SCIDATA_ENV}' (RDKit/torch).\n"
            "  conda activate scidata && pytest ...\n"
            "  ./scripts/run_pytest.sh ...\n"
            "Override (not recommended): CHEMDB_ALLOW_NON_SCIDATA_TEST=1 pytest ...",
            returncode=2,
        )
    # When the user passes `-m ...`, do not apply the default exclusion filter.
    if config.option.markexpr:
        return
    config.option.markexpr = _DEFAULT_MARKEXPR
