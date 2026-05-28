#!/usr/bin/env bash
# Canonical test runner: always uses conda env "scidata" (RDKit, torch, etc.).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

ENV_NAME="${CHEMDB_TEST_CONDA_ENV:-scidata}"

if ! command -v conda >/dev/null 2>&1; then
  echo "error: conda not found; activate scidata manually and run pytest" >&2
  exit 1
fi

if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "error: conda env '$ENV_NAME' not found (create or set CHEMDB_TEST_CONDA_ENV)" >&2
  exit 1
fi

exec conda run --no-capture-output -n "$ENV_NAME" python -m pytest "$@"
