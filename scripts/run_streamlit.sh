#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DB="${CHEMDB_DB_PATH:-$ROOT/workspace/chemdb.sqlite}"
WS="${CHEMDB_WORKSPACE_ROOT:-$ROOT/workspace}"

if [[ ! -f "$DB" ]]; then
  echo "warning: sqlite db not found at: $DB"
  echo "rebuild it first:"
  echo "  python -m app.storage.cli rebuild-index --workspace-root \"$WS\" --db-path \"$DB\""
fi

exec env PYTHONPATH=. streamlit run app/streamlit_ui/Home.py
