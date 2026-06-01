#!/usr/bin/env bash
# Phase 2A SQLite registry validation
# Run from ChemDBWebVersion project root.

set +e

DB="/tmp/chemdb_phase2a_validation.sqlite"
WORKSPACE="workspace"

echo "============================================================"
echo "0. Environment"
echo "============================================================"
pwd
python --version
sqlite3 --version
echo "DB=$DB"
echo "WORKSPACE=$WORKSPACE"

echo
echo "============================================================"
echo "1. Rebuild fresh temporary SQLite DB"
echo "============================================================"
rm -f "$DB" "$DB-wal" "$DB-shm"

python -m app.storage.cli rebuild-index \
  --workspace-root "$WORKSPACE" \
  --db-path "$DB"

echo
echo "DB file:"
ls -lh "$DB" "$DB-wal" "$DB-shm" 2>/dev/null

echo
echo "============================================================"
echo "2. Tables and schema"
echo "============================================================"
sqlite3 "$DB" ".tables"

echo
echo "---- schema: models ----"
sqlite3 "$DB" ".schema models"

echo
echo "---- schema: model_after_model_results ----"
sqlite3 "$DB" ".schema model_after_model_results"

echo
echo "---- schema: run_ga_bindings ----"
sqlite3 "$DB" ".schema run_ga_bindings"

echo
echo "---- PRAGMA table_info(models) ----"
sqlite3 -header -column "$DB" "PRAGMA table_info(models);"

echo
echo "---- PRAGMA table_info(model_after_model_results) ----"
sqlite3 -header -column "$DB" "PRAGMA table_info(model_after_model_results);"

echo
echo "============================================================"
echo "3. Row counts"
echo "============================================================"
sqlite3 -header -column "$DB" "
SELECT 'projects' AS table_name, COUNT(*) AS n FROM projects
UNION ALL SELECT 'runs', COUNT(*) FROM runs
UNION ALL SELECT 'ga_sets', COUNT(*) FROM ga_sets
UNION ALL SELECT 'ga_versions', COUNT(*) FROM ga_versions
UNION ALL SELECT 'run_ga_bindings', COUNT(*) FROM run_ga_bindings
UNION ALL SELECT 'step_executions', COUNT(*) FROM step_executions
UNION ALL SELECT 'artifacts', COUNT(*) FROM artifacts
UNION ALL SELECT 'models', COUNT(*) FROM models
UNION ALL SELECT 'model_after_tasks', COUNT(*) FROM model_after_tasks
UNION ALL SELECT 'model_after_batches', COUNT(*) FROM model_after_batches
UNION ALL SELECT 'model_after_model_results', COUNT(*) FROM model_after_model_results;
"

echo
echo "============================================================"
echo "4. Foreign key check"
echo "============================================================"
sqlite3 -header -column "$DB" "PRAGMA foreign_key_check;"

echo
echo "============================================================"
echo "5. Projects and runs"
echo "============================================================"
sqlite3 -header -column "$DB" "
SELECT *
FROM projects;
"

echo
sqlite3 -header -column "$DB" "
SELECT *
FROM runs;
"

echo
echo "============================================================"
echo "6. GA sets / versions / run bindings"
echo "============================================================"
sqlite3 -header -column "$DB" "
SELECT *
FROM ga_sets;
"

echo
sqlite3 -header -column "$DB" "
SELECT *
FROM ga_versions;
"

echo
sqlite3 -header -column "$DB" "
SELECT *
FROM run_ga_bindings;
"

echo
echo "---- run -> GA binding summary ----"
sqlite3 -header -column "$DB" "
SELECT
  r.project_id,
  r.run_id,
  b.status AS ga_binding_status,
  b.ga_set_id,
  b.ga_version_id,
  b.checksum
FROM runs r
LEFT JOIN run_ga_bindings b
  ON r.project_id = b.project_id
 AND r.run_id = b.run_id
ORDER BY r.project_id, r.run_id;
"

echo
echo "============================================================"
echo "7. Models: raw rows"
echo "============================================================"
sqlite3 -header -column "$DB" "
SELECT *
FROM models
LIMIT 20;
"

echo
echo "============================================================"
echo "8. Model-after tasks / batches"
echo "============================================================"
sqlite3 -header -column "$DB" "
SELECT *
FROM model_after_tasks
ORDER BY task_id;
"

echo
sqlite3 -header -column "$DB" "
SELECT *
FROM model_after_batches
ORDER BY task_id, batch_id;
"

echo
echo "============================================================"
echo "9. Model-after results: raw rows"
echo "============================================================"
sqlite3 -header -column "$DB" "
SELECT *
FROM model_after_model_results
ORDER BY task_id, batch_id, rank_among_models
LIMIT 50;
"

echo
echo "============================================================"
echo "10. Model-after results joined to models by raw model_id"
echo "============================================================"
sqlite3 -header -column "$DB" "
SELECT
  r.task_id,
  r.batch_id,
  r.model_id AS result_model_id,
  CASE WHEN m.model_id IS NULL THEN 'missing_model' ELSE 'matched_model' END AS model_join,
  m.model_id AS registry_model_id,
  m.project_id,
  m.run_id,
  r.status,
  r.mrr,
  r.hit_at_5,
  r.selection_score,
  r.rank_among_models
FROM model_after_model_results r
LEFT JOIN models m
  ON r.model_id = m.model_id
ORDER BY r.task_id, r.batch_id, r.rank_among_models;
"

echo
echo "============================================================"
echo "11. If registry_model_id exists, test join through it"
echo "============================================================"
sqlite3 -header -column "$DB" "
SELECT
  r.task_id,
  r.batch_id,
  r.model_id AS raw_model_id,
  r.registry_model_id,
  CASE WHEN m.model_id IS NULL THEN 'missing_model' ELSE 'matched_model' END AS model_join,
  m.project_id,
  m.run_id,
  r.status,
  r.mrr,
  r.hit_at_5,
  r.selection_score,
  r.rank_among_models
FROM model_after_model_results r
LEFT JOIN models m
  ON r.registry_model_id = m.model_id
ORDER BY r.task_id, r.batch_id, r.rank_among_models;
"

echo
echo "============================================================"
echo "12. Artifact samples"
echo "============================================================"
sqlite3 -header -column "$DB" "
SELECT *
FROM artifacts
ORDER BY project_id, run_id, step_id, role, rel_path
LIMIT 50;
"

echo
echo "============================================================"
echo "13. Step execution samples"
echo "============================================================"
sqlite3 -header -column "$DB" "
SELECT *
FROM step_executions
ORDER BY project_id, run_id, started_at
LIMIT 30;
"

echo
echo "============================================================"
echo "14. Idempotence check: counts before and after rebuild-index"
echo "============================================================"
echo "---- before ----"
sqlite3 -header -column "$DB" "
SELECT 'models' AS table_name, COUNT(*) AS n FROM models
UNION ALL SELECT 'results', COUNT(*) FROM model_after_model_results
UNION ALL SELECT 'artifacts', COUNT(*) FROM artifacts
UNION ALL SELECT 'steps', COUNT(*) FROM step_executions;
"

python -m app.storage.cli rebuild-index \
  --workspace-root "$WORKSPACE" \
  --db-path "$DB"

echo "---- after ----"
sqlite3 -header -column "$DB" "
SELECT 'models' AS table_name, COUNT(*) AS n FROM models
UNION ALL SELECT 'results', COUNT(*) FROM model_after_model_results
UNION ALL SELECT 'artifacts', COUNT(*) FROM artifacts
UNION ALL SELECT 'steps', COUNT(*) FROM step_executions;
"

echo
echo "============================================================"
echo "15. Repository query smoke test"
echo "============================================================"
python - <<'PY'
from pathlib import Path
from app.storage.db import connect_db
from app.storage import repositories as repo

db_path = Path("/tmp/chemdb_phase2a_validation.sqlite")
conn = connect_db(db_path)

def show(name, fn):
    print(f"\n---- {name} ----")
    try:
        out = fn()
        print(out)
    except Exception as e:
        print(type(e).__name__, str(e))

show("list_projects", lambda: repo.list_projects(conn))
show("list_runs", lambda: repo.list_runs(conn))
show("list_ga_sets", lambda: repo.list_ga_sets(conn))
show("list_models", lambda: repo.list_models(conn))
show("list_model_after_tasks", lambda: repo.list_model_after_tasks(conn))

try:
    tasks = repo.list_model_after_tasks(conn)
    for t in tasks:
        task_id = t["task_id"] if isinstance(t, dict) else t["task_id"]
        show(f"list_batches({task_id})", lambda task_id=task_id: repo.list_batches(conn, task_id))
except Exception as e:
    print("\n---- list_batches smoke skipped/error ----")
    print(type(e).__name__, str(e))
PY

echo
echo "============================================================"
echo "16. Check real workspace DB was not created/modified by this validation"
echo "============================================================"
ls -lh workspace/chemdb.sqlite workspace/chemdb.sqlite-wal workspace/chemdb.sqlite-shm 2>/dev/null || true

echo
echo "============================================================"
echo "VALIDATION DONE"
echo "Please paste the complete output back."
echo "============================================================"
