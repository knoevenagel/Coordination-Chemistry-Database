#!/usr/bin/env bash
# Validate Phase 2A against a bound-GA run using a copied workspace.
# Run from ChemDBWebVersion project root.

set +e

SRC_WORKSPACE="workspace"
TMP_WORKSPACE="/tmp/chemdb_workspace_bound_validation"
DB="/tmp/chemdb_phase2a_bound_validation.sqlite"

echo "============================================================"
echo "0. Prepare temporary copied workspace"
echo "============================================================"
rm -rf "$TMP_WORKSPACE"
rm -f "$DB" "$DB-wal" "$DB-shm"

cp -a "$SRC_WORKSPACE" "$TMP_WORKSPACE"

RUN_ROOT="$TMP_WORKSPACE/projects/p001/runs/run_demo"
GA_SET_ID="fixture_run_demo"
GA_VERSION_ID="v001"

echo "TMP_WORKSPACE=$TMP_WORKSPACE"
echo "RUN_ROOT=$RUN_ROOT"
echo "GA_SET_ID=$GA_SET_ID"
echo "GA_VERSION_ID=$GA_VERSION_ID"

echo
echo "============================================================"
echo "1. Bind workspace GA version to copied run"
echo "============================================================"

python -m app.services.ga_registry_manager bind-ga-version-to-run \
  --workspace-root "$TMP_WORKSPACE" \
  --run-root "$RUN_ROOT" \
  --ga-set-id "$GA_SET_ID" \
  --ga-version-id "$GA_VERSION_ID"

echo
echo "---- ga_binding.json ----"
cat "$RUN_ROOT/manifests/ga_binding.json" 2>/dev/null || true

echo
echo "---- materialized GA_with_id.csv ----"
ls -lh "$RUN_ROOT/tmp/GA_with_id.csv"
head -5 "$RUN_ROOT/tmp/GA_with_id.csv"

echo
echo "---- stale report, if created ----"
cat "$RUN_ROOT/manifests/ga_stale_report.json" 2>/dev/null || true

echo
echo "============================================================"
echo "2. Rebuild SQLite from copied workspace"
echo "============================================================"

python -m app.storage.cli rebuild-index \
  --workspace-root "$TMP_WORKSPACE" \
  --db-path "$DB"

echo
echo "============================================================"
echo "3. Check run GA binding status"
echo "============================================================"

sqlite3 -header -column "$DB" "
SELECT
  project_id,
  run_id,
  status,
  ga_set_id,
  ga_version_id,
  checksum,
  num_ga,
  materialized_path
FROM run_ga_bindings;
"

echo
echo "============================================================"
echo "4. Check runs table GA status"
echo "============================================================"

sqlite3 -header -column "$DB" "
SELECT
  project_id,
  run_id,
  ga_binding_status,
  status,
  pipeline_core_ok,
  pipeline_training_ok
FROM runs;
"

echo
echo "============================================================"
echo "5. Check models inherit GA metadata"
echo "============================================================"

sqlite3 -header -column "$DB" "
SELECT
  model_id,
  project_id,
  run_id,
  checkpoint_stem,
  ga_set_id,
  ga_version_id,
  ga_binding_checksum
FROM models
ORDER BY model_id;
"

echo
echo "============================================================"
echo "6. Check model-after result registry join still works"
echo "============================================================"

sqlite3 -header -column "$DB" "
SELECT
  r.task_id,
  r.batch_id,
  r.model_id AS raw_model_id,
  r.registry_model_id,
  CASE WHEN m.model_id IS NULL THEN 'missing_model' ELSE 'matched_model' END AS model_join,
  m.ga_set_id,
  m.ga_version_id,
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
echo "7. Idempotence check"
echo "============================================================"

echo "---- before ----"
sqlite3 -header -column "$DB" "
SELECT 'run_ga_bindings' AS table_name, COUNT(*) AS n FROM run_ga_bindings
UNION ALL SELECT 'models', COUNT(*) FROM models
UNION ALL SELECT 'results', COUNT(*) FROM model_after_model_results
UNION ALL SELECT 'artifacts', COUNT(*) FROM artifacts
UNION ALL SELECT 'steps', COUNT(*) FROM step_executions;
"

python -m app.storage.cli rebuild-index \
  --workspace-root "$TMP_WORKSPACE" \
  --db-path "$DB"

echo "---- after ----"
sqlite3 -header -column "$DB" "
SELECT 'run_ga_bindings' AS table_name, COUNT(*) AS n FROM run_ga_bindings
UNION ALL SELECT 'models', COUNT(*) FROM models
UNION ALL SELECT 'results', COUNT(*) FROM model_after_model_results
UNION ALL SELECT 'artifacts', COUNT(*) FROM artifacts
UNION ALL SELECT 'steps', COUNT(*) FROM step_executions;
"

echo
echo "============================================================"
echo "8. Confirm real workspace DB not created"
echo "============================================================"
ls -lh workspace/chemdb.sqlite workspace/chemdb.sqlite-wal workspace/chemdb.sqlite-shm 2>/dev/null || true

echo
echo "============================================================"
echo "BOUND-GA VALIDATION DONE"
echo "============================================================"
