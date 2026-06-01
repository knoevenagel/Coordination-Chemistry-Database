# Phase 2A SQLite Registry 验证输出

| 字段 | 值 |
|------|-----|
| 日期 | 2026-05-28 |
| 脚本 | `scripts/phase2a_validation.sh` |
| DB | `/tmp/chemdb_phase2a_validation.sqlite` |
| Workspace | `workspace/` |

## 完整终端输出

```text
============================================================
0. Environment
============================================================
/home/chenghua/ZhuoLi/ChemDBWebVersion
Python 3.13.2
3.50.2 2025-06-28 14:00:48 2af157d77fb1304a74176eaee7fbc7c7e932d946bf25325e9c26c91db19e3079 (64-bit)
DB=/tmp/chemdb_phase2a_validation.sqlite
WORKSPACE=workspace

============================================================
1. Rebuild fresh temporary SQLite DB
============================================================
{
  "ok": true,
  "summary": {
    "projects": 1,
    "runs": 1,
    "ga_sets": 1,
    "tasks": 2,
    "batches": 3,
    "errors": []
  }
}

DB file:
-rw-r--r-- 1 root root 168K May 28 10:07 /tmp/chemdb_phase2a_validation.sqlite

============================================================
2. Tables and schema
============================================================
artifacts                  models                   
ga_sets                    projects                 
ga_versions                run_ga_bindings          
model_after_batches        runs                     
model_after_model_results  step_executions          
model_after_tasks        

---- schema: models ----
CREATE TABLE models (
    model_id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    checkpoint_path TEXT NOT NULL,
    checkpoint_stem TEXT,
    config_path TEXT,
    index_path TEXT,
    embedding_backend TEXT,
    l3_embedding_path TEXT,
    metal_embedding_path TEXT,
    history_path TEXT,
    best_epoch INTEGER,
    val_mrr REAL,
    val_recall_at_1 REAL,
    val_recall_at_5 REAL,
    ga_set_id TEXT,
    ga_version_id TEXT,
    ga_binding_checksum TEXT,
    size_bytes INTEGER,
    ingest_status TEXT NOT NULL DEFAULT 'ok',
    ingest_error TEXT,
    ingested_at TEXT NOT NULL DEFAULT (datetime('now')),
    FOREIGN KEY (project_id, run_id) REFERENCES runs(project_id, run_id) ON DELETE CASCADE,
    UNIQUE (project_id, run_id, checkpoint_path)
);
CREATE INDEX idx_models_run ON models(project_id, run_id);

---- schema: model_after_model_results ----
CREATE TABLE model_after_model_results (
    task_id TEXT NOT NULL,
    batch_id TEXT NOT NULL,
    model_id TEXT NOT NULL,
    model_run_root TEXT,
    checkpoint_path TEXT,
    project_id TEXT,
    run_id TEXT,
    registry_model_id TEXT,
    model_name TEXT,
    mrr REAL,
    hit_at_5 REAL,
    hit_at_10 REAL,
    hit_at_20 REAL,
    hit_at_50 REAL,
    selection_score REAL,
    rank_among_models INTEGER,
    status TEXT,
    ranking_path TEXT,
    metrics_path TEXT,
    report_path TEXT,
    ingest_status TEXT NOT NULL DEFAULT 'ok',
    ingest_error TEXT,
    ingested_at TEXT NOT NULL DEFAULT (datetime('now')),
    PRIMARY KEY (task_id, batch_id, model_id),
    FOREIGN KEY (task_id, batch_id) REFERENCES model_after_batches(task_id, batch_id) ON DELETE CASCADE
);
CREATE INDEX idx_model_results_task ON model_after_model_results(task_id, batch_id);

---- schema: run_ga_bindings ----
CREATE TABLE run_ga_bindings (
    project_id TEXT NOT NULL,
    run_id TEXT NOT NULL,
    ga_set_id TEXT,
    ga_version_id TEXT,
    checksum TEXT,
    num_ga INTEGER,
    source_ga_csv TEXT,
    materialized_path TEXT,
    bound_at TEXT,
    bound_by TEXT,
    status TEXT NOT NULL,
    is_current INTEGER NOT NULL DEFAULT 1,
    ingest_status TEXT NOT NULL DEFAULT 'ok',
    ingest_error TEXT,
    ingested_at TEXT NOT NULL DEFAULT (datetime('now')),
    PRIMARY KEY (project_id, run_id),
    FOREIGN KEY (project_id, run_id) REFERENCES runs(project_id, run_id) ON DELETE CASCADE,
    FOREIGN KEY (ga_set_id) REFERENCES ga_sets(ga_set_id),
    CHECK (is_current IN (0, 1))
);

---- PRAGMA table_info(models) ----
cid  name                  type     notnull  dflt_value       pk
---  --------------------  -------  -------  ---------------  --
0    model_id              TEXT     0                         1 
1    project_id            TEXT     1                         0 
2    run_id                TEXT     1                         0 
3    checkpoint_path       TEXT     1                         0 
4    checkpoint_stem       TEXT     0                         0 
5    config_path           TEXT     0                         0 
6    index_path            TEXT     0                         0 
7    embedding_backend     TEXT     0                         0 
8    l3_embedding_path     TEXT     0                         0 
9    metal_embedding_path  TEXT     0                         0 
10   history_path          TEXT     0                         0 
11   best_epoch            INTEGER  0                         0 
12   val_mrr               REAL     0                         0 
13   val_recall_at_1       REAL     0                         0 
14   val_recall_at_5       REAL     0                         0 
15   ga_set_id             TEXT     0                         0 
16   ga_version_id         TEXT     0                         0 
17   ga_binding_checksum   TEXT     0                         0 
18   size_bytes            INTEGER  0                         0 
19   ingest_status         TEXT     1        'ok'             0 
20   ingest_error          TEXT     0                         0 
21   ingested_at           TEXT     1        datetime('now')  0 

---- PRAGMA table_info(model_after_model_results) ----
cid  name               type     notnull  dflt_value       pk
---  -----------------  -------  -------  ---------------  --
0    task_id            TEXT     1                         1 
1    batch_id           TEXT     1                         2 
2    model_id           TEXT     1                         3 
3    model_run_root     TEXT     0                         0 
4    checkpoint_path    TEXT     0                         0 
5    project_id         TEXT     0                         0 
6    run_id             TEXT     0                         0 
7    registry_model_id  TEXT     0                         0 
8    model_name         TEXT     0                         0 
9    mrr                REAL     0                         0 
10   hit_at_5           REAL     0                         0 
11   hit_at_10          REAL     0                         0 
12   hit_at_20          REAL     0                         0 
13   hit_at_50          REAL     0                         0 
14   selection_score    REAL     0                         0 
15   rank_among_models  INTEGER  0                         0 
16   status             TEXT     0                         0 
17   ranking_path       TEXT     0                         0 
18   metrics_path       TEXT     0                         0 
19   report_path        TEXT     0                         0 
20   ingest_status      TEXT     1        'ok'             0 
21   ingest_error       TEXT     0                         0 
22   ingested_at        TEXT     1        datetime('now')  0 

============================================================
3. Row counts
============================================================
table_name                 n 
-------------------------  --
projects                   1 
runs                       1 
ga_sets                    1 
ga_versions                1 
run_ga_bindings            1 
step_executions            13
artifacts                  80
models                     2 
model_after_tasks          2 
model_after_batches        3 
model_after_model_results  4 

============================================================
4. Foreign key check
============================================================

============================================================
5. Projects and runs
============================================================
project_id  project_root                                                    name  created_at  updated_at  notes  ingest_status  ingest_error  ingested_at                     
----------  --------------------------------------------------------------  ----  ----------  ----------  -----  -------------  ------------  --------------------------------
p001        /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001  p001                                 ok                           2026-05-28T10:07:26.591744+00:00

project_id  run_id    run_root                                                                      created_at                        chemdb_repo_root                               status   pubchem_file_count  pipeline_core_ok  pipeline_training_ok  last_finished_at                  ga_binding_status  ingest_status  ingest_error  ingested_at                     
----------  --------  ----------------------------------------------------------------------------  --------------------------------  ---------------------------------------------  -------  ------------------  ----------------  --------------------  --------------------------------  -----------------  -------------  ------------  --------------------------------
p001        run_demo  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo  2026-05-21T06:13:47.371952+00:00  /home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB  success  79                  1                 1                     2026-05-21T06:43:08.953414+00:00  legacy_auto_ga     ok                           2026-05-28T10:07:26.592137+00:00

============================================================
6. GA sets / versions / run bindings
============================================================
ga_set_id         ga_set_root                                                                name                 description     created_at                 updated_at                 tags_json    source_json          ingest_status  ingest_error  ingested_at                     
----------------  -------------------------------------------------------------------------  -------------------  --------------  -------------------------  -------------------------  -----------  -------------------  -------------  ------------  --------------------------------
fixture_run_demo  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/ga_sets/fixture_run_demo  run_demo GA (small)  pytest fixture  2026-05-21T00:00:00+00:00  2026-05-21T00:00:00+00:00  ["fixture"]  {"type": "fixture"}  ok                           2026-05-28T10:07:26.607906+00:00

ga_set_id         ga_version_id  parent_version_id  created_at                 created_by  num_ga  checksum                                                                 ga_csv_path                                                                                             source_json          notes  ingest_status  ingest_error  ingested_at                     
----------------  -------------  -----------------  -------------------------  ----------  ------  -----------------------------------------------------------------------  ------------------------------------------------------------------------------------------------------  -------------------  -----  -------------  ------------  --------------------------------
fixture_run_demo  v001                              2026-05-21T00:00:00+00:00  fixture     126     sha256:edc18600d0819d62511afa6de4a410c27b0830f7a0d6b8f92f16df99921e56ce  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/ga_sets/fixture_run_demo/versions/v001/GA_with_id.csv  {"type": "fixture"}         ok                           2026-05-28T10:07:26.608098+00:00

project_id  run_id    ga_set_id  ga_version_id  checksum  num_ga  source_ga_csv  materialized_path  bound_at  bound_by  status          is_current  ingest_status  ingest_error  ingested_at                     
----------  --------  ---------  -------------  --------  ------  -------------  -----------------  --------  --------  --------------  ----------  -------------  ------------  --------------------------------
p001        run_demo                                                                                                    legacy_auto_ga  1           ok                           2026-05-28T10:07:26.592515+00:00

---- run -> GA binding summary ----
project_id  run_id    ga_binding_status  ga_set_id  ga_version_id  checksum
----------  --------  -----------------  ---------  -------------  --------
p001        run_demo  legacy_auto_ga                                       

============================================================
7. Models: raw rows
============================================================
model_id            project_id  run_id    checkpoint_path                                                                                      checkpoint_stem  config_path                                                                                        index_path                                                                                        embedding_backend  l3_embedding_path                                                                                                     metal_embedding_path                                                                                                           history_path                                                                                              best_epoch  val_mrr             val_recall_at_1  val_recall_at_5    ga_set_id  ga_version_id  ga_binding_checksum  size_bytes  ingest_status  ingest_error  ingested_at                     
------------------  ----------  --------  ---------------------------------------------------------------------------------------------------  ---------------  -------------------------------------------------------------------------------------------------  ------------------------------------------------------------------------------------------------  -----------------  --------------------------------------------------------------------------------------------------------------------  -----------------------------------------------------------------------------------------------------------------------------  --------------------------------------------------------------------------------------------------------  ----------  ------------------  ---------------  -----------------  ---------  -------------  -------------------  ----------  -------------  ------------  --------------------------------
p001_run_demo_best  p001        run_demo  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/ckpts/best.pt  best             /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/config.yaml  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/index.json  ecfp               /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/L3_embedding/L3_embedding_ecfp.npz  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/metal_embedding/element_features_zscore.csv  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/ckpts/history.json  3           0.0760487309346596  0.0              0.133333333333333                                                 236892252   ok                           2026-05-28T10:07:26.607211+00:00
p001_run_demo_last  p001        run_demo  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/ckpts/last.pt  last             /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/config.yaml  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/index.json  ecfp               /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/L3_embedding/L3_embedding_ecfp.npz  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/metal_embedding/element_features_zscore.csv  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/ckpts/history.json  3           0.0760487309346596  0.0              0.133333333333333                                                 236892252   ok                           2026-05-28T10:07:26.607452+00:00

============================================================
8. Model-after tasks / batches
============================================================
task_id           task_dir                                                                             metal  embedding  candidate_file  positive_file  negative_file  notes                                                                                                                   candidate_count  positive_count  negative_count  ingest_status  ingest_error  ingested_at                     
----------------  -----------------------------------------------------------------------------------  -----  ---------  --------------  -------------  -------------  ----------------------------------------------------------------------------------------------------------------------  ---------------  --------------  --------------  -------------  ------------  --------------------------------
ni_lkb_p_ni       /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_tasks/ni_lkb_p_ni       Ni     ecfp       candidates.csv  positives.csv  negatives.csv  Converted from ChemDB_restructured legacy quest ni_lkb_p_ni.csv (344 candidates, 11 positives, no explicit negatives).  344              11              0               ok                           2026-05-28T10:07:26.611802+00:00
pd_lkb_p_cluster  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_tasks/pd_lkb_p_cluster  Pd     ecfp       candidates.csv  positives.csv  negatives.csv  Converted from ChemDB_restructured legacy quest pd_lkb_p_cluster.csv                                                    343              10              1               ok                           2026-05-28T10:07:26.613054+00:00

task_id           batch_id           output_dir                                                                                           command          model_runs_file                                                                                                    finished_at                       n_models  n_success  best_model_id  summary_path                                                                                                                    ingest_status  ingest_error  ingested_at                     
----------------  -----------------  ---------------------------------------------------------------------------------------------------  ---------------  -----------------------------------------------------------------------------------------------------------------  --------------------------------  --------  ---------  -------------  ------------------------------------------------------------------------------------------------------------------------------  -------------  ------------  --------------------------------
ni_lkb_p_ni       batch_dummy_ckpts  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/ni_lkb_p_ni/batch_dummy_ckpts   evaluate-models  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/ni_lkb_p_ni/batch_dummy_ckpts/model_runs.csv  2026-05-21T09:17:57.873556+00:00  2         2          run_demo_last  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/ni_lkb_p_ni/batch_dummy_ckpts/model_selection_summary.csv  ok                           2026-05-28T10:07:26.610452+00:00
ni_lkb_p_ni       p001_run_demo      /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/ni_lkb_p_ni/p001_run_demo       evaluate-model                                                                                                                      2026-05-21T09:05:05.892604+00:00                                                                                                                                                                      ok                           2026-05-28T10:07:26.612153+00:00
pd_lkb_p_cluster  p001_run_demo      /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/pd_lkb_p_cluster/p001_run_demo  evaluate-model                                                                                                                      2026-05-21T09:04:32.709394+00:00                                                                                                                                                                      ok                           2026-05-28T10:07:26.613409+00:00

============================================================
9. Model-after results: raw rows
============================================================
task_id           batch_id           model_id       model_run_root                                                                checkpoint_path                                                                                      project_id  run_id    registry_model_id   model_name       mrr                 hit_at_5           hit_at_10          hit_at_20             hit_at_50          selection_score     rank_among_models  status   ranking_path                                                                                                                         metrics_path                                                                                                                          report_path                                                                                                                          ingest_status  ingest_error  ingested_at                     
----------------  -----------------  -------------  ----------------------------------------------------------------------------  ---------------------------------------------------------------------------------------------------  ----------  --------  ------------------  ---------------  ------------------  -----------------  -----------------  --------------------  -----------------  ------------------  -----------------  -------  -----------------------------------------------------------------------------------------------------------------------------------  ------------------------------------------------------------------------------------------------------------------------------------  -----------------------------------------------------------------------------------------------------------------------------------  -------------  ------------  --------------------------------
ni_lkb_p_ni       batch_dummy_ckpts  run_demo_last  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/ckpts/last.pt  p001        run_demo  p001_run_demo_last  last checkpoint  0.171211730938318   0.187195689615044  0.272575214510698  0.272900277738987     0.44580947726109   0.171211730938318   1                  success  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/ni_lkb_p_ni/batch_dummy_ckpts/models/run_demo_last/ranking.csv  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/ni_lkb_p_ni/batch_dummy_ckpts/models/run_demo_last/metrics.json  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/ni_lkb_p_ni/batch_dummy_ckpts/models/run_demo_last/report.json  ok                           2026-05-28T10:07:26.611098+00:00
ni_lkb_p_ni       batch_dummy_ckpts  run_demo_best  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/ckpts/best.pt  p001        run_demo  p001_run_demo_best  best checkpoint  0.0155332868367695  0.0                0.0                0.00210651057425251   0.262319236916011  0.0155332868367695  2                  success  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/ni_lkb_p_ni/batch_dummy_ckpts/models/run_demo_best/ranking.csv  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/ni_lkb_p_ni/batch_dummy_ckpts/models/run_demo_best/metrics.json  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/ni_lkb_p_ni/batch_dummy_ckpts/models/run_demo_best/report.json  ok                           2026-05-28T10:07:26.610837+00:00
ni_lkb_p_ni       p001_run_demo      run_demo       /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/ckpts/best.pt  p001        run_demo  p001_run_demo_best                   0.0155332868367695  0.0                0.0                0.00210651057425251   0.262319236916011                                         success  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/ni_lkb_p_ni/p001_run_demo/ranking.csv                           /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/ni_lkb_p_ni/p001_run_demo/metrics.json                           /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/ni_lkb_p_ni/p001_run_demo/report.json                           ok                           2026-05-28T10:07:26.612375+00:00
pd_lkb_p_cluster  p001_run_demo      run_demo       /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/ckpts/best.pt  p001        run_demo  p001_run_demo_best                   0.0159977312363552  0.0                0.0                0.000231028484453142  0.30264537321778                                          success  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/pd_lkb_p_cluster/p001_run_demo/ranking.csv                      /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/pd_lkb_p_cluster/p001_run_demo/metrics.json                      /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/pd_lkb_p_cluster/p001_run_demo/report.json                      ok                           2026-05-28T10:07:26.613636+00:00

============================================================
10. Model-after results joined to models by raw model_id
============================================================
task_id           batch_id           result_model_id  model_join     registry_model_id  project_id  run_id  status   mrr                 hit_at_5           selection_score     rank_among_models
----------------  -----------------  ---------------  -------------  -----------------  ----------  ------  -------  ------------------  -----------------  ------------------  -----------------
ni_lkb_p_ni       batch_dummy_ckpts  run_demo_last    missing_model                                         success  0.171211730938318   0.187195689615044  0.171211730938318   1                
ni_lkb_p_ni       batch_dummy_ckpts  run_demo_best    missing_model                                         success  0.0155332868367695  0.0                0.0155332868367695  2                
ni_lkb_p_ni       p001_run_demo      run_demo         missing_model                                         success  0.0155332868367695  0.0                                                     
pd_lkb_p_cluster  p001_run_demo      run_demo         missing_model                                         success  0.0159977312363552  0.0                                                     

============================================================
11. If registry_model_id exists, test join through it
============================================================
task_id           batch_id           raw_model_id   registry_model_id   model_join     project_id  run_id    status   mrr                 hit_at_5           selection_score     rank_among_models
----------------  -----------------  -------------  ------------------  -------------  ----------  --------  -------  ------------------  -----------------  ------------------  -----------------
ni_lkb_p_ni       batch_dummy_ckpts  run_demo_last  p001_run_demo_last  matched_model  p001        run_demo  success  0.171211730938318   0.187195689615044  0.171211730938318   1                
ni_lkb_p_ni       batch_dummy_ckpts  run_demo_best  p001_run_demo_best  matched_model  p001        run_demo  success  0.0155332868367695  0.0                0.0155332868367695  2                
ni_lkb_p_ni       p001_run_demo      run_demo       p001_run_demo_best  matched_model  p001        run_demo  success  0.0155332868367695  0.0                                                     
pd_lkb_p_cluster  p001_run_demo      run_demo       p001_run_demo_best  matched_model  p001        run_demo  success  0.0159977312363552  0.0                                                     

============================================================
12. Artifact samples
============================================================
project_id  run_id    step_id                  rel_path                                          role    abs_path                                                                                                                       size_bytes  exists_flag  artifact_type  ingest_status  ingest_error  ingested_at                     
----------  --------  -----------------------  ------------------------------------------------  ------  -----------------------------------------------------------------------------------------------------------------------------  ----------  -----------  -------------  -------------  ------------  --------------------------------
p001        run_demo                           data/L3_embedding/L3_embedding_ecfp.npz           key     /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/L3_embedding/L3_embedding_ecfp.npz           220359965   1            npz            ok                           2026-05-28T10:07:26.604736+00:00
p001        run_demo                           tmp/GA_with_id.csv                                key     /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/GA_with_id.csv                                2352        1            csv            ok                           2026-05-28T10:07:26.604514+00:00
p001        run_demo                           tmp/step13_kl_nl_samples.csv                      key     /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/step13_kl_nl_samples.csv                      1894037     1            csv            ok                           2026-05-28T10:07:26.604620+00:00
p001        run_demo                           training/ckpts/best.pt                            key     /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/ckpts/best.pt                            236892252   1            pt             ok                           2026-05-28T10:07:26.604165+00:00
p001        run_demo                           training/ckpts/history.json                       key     /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/ckpts/history.json                       567         1            json           ok                           2026-05-28T10:07:26.604406+00:00
p001        run_demo                           training/ckpts/last.pt                            key     /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/ckpts/last.pt                            236892252   1            pt             ok                           2026-05-28T10:07:26.604281+00:00
p001        run_demo                           training/config.yaml                              key     /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/config.yaml                              796         1            yaml           ok                           2026-05-28T10:07:26.604045+00:00
p001        run_demo                           training/index.json                               key     /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/index.json                               1584        1            json           ok                           2026-05-28T10:07:26.603829+00:00
p001        run_demo  build_l3_embedding_ecfp  tmp/metal_l3_index.csv                            input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/metal_l3_index.csv                            34164398    1            csv            ok                           2026-05-28T10:07:26.593639+00:00
p001        run_demo  build_l3_embedding_ecfp  tmp/repaired_ligand_data.csv                      input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/repaired_ligand_data.csv                      5011954     1            csv            ok                           2026-05-28T10:07:26.593438+00:00
p001        run_demo  build_l3_embedding_ecfp  data/L3_embedding/L3_embedding_ecfp.npz           output  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/L3_embedding/L3_embedding_ecfp.npz           220359965   1            npz            ok                           2026-05-28T10:07:26.593758+00:00
p001        run_demo  prepare_training_config  training/index.json                               input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/index.json                               1584        1            json           ok                           2026-05-28T10:07:26.593964+00:00
p001        run_demo  prepare_training_config  training/config.yaml                              output  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/config.yaml                              796         1            yaml           ok                           2026-05-28T10:07:26.594064+00:00
p001        run_demo  prepare_training_index   data/L3_embedding/L3_embedding_ecfp.npz           input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/L3_embedding/L3_embedding_ecfp.npz           220359965   1            npz            ok                           2026-05-28T10:07:26.594561+00:00
p001        run_demo  prepare_training_index   data/metal_embedding/element_features_zscore.csv  input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/metal_embedding/element_features_zscore.csv  93702       1            csv            ok                           2026-05-28T10:07:26.594669+00:00
p001        run_demo  prepare_training_index   tmp/l3_gac.json                                   input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/l3_gac.json                                   815892      1            json           ok                           2026-05-28T10:07:26.594456+00:00
p001        run_demo  prepare_training_index   tmp/m_l3_pairs.csv                                input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/m_l3_pairs.csv                                647764      1            csv            ok                           2026-05-28T10:07:26.594357+00:00
p001        run_demo  prepare_training_index   tmp/step13_kl_nl_samples.csv                      input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/step13_kl_nl_samples.csv                      1894037     1            csv            ok                           2026-05-28T10:07:26.594256+00:00
p001        run_demo  prepare_training_index   training/index.json                               output  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/index.json                               1584        1            json           ok                           2026-05-28T10:07:26.594765+00:00
p001        run_demo  step1                    data/metal_list.txt                               input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/metal_list.txt                               258         1            other          ok                           2026-05-28T10:07:26.595041+00:00
p001        run_demo  step1                    data/pubchem                                      input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/pubchem                                                  1            other          ok                           2026-05-28T10:07:26.594945+00:00
p001        run_demo  step1                    tmp/complex_data.csv                              output  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/complex_data.csv                              3211377     1            csv            ok                           2026-05-28T10:07:26.595142+00:00
p001        run_demo  step1                    tmp/ligand_data.csv                               output  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/ligand_data.csv                               3759783     1            csv            ok                           2026-05-28T10:07:26.595239+00:00
p001        run_demo  step12                   tmp/ligand_with_gac.csv                           input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/ligand_with_gac.csv                           1727303     1            csv            ok                           2026-05-28T10:07:26.595612+00:00
p001        run_demo  step12                   tmp/neo4j_l1_l3_relationships.csv                 input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/neo4j_l1_l3_relationships.csv                 2036177     1            csv            ok                           2026-05-28T10:07:26.595824+00:00
p001        run_demo  step12                   tmp/neo4j_l3_l4_relationships.csv                 input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/neo4j_l3_l4_relationships.csv                 11339633    1            csv            ok                           2026-05-28T10:07:26.595416+00:00
p001        run_demo  step12                   tmp/neo4j_l4_l5_relationships.csv                 input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/neo4j_l4_l5_relationships.csv                 2324668     1            csv            ok                           2026-05-28T10:07:26.595513+00:00
p001        run_demo  step12                   tmp/neo4j_m_l1_relationships.csv                  input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/neo4j_m_l1_relationships.csv                  680389      1            csv            ok                           2026-05-28T10:07:26.595710+00:00
p001        run_demo  step12                   tmp/neo4j_repaired_ligands.csv                    input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/neo4j_repaired_ligands.csv                    2189836     1            csv            ok                           2026-05-28T10:07:26.595921+00:00
p001        run_demo  step12                   tmp/l3_gac.json                                   output  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/l3_gac.json                                   815892      1            json           ok                           2026-05-28T10:07:26.596304+00:00
p001        run_demo  step12                   tmp/l3_l5.json                                    output  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/l3_l5.json                                    4860369     1            json           ok                           2026-05-28T10:07:26.596019+00:00
p001        run_demo  step12                   tmp/l5_freq_weight.json                           output  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/l5_freq_weight.json                           19928       1            json           ok                           2026-05-28T10:07:26.596209+00:00
p001        run_demo  step12                   tmp/l5_l3.json                                    output  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/l5_l3.json                                    4045999     1            json           ok                           2026-05-28T10:07:26.596114+00:00
p001        run_demo  step12                   tmp/l5_l3_filtered_K30.json                       output  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/l5_l3_filtered_K30.json                       3651772     1            json           ok                           2026-05-28T10:07:26.596404+00:00
p001        run_demo  step12                   tmp/m_l3_pairs.csv                                output  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/m_l3_pairs.csv                                647764      1            csv            ok                           2026-05-28T10:07:26.596501+00:00
p001        run_demo  step13                   tmp/l3_gac.json                                   input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/l3_gac.json                                   815892      1            json           ok                           2026-05-28T10:07:26.597096+00:00
p001        run_demo  step13                   tmp/l3_l5.json                                    input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/l3_l5.json                                    4860369     1            json           ok                           2026-05-28T10:07:26.596807+00:00
p001        run_demo  step13                   tmp/l5_freq_weight.json                           input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/l5_freq_weight.json                           19928       1            json           ok                           2026-05-28T10:07:26.597001+00:00
p001        run_demo  step13                   tmp/l5_l3.json                                    input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/l5_l3.json                                    4045999     1            json           ok                           2026-05-28T10:07:26.596904+00:00
p001        run_demo  step13                   tmp/m_l3_pairs.csv                                input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/m_l3_pairs.csv                                647764      1            csv            ok                           2026-05-28T10:07:26.596676+00:00
p001        run_demo  step13                   tmp/step13_kl_nl_samples.csv                      output  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/step13_kl_nl_samples.csv                      1894037     1            csv            ok                           2026-05-28T10:07:26.597196+00:00
p001        run_demo  step2                    tmp/ligand_data.csv                               input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/ligand_data.csv                               3759783     1            csv            ok                           2026-05-28T10:07:26.597360+00:00
p001        run_demo  step2                    tmp/repaired_ligand_data.csv                      output  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/repaired_ligand_data.csv                      5011954     1            csv            ok                           2026-05-28T10:07:26.597458+00:00
p001        run_demo  step4_5                  tmp/repaired_ligand_data.csv                      input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/repaired_ligand_data.csv                      5011954     1            csv            ok                           2026-05-28T10:07:26.597621+00:00
p001        run_demo  step4_5                  tmp/GA_with_id.csv                                output  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/GA_with_id.csv                                2352        1            csv            ok                           2026-05-28T10:07:26.597721+00:00
p001        run_demo  step4_5                  tmp/IRL_filtered.csv                              output  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/IRL_filtered.csv                              7691        1            csv            ok                           2026-05-28T10:07:26.597912+00:00
p001        run_demo  step4_5                  tmp/ligand_with_gac.csv                           output  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/ligand_with_gac.csv                           1727303     1            csv            ok                           2026-05-28T10:07:26.597817+00:00
p001        run_demo  step6_7                  tmp/GA_with_id.csv                                input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/GA_with_id.csv                                2352        1            csv            ok                           2026-05-28T10:07:26.598171+00:00
p001        run_demo  step6_7                  tmp/IRL_filtered.csv                              input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/IRL_filtered.csv                              7691        1            csv            ok                           2026-05-28T10:07:26.598266+00:00
p001        run_demo  step6_7                  tmp/repaired_ligand_data.csv                      input   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp/repaired_ligand_data.csv                      5011954     1            csv            ok                           2026-05-28T10:07:26.598075+00:00

============================================================
13. Step execution samples
============================================================
project_id  run_id    step_id                  started_at                        status   command_json                                                                                                                                                                                                                                                                                                                                                                                                                                            cwd                                            finished_at                       duration_seconds  exit_code  error  log_path                                                                                                       manifest_path                                                                                                        in_process_result_json                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  ingest_status  ingest_error  ingested_at                     
----------  --------  -----------------------  --------------------------------  -------  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  ---------------------------------------------  --------------------------------  ----------------  ---------  -----  -------------------------------------------------------------------------------------------------------------  -------------------------------------------------------------------------------------------------------------------  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  -------------  ------------  --------------------------------
p001        run_demo  step1                    2026-05-21T06:15:33.208984+00:00  success  ["/home/chenghua/miniconda3/envs/scidata/bin/python", "/home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB/src/step1.py", "--pubchem-dir", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/pubchem", "--metal-list", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/metal_list.txt", "--tmp-dir", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp"]  /home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB  2026-05-21T06:29:56.899742+00:00  863.691           0                 /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/logs/step1.log                    /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/manifests/step1.json                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ok                           2026-05-28T10:07:26.594845+00:00
p001        run_demo  step2                    2026-05-21T06:29:56.900610+00:00  success  ["/home/chenghua/miniconda3/envs/scidata/bin/python", "/home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB/src/step2.py", "--input-dir", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp", "--output-dir", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp"]                                                                                                                              /home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB  2026-05-21T06:30:05.706543+00:00  8.806             0                 /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/logs/step2.log                    /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/manifests/step2.json                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ok                           2026-05-28T10:07:26.597260+00:00
p001        run_demo  step4_5                  2026-05-21T06:30:05.707937+00:00  success  ["/home/chenghua/miniconda3/envs/scidata/bin/python", "/home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB/src/step4_5.py", "--input-dir", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp", "--output-dir", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp"]                                                                                                                            /home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB  2026-05-21T06:30:32.874741+00:00  27.167            0                 /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/logs/step4_5.log                  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/manifests/step4_5.json                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          ok                           2026-05-28T10:07:26.597522+00:00
p001        run_demo  step6_7                  2026-05-21T06:30:32.875693+00:00  success  ["/home/chenghua/miniconda3/envs/scidata/bin/python", "/home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB/src/step6_7.py", "--input-dir", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp", "--output-dir", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp"]                                                                                                                            /home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB  2026-05-21T06:30:51.286684+00:00  18.411            0                 /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/logs/step6_7.log                  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/manifests/step6_7.json                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          ok                           2026-05-28T10:07:26.597977+00:00
p001        run_demo  step8                    2026-05-21T06:30:51.288356+00:00  success  ["/home/chenghua/miniconda3/envs/scidata/bin/python", "/home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB/src/step8.py", "--input-dir", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp", "--output-dir", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp", "--metal-list", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/metal_list.txt"]          /home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB  2026-05-21T06:30:55.725685+00:00  4.437             0                 /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/logs/step8.log                    /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/manifests/step8.json                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ok                           2026-05-28T10:07:26.598432+00:00
p001        run_demo  step12                   2026-05-21T06:30:55.727899+00:00  success  ["/home/chenghua/miniconda3/envs/scidata/bin/python", "/home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB/src/step12.py", "--input-dir", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp", "--output-dir", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp", "--K", "30"]                                                                                                                /home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB  2026-05-21T06:31:07.339510+00:00  11.612            0                 /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/logs/step12.log                   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/manifests/step12.json                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           ok                           2026-05-28T10:07:26.595313+00:00
p001        run_demo  step13                   2026-05-21T06:31:07.342185+00:00  success  ["/home/chenghua/miniconda3/envs/scidata/bin/python", "/home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB/src/step_13_complete.py", "--input-dir", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp", "--output-dir", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp", "--kl-nl-only"]                                                                                                   /home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB  2026-05-21T06:31:43.075782+00:00  35.734            0                 /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/logs/step13.log                   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/manifests/step13.json                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           ok                           2026-05-28T10:07:26.596573+00:00
p001        run_demo  build_l3_embedding_ecfp  2026-05-21T06:42:53.337647+00:00  success  ["/home/chenghua/miniconda3/envs/scidata/bin/python", "/home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB/src/tools/build_L3_embedding_index.py", "--tmp-dir", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/tmp", "--out-dir", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/L3_embedding", "--backends", "ecfp"]                                                                      /home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB  2026-05-21T06:43:02.332769+00:00  8.995             0                 /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/logs/build_l3_embedding_ecfp.log  /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/manifests/build_l3_embedding_ecfp.json                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          ok                           2026-05-28T10:07:26.593202+00:00
p001        run_demo  zscore_metal_embedding   2026-05-21T06:43:02.333887+00:00  success  ["/home/chenghua/miniconda3/envs/scidata/bin/python", "/home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB/data/metal_embedding/zscore_element_features.py", "--input", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/metal_embedding/element_features.csv", "--output", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/metal_embedding/element_features_zscore.csv"]                /home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB  2026-05-21T06:43:02.775181+00:00  0.441             0                 /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/logs/zscore_metal_embedding.log   /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/manifests/zscore_metal_embedding.json                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           ok                           2026-05-28T10:07:26.602662+00:00
p001        run_demo  prepare_training_index   2026-05-21T06:43:02.776263+00:00  success  ["in_process:prepare_training_index"]                                                                                                                                                                                                                                                                                                                                                                                                                   /home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB  2026-05-21T06:43:02.778479+00:00  0.002             0                                                                                                                                /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/manifests/prepare_training_index.json   {"ok": true, "index_path": "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/index.json", "ligand_embedding": {"backend": "ecfp", "mode": null, "path": "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/L3_embedding/L3_embedding_ecfp.npz", "dir": "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/L3_embedding", "format": "npz", "variants": {"ecfp": {"path": "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/L3_embedding/L3_embedding_ecfp.npz", "file": "L3_embedding_ecfp.npz"}}, "npz_keys": ["dids", "smiles", "embeddings"]}}  ok                           2026-05-28T10:07:26.594147+00:00
p001        run_demo  prepare_training_config  2026-05-21T06:43:02.779567+00:00  success  ["in_process:prepare_training_config"]                                                                                                                                                                                                                                                                                                                                                                                                                  /home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB  2026-05-21T06:43:02.780160+00:00  0.001             0                                                                                                                                /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/manifests/prepare_training_config.json  {"ok": true, "config_path": "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/config.yaml", "integration": true}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   ok                           2026-05-28T10:07:26.593856+00:00
p001        run_demo  training_data            2026-05-21T06:43:02.780639+00:00  success  ["/home/chenghua/miniconda3/envs/scidata/bin/python", "-m", "training.data", "--training-dir", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training", "--index", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/index.json", "--seed", "42"]                                                                                                                                  /home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB  2026-05-21T06:43:04.109885+00:00  1.329             0                 /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/logs/training_data.log            /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/manifests/training_data.json                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    ok                           2026-05-28T10:07:26.600096+00:00
p001        run_demo  training_train           2026-05-21T06:43:04.112720+00:00  success  ["/home/chenghua/miniconda3/envs/scidata/bin/python", "-m", "training.train", "--config", "/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/config.yaml"]                                                                                                                                                                                                                                                          /home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB  2026-05-21T06:43:08.951662+00:00  4.839             0                 /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/logs/training_train.log           /home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/manifests/training_train.json                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   ok                           2026-05-28T10:07:26.601601+00:00

============================================================
14. Idempotence check: counts before and after rebuild-index
============================================================
---- before ----
table_name  n 
----------  --
models      2 
results     4 
artifacts   80
steps       13
{
  "ok": true,
  "summary": {
    "projects": 1,
    "runs": 1,
    "ga_sets": 1,
    "tasks": 2,
    "batches": 3,
    "errors": []
  }
}
---- after ----
table_name  n 
----------  --
models      2 
results     4 
artifacts   80
steps       13

============================================================
15. Repository query smoke test
============================================================

---- list_projects ----
[{'project_id': 'p001', 'project_root': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001', 'name': 'p001', 'created_at': None, 'updated_at': None, 'notes': None, 'ingest_status': 'ok', 'ingest_error': None, 'ingested_at': '2026-05-28T10:07:26.748757+00:00'}]

---- list_runs ----
[{'project_id': 'p001', 'run_id': 'run_demo', 'run_root': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo', 'created_at': '2026-05-21T06:13:47.371952+00:00', 'chemdb_repo_root': '/home/chenghua/ZhuoLi/ChemDBWebVersion/ChemDB', 'status': 'success', 'pubchem_file_count': 79, 'pipeline_core_ok': 1, 'pipeline_training_ok': 1, 'last_finished_at': '2026-05-21T06:43:08.953414+00:00', 'ga_binding_status': 'legacy_auto_ga', 'ingest_status': 'ok', 'ingest_error': None, 'ingested_at': '2026-05-28T10:07:26.749353+00:00'}]

---- list_ga_sets ----
[{'ga_set_id': 'fixture_run_demo', 'ga_set_root': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/ga_sets/fixture_run_demo', 'name': 'run_demo GA (small)', 'description': 'pytest fixture', 'created_at': '2026-05-21T00:00:00+00:00', 'updated_at': '2026-05-21T00:00:00+00:00', 'tags_json': '["fixture"]', 'source_json': '{"type": "fixture"}', 'ingest_status': 'ok', 'ingest_error': None, 'ingested_at': '2026-05-28T10:07:26.766175+00:00'}]

---- list_models ----
[{'model_id': 'p001_run_demo_best', 'project_id': 'p001', 'run_id': 'run_demo', 'checkpoint_path': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/ckpts/best.pt', 'checkpoint_stem': 'best', 'config_path': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/config.yaml', 'index_path': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/index.json', 'embedding_backend': 'ecfp', 'l3_embedding_path': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/L3_embedding/L3_embedding_ecfp.npz', 'metal_embedding_path': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/metal_embedding/element_features_zscore.csv', 'history_path': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/ckpts/history.json', 'best_epoch': 3, 'val_mrr': 0.07604873093465964, 'val_recall_at_1': 0.0, 'val_recall_at_5': 0.13333333333333333, 'ga_set_id': None, 'ga_version_id': None, 'ga_binding_checksum': None, 'size_bytes': 236892252, 'ingest_status': 'ok', 'ingest_error': None, 'ingested_at': '2026-05-28T10:07:26.765574+00:00'}, {'model_id': 'p001_run_demo_last', 'project_id': 'p001', 'run_id': 'run_demo', 'checkpoint_path': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/ckpts/last.pt', 'checkpoint_stem': 'last', 'config_path': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/config.yaml', 'index_path': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/index.json', 'embedding_backend': 'ecfp', 'l3_embedding_path': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/L3_embedding/L3_embedding_ecfp.npz', 'metal_embedding_path': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/data/metal_embedding/element_features_zscore.csv', 'history_path': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo/training/ckpts/history.json', 'best_epoch': 3, 'val_mrr': 0.07604873093465964, 'val_recall_at_1': 0.0, 'val_recall_at_5': 0.13333333333333333, 'ga_set_id': None, 'ga_version_id': None, 'ga_binding_checksum': None, 'size_bytes': 236892252, 'ingest_status': 'ok', 'ingest_error': None, 'ingested_at': '2026-05-28T10:07:26.765799+00:00'}]

---- list_model_after_tasks ----
[{'task_id': 'ni_lkb_p_ni', 'task_dir': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_tasks/ni_lkb_p_ni', 'metal': 'Ni', 'embedding': 'ecfp', 'candidate_file': 'candidates.csv', 'positive_file': 'positives.csv', 'negative_file': 'negatives.csv', 'notes': 'Converted from ChemDB_restructured legacy quest ni_lkb_p_ni.csv (344 candidates, 11 positives, no explicit negatives).', 'candidate_count': 344, 'positive_count': 11, 'negative_count': 0, 'ingest_status': 'ok', 'ingest_error': None, 'ingested_at': '2026-05-28T10:07:26.769942+00:00'}, {'task_id': 'pd_lkb_p_cluster', 'task_dir': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_tasks/pd_lkb_p_cluster', 'metal': 'Pd', 'embedding': 'ecfp', 'candidate_file': 'candidates.csv', 'positive_file': 'positives.csv', 'negative_file': 'negatives.csv', 'notes': 'Converted from ChemDB_restructured legacy quest pd_lkb_p_cluster.csv', 'candidate_count': 343, 'positive_count': 10, 'negative_count': 1, 'ingest_status': 'ok', 'ingest_error': None, 'ingested_at': '2026-05-28T10:07:26.771137+00:00'}]

---- list_batches(ni_lkb_p_ni) ----
[{'task_id': 'ni_lkb_p_ni', 'batch_id': 'batch_dummy_ckpts', 'output_dir': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/ni_lkb_p_ni/batch_dummy_ckpts', 'command': 'evaluate-models', 'model_runs_file': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/ni_lkb_p_ni/batch_dummy_ckpts/model_runs.csv', 'finished_at': '2026-05-21T09:17:57.873556+00:00', 'n_models': 2, 'n_success': 2, 'best_model_id': 'run_demo_last', 'summary_path': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/ni_lkb_p_ni/batch_dummy_ckpts/model_selection_summary.csv', 'ingest_status': 'ok', 'ingest_error': None, 'ingested_at': '2026-05-28T10:07:26.768704+00:00'}, {'task_id': 'ni_lkb_p_ni', 'batch_id': 'p001_run_demo', 'output_dir': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/ni_lkb_p_ni/p001_run_demo', 'command': 'evaluate-model', 'model_runs_file': None, 'finished_at': '2026-05-21T09:05:05.892604+00:00', 'n_models': None, 'n_success': None, 'best_model_id': None, 'summary_path': None, 'ingest_status': 'ok', 'ingest_error': None, 'ingested_at': '2026-05-28T10:07:26.770278+00:00'}]

---- list_batches(pd_lkb_p_cluster) ----
[{'task_id': 'pd_lkb_p_cluster', 'batch_id': 'p001_run_demo', 'output_dir': '/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/model_after_results/pd_lkb_p_cluster/p001_run_demo', 'command': 'evaluate-model', 'model_runs_file': None, 'finished_at': '2026-05-21T09:04:32.709394+00:00', 'n_models': None, 'n_success': None, 'best_model_id': None, 'summary_path': None, 'ingest_status': 'ok', 'ingest_error': None, 'ingested_at': '2026-05-28T10:07:26.771515+00:00'}]

============================================================
16. Check real workspace DB was not created/modified by this validation
============================================================

============================================================
VALIDATION DONE
Please paste the complete output back.
============================================================
```
