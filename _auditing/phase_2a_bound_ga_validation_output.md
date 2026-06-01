# Phase 2A Bound-GA 验证输出

| 字段 | 值 |
|------|-----|
| 日期 | 2026-05-28 |
| 脚本 | `scripts/phase2a_bound_ga_validation.sh` |
| DB | `/tmp/chemdb_phase2a_bound_validation.sqlite` |
| Workspace（复制） | `/tmp/chemdb_workspace_bound_validation` |
| Run | `p001/run_demo` |
| GA 绑定 | `fixture_run_demo` / `v001` |

## 结论摘要

| 检查项 | 结果 |
|--------|------|
| bind-ga-version-to-run | 成功，`status=bound`，126 GA |
| `run_ga_bindings` | 1 行，`bound`，checksum 正确 |
| `runs.ga_binding_status` | `bound`（非 legacy `legacy_auto_ga`） |
| `models` GA 元数据 | 2 个 checkpoint 均继承 `fixture_run_demo/v001` + checksum |
| model-after join | 4 行结果全部 `matched_model`，GA 字段可 join |
| rebuild-index 幂等 | 行数不变（bindings 1, models 2, results 4, artifacts 80, steps 13） |
| 真实 workspace DB | 未创建 `workspace/chemdb.sqlite` |

### 说明

- bind 使用 `--workspace-root "$TMP_WORKSPACE"`；`source_ga_csv` 指向复制 workspace 内 GA（非 repo `workspace/ga_sets/`）。
- bind 后生成 `ga_stale_report.json`，列出 39 个 downstream 产物为 stale（预期行为，未在本验证中 re-run pipeline）。
- 复制 workspace 上的 `rebuild-index` 正确 ingest 了 `manifests/ga_binding.json`。

## 完整终端输出

```text
============================================================
0. Prepare temporary copied workspace
============================================================
TMP_WORKSPACE=/tmp/chemdb_workspace_bound_validation
RUN_ROOT=/tmp/chemdb_workspace_bound_validation/projects/p001/runs/run_demo
GA_SET_ID=fixture_run_demo
GA_VERSION_ID=v001

============================================================
1. Bind workspace GA version to copied run
============================================================
{
  "run_root": "/tmp/chemdb_workspace_bound_validation/projects/p001/runs/run_demo",
  "ga_set_id": "fixture_run_demo",
  "ga_version_id": "v001",
  "source_ga_csv": "/tmp/chemdb_workspace_bound_validation/ga_sets/fixture_run_demo/versions/v001/GA_with_id.csv",
  "materialized_path": "tmp/GA_with_id.csv",
  "checksum": "sha256:edc18600d0819d62511afa6de4a410c27b0830f7a0d6b8f92f16df99921e56ce",
  "num_ga": 126,
  "bound_at": "2026-05-28T10:11:33.278961+00:00",
  "bound_by": "ga_registry_manager",
  "status": "bound",
  "stale_report": "/tmp/chemdb_workspace_bound_validation/projects/p001/runs/run_demo/manifests/ga_stale_report.json"
}

---- ga_binding.json ----
{
  "run_root": "/tmp/chemdb_workspace_bound_validation/projects/p001/runs/run_demo",
  "ga_set_id": "fixture_run_demo",
  "ga_version_id": "v001",
  "source_ga_csv": "/tmp/chemdb_workspace_bound_validation/ga_sets/fixture_run_demo/versions/v001/GA_with_id.csv",
  "materialized_path": "tmp/GA_with_id.csv",
  "checksum": "sha256:edc18600d0819d62511afa6de4a410c27b0830f7a0d6b8f92f16df99921e56ce",
  "num_ga": 126,
  "bound_at": "2026-05-28T10:11:33.278961+00:00",
  "bound_by": "ga_registry_manager",
  "status": "bound"
}
---- materialized GA_with_id.csv ----
-rw-rw-r-- 1 chenghua chenghua 2.5K May 28 10:11 /tmp/chemdb_workspace_bound_validation/projects/p001/runs/run_demo/tmp/GA_with_id.csv
GA_SMILES,GA_ID
c1cc[cH-]c1,G218265
c1ccncc1,G289719
c1ccccc1,G874135
c1cnoc1,G296672

---- stale report, if created ----
{
  "generated_at": "2026-05-28T10:11:33.280390+00:00",
  "trigger": "bind-ga-version-to-run",
  "previous_binding": null,
  "new_binding": {
    "ga_set_id": "fixture_run_demo",
    "ga_version_id": "v001",
    "checksum": "sha256:edc18600d0819d62511afa6de4a410c27b0830f7a0d6b8f92f16df99921e56ce"
  },
  "stale_files": [
    {
      "path": "tmp/ligand_with_gac.csv",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/IRL_filtered.csv",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/IRL_filtered_cleaned.csv",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/step4_5_stats.json",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/fragments.csv",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/step6_7_stats.json",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/l3_l5.json",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/l5_l3.json",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/l5_freq_weight.json",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/l3_gac.json",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/m_l3_pairs.csv",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/metal_l3_index.csv",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/step12_stats.json",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/step13_kl_nl_samples.csv",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/step13_kl_nl_samples.stats.csv",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/step13_stats.json",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "data/L3_embedding/L3_embedding_ecfp.npz",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "training/index.json",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "training/config.yaml",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "training/split_index.json",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "training/train_records.pkl",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "training/val_records.pkl",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "training/test_records.pkl",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "training/ckpts/best.pt",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "training/ckpts/last.pt",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "training/ckpts/history.json",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/neo4j_complexes.csv",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/neo4j_fragments.csv",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/neo4j_irl.csv",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/neo4j_l1_l2_relationships.csv",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/neo4j_l1_l3_relationships.csv",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/neo4j_l2_l3_relationships.csv",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/neo4j_l3_l4_relationships.csv",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/neo4j_l4_l5_relationships.csv",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/neo4j_ligands.csv",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/neo4j_m_l1_relationships.csv",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/neo4j_metals.csv",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "tmp/neo4j_repaired_ligands.csv",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "training/ckpts/best.pt",
      "exists": true,
      "reason": "downstream_of_ga"
    },
    {
      "path": "training/ckpts/last.pt",
      "exists": true,
      "reason": "downstream_of_ga"
    }
  ],
  "recommended_actions": [
    "python -m app.services.ga_registry_manager apply-ga-to-run --run-root <run>",
    "python -m app.services.orchestrator run-pipeline --run-root <run> --through step13"
  ]
}
============================================================
2. Rebuild SQLite from copied workspace
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

============================================================
3. Check run GA binding status
============================================================
project_id  run_id    status  ga_set_id         ga_version_id  checksum                                                                 num_ga  materialized_path
----------  --------  ------  ----------------  -------------  -----------------------------------------------------------------------  ------  ------------------
p001        run_demo  bound   fixture_run_demo  v001           sha256:edc18600d0819d62511afa6de4a410c27b0830f7a0d6b8f92f16df99921e56ce  126     tmp/GA_with_id.csv

============================================================
4. Check runs table GA status
============================================================
project_id  run_id    ga_binding_status  status   pipeline_core_ok  pipeline_training_ok
----------  --------  -----------------  -------  ----------------  --------------------
p001        run_demo  bound              success  1                 1

============================================================
5. Check models inherit GA metadata
============================================================
model_id            project_id  run_id    checkpoint_stem  ga_set_id         ga_version_id  ga_binding_checksum
------------------  ----------  --------  ---------------  ----------------  -------------  -----------------------------------------------------------------------
p001_run_demo_best  p001        run_demo  best             fixture_run_demo  v001           sha256:edc18600d0819d62511afa6de4a410c27b0830f7a0d6b8f92f16df99921e56ce
p001_run_demo_last  p001        run_demo  last             fixture_run_demo  v001           sha256:edc18600d0819d62511afa6de4a410c27b0830f7a0d6b8f92f16df99921e56ce

============================================================
6. Check model-after result registry join still works
============================================================
task_id           batch_id           raw_model_id   registry_model_id   model_join     ga_set_id         ga_version_id  status   mrr                 hit_at_5           selection_score     rank_among_models
----------------  -----------------  -------------  ------------------  -------------  ----------------  -------------  -------  ------------------  -----------------  ------------------  -----------------
ni_lkb_p_ni       batch_dummy_ckpts  run_demo_last  p001_run_demo_last  matched_model  fixture_run_demo  v001           success  0.171211730938318   0.187195689615044  0.171211730938318   1
ni_lkb_p_ni       batch_dummy_ckpts  run_demo_best  p001_run_demo_best  matched_model  fixture_run_demo  v001           success  0.0155332868367695  0.0                0.0155332868367695  2
ni_lkb_p_ni       p001_run_demo      run_demo       p001_run_demo_best  matched_model  fixture_run_demo  v001           success  0.0155332868367695  0.0
pd_lkb_p_cluster  p001_run_demo      run_demo       p001_run_demo_best  matched_model  fixture_run_demo  v001           success  0.0159977312363552  0.0

============================================================
7. Idempotence check
============================================================
---- before ----
table_name       n
---------------  --
run_ga_bindings  1
models           2
results          4
artifacts        80
steps            13
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
table_name       n
---------------  --
run_ga_bindings  1
models           2
results          4
artifacts        80
steps            13

============================================================
8. Confirm real workspace DB not created
============================================================

============================================================
BOUND-GA VALIDATION DONE
============================================================
```
