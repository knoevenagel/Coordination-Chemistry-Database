# Phase 2A 实施记录：SQLite Read-only Registry

**日期**：2026-05-28  
**状态**：已完成

## 目标

实现 Phase 2A **只读 ingest** SQLite metadata layer，登记 workspace 中 Project、Run、GA、Step、Artifact、Model、ModelAfter 元数据；**不改变** orchestrator / ga_registry / model_after pipeline 运行时行为。

## 新增文件

| 路径 | 说明 |
|------|------|
| `app/storage/schema.sql` | 11 表 DDL + 索引 |
| `app/storage/db.py` | `connect_db`, `init_db`, `execute_schema` |
| `app/storage/ingest.py` | filesystem → SQLite 纯函数 + `rebuild_index` |
| `app/storage/repositories.py` | 查询 helper（供 Phase 3 Streamlit） |
| `app/storage/cli.py` | `init-db`, `register-*`, `rebuild-index` |
| `app/storage/__init__.py`, `__main__.py` | 包入口 |
| `tests/storage/test_phase2a_sqlite_registry.py` | 12 项单元测试 |
| `.gitignore` | 增加 `workspace/chemdb.sqlite` |

## Schema 摘要

| 表 | 主键 | 主要外键 |
|----|------|----------|
| `projects` | `project_id` | — |
| `runs` | `(project_id, run_id)` | → `projects` |
| `ga_sets` | `ga_set_id` | — |
| `ga_versions` | `(ga_set_id, ga_version_id)` | → `ga_sets` |
| `run_ga_bindings` | `(project_id, run_id)` | → `runs`；`ga_set_id` 可空 |
| `step_executions` | `(project_id, run_id, step_id, started_at)` | → `runs` |
| `artifacts` | `(project_id, run_id, step_id, rel_path, role)` | → `runs` |
| `models` | `model_id` | → `runs`；`UNIQUE(project_id, run_id, checkpoint_path)` |
| `model_after_tasks` | `task_id` | — |
| `model_after_batches` | `(task_id, batch_id)` | → `model_after_tasks` |
| `model_after_model_results` | `(task_id, batch_id, model_id)` | → `batches` |

所有表使用 `INSERT ... ON CONFLICT DO UPDATE` 支持幂等 ingest。DB **只存路径、状态、摘要指标**，不存文件内容。

## Ingest 字段来源

| 实体 | 主要来源 |
|------|----------|
| `projects` | `workspace/projects/{id}/`；可选 `project.json` |
| `runs` | `manifests/run.json`, `pipeline_core.json`, `pipeline_training.json`；`project_id`/`run_id` 从路径解析 |
| `ga_sets` / `ga_versions` | `ga_sets/{id}/ga_set.json`, `versions/{vid}/ga_version.json` |
| `run_ga_bindings` | `manifests/ga_binding.json`；缺失时见 legacy 策略 |
| `step_executions` | `manifests/{step_id}.json`（跳过 run/pipeline/ga manifest） |
| `artifacts` | step manifest `inputs`/`outputs` + `ingest_key_artifacts` 关键路径 |
| `models` | `training/index.json`, `config.yaml`, `ckpts/*.pt`, `history.json` + join binding |
| `model_after_tasks` | `task.json` + CSV 行数 |
| `model_after_batches` | batch 目录 manifest + `model_selection_manifest.json` |
| `model_after_model_results` | `metrics.json`, `model_selection_summary.csv`；多模型在 `models/{id}/` |

## Legacy run 策略

| 条件 | `ga_binding_status` / binding `status` |
|------|----------------------------------------|
| 有 `ga_binding.json` | `bound` |
| 无 binding，有 `step4_5.json` | `legacy_auto_ga`（如 `run_demo`） |
| 无 binding，无 step4_5 | `missing` |

不抛异常；`ga_set_id` / `ga_version_id` 可为 NULL。

## `model_id` 合成规则

`{project_id}_{run_id}_{checkpoint_stem}`，例如 `p001_run_demo_best`。

## CLI 示例

```bash
python -m app.storage.cli init-db --db-path workspace/chemdb.sqlite

python -m app.storage.cli register-run \
  --run-root workspace/projects/p001/runs/run_demo \
  --db-path workspace/chemdb.sqlite

python -m app.storage.cli rebuild-index \
  --workspace-root workspace \
  --db-path workspace/chemdb.sqlite
```

## 测试结果

```bash
./scripts/run_pytest.sh tests/storage/test_phase2a_sqlite_registry.py -v   # 12 passed
./scripts/run_pytest.sh -q                                                  # 61 passed, 12 deselected
```

真实 workspace ingest（**输出到临时 DB，未写 workspace/chemdb.sqlite**）：

```bash
python -m app.storage.cli rebuild-index \
  --workspace-root workspace \
  --db-path /tmp/chemdb_phase2a_test.sqlite
```

结果：`projects=1`, `runs=1`, `ga_sets=1`, `tasks=2`, `batches=3`, `errors=[]`。

`run_demo`：`ga_binding_status=legacy_auto_ga`；models 含 `p001_run_demo_best` / `p001_run_demo_last`（若存在 ckpt）。

## 已知限制

1. **绝对路径**：ingest 原样存储 manifest/index 中的绝对路径；搬迁 workspace 需 re-ingest。
2. **无 runtime hook**：pipeline 成功/失败不会自动写 DB。
3. **stub 行**：binding 引用缺失的 ga_set、或 batch 无 task 目录时，可写入 `ingest_status=stub` 的占位行以满足 FK。
4. **`registry_model_id`**：model-after result 尽量 join `models`；join 失败仍保留 `model_run_root` / `checkpoint_path`。
5. **未 ingest 历史 ga_stale_report** 为多行 binding 事件（Phase 2B 可选）。

## Phase 2B 建议

1. 可选 `--write-db` hook（orchestrator 每步完成后 upsert）。
2. `run.json` 写入 `project_id` / `run_id`。
3. `training/model.json` 或 post-train manifest。
4. model-after `metrics.json` 增加 `ga_set_id` / `ga_version_id`。
5. Streamlit 只读 dashboards 基于 `repositories.py`。

## 关联

- [phase_2_readiness_after_ga_workspace.md](./phase_2_readiness_after_ga_workspace.md)
