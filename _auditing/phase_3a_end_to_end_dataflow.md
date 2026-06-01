# Phase 3A 端到端数据流与脚本依赖说明

| 字段 | 值 |
|------|-----|
| 类型 | 架构/数据流说明 |
| 日期 | 2026-05-29 |
| 目标 | 说明从原始数据到 Streamlit 展示的完整链路、脚本职责、文件内容与依赖关系 |
| 范围 | `app/services`、`ChemDB/src`、`app/storage`、`app/streamlit_ui`、`workspace` |
| 边界 | Phase 3A 只读展示，不改变科研逻辑 |

---

## 1. 一张图看全链路

```text
[原始输入数据]
  workspace/projects/<project>/runs/<run>/data/pubchem/*.csv
  + metal_list.txt + p_elements_list.txt + metal_embedding/element_features.csv
        |
        v
[Run 编排层] app/services/orchestrator.py + step_registry.py
  step1 -> step2 -> require_bound_ga -> apply_ga_to_run -> step6_7 -> step8 -> step12 -> step13
        |
        v
[训练链路]
  training_manager.prepare_training_index/config
  -> python -m training.data
  -> python -m training.train
        |
        v
[Model-After 评估结果]
  app/services/model_after_manager.py
  -> workspace/model_after_tasks/*
  -> workspace/model_after_results/*/*
        |
        v
[SQLite 只读登记]
  python -m app.storage.cli rebuild-index
  (app/storage/ingest.py -> schema.sql)
        |
        v
[Phase 3A 网页]
  app/streamlit_ui/Home.py + pages/*
  -> app/storage/repositories.py + previews.py
  -> workspace/chemdb.sqlite + 文件预览
```

---

## 2. 阶段 0：目录与资源分层（谁是“源数据”）

当前系统有三类关键目录：

1. **Run 工作区（每次运行隔离）**
   - `workspace/projects/<project_id>/runs/<run_id>/`
   - 核心子目录：
     - `data/`（输入与 embedding 文件）
     - `tmp/`（流程中间产物）
     - `training/`（训练索引、配置、ckpt）
     - `manifests/`（步骤执行与状态）
     - `logs/`

2. **Workspace 级资源（跨 run 复用）**
   - `workspace/ga_sets/<ga_set_id>/versions/<ga_version_id>/GA_with_id.csv`
   - `workspace/model_after_tasks/<task_id>/`
   - `workspace/model_after_results/<task_id>/<batch_id>/`

3. **Phase 2A/3A 元数据层**
   - `workspace/chemdb.sqlite`
   - 用于前端只读查询，不替代原始文件

---

## 3. 阶段 1：从原始数据到 run 产物（orchestrator + ChemDB/src）

### 3.1 编排入口与运行上下文

- 入口：`app/services/orchestrator.py`
- 运行上下文：`app/services/run_context.py`
  - 将 `run_root` 解析为 `data_dir/tmp_dir/training_dir/manifest_dir` 等
- 步骤定义：`app/services/step_registry.py`
  - 明确每一步的：
    - `required_inputs`
    - `expected_outputs`
    - `optional_outputs`
    - 执行命令模板

### 3.2 core pipeline 的步骤与依赖（Phase 1E 后）

`PIPELINE_ORDER`：

1. `step1`（`ChemDB/src/step1.py`）
   - 读：`data/pubchem/*.csv`, `data/metal_list.txt`
   - 写：`tmp/complex_data.csv`, `tmp/ligand_data.csv`

2. `step2`（`ChemDB/src/step2.py`）
   - 读：`tmp/ligand_data.csv`
   - 写：`tmp/repaired_ligand_data.csv`

3. `require_bound_ga`（in-process，`ga_registry_manager.require_bound_ga`）
   - 校验必须存在：
     - `manifests/ga_binding.json`
     - `tmp/GA_with_id.csv`

4. `apply_ga_to_run`（`ChemDB/src/step4_5.py --mode apply-ga`）
   - 读：`tmp/repaired_ligand_data.csv`, `tmp/GA_with_id.csv`
   - 写：`tmp/ligand_with_gac.csv`, `tmp/IRL_filtered.csv`

5. `step6_7`（`ChemDB/src/step6_7.py`）
   - 读：`repaired_ligand_data.csv`, `GA_with_id.csv`, `IRL_filtered.csv`
   - 写：`tmp/fragments.csv`

6. `step8`（`ChemDB/src/step8.py`）
   - 读：前序 csv + `metal_list.txt`
   - 写：`tmp/neo4j_*.csv`

7. `step12`（`ChemDB/src/step12.py`）
   - 读：neo4j 关系 csv + GAC 数据
   - 写：`tmp/l3_l5.json`, `tmp/l5_l3.json`, `tmp/l3_gac.json`, `tmp/m_l3_pairs.csv` 等

8. `step13`（`ChemDB/src/step_13_complete.py`）
   - 读：step12 的 json/csv
   - 写：`tmp/step13_kl_nl_samples.csv`（训练关键输入）

### 3.3 每一步的状态如何落盘

`orchestrator.run_step()` 会为每一步写 `manifests/<step_id>.json`，包含：

- `status` / `started_at` / `finished_at` / `duration_seconds`
- `command` / `cwd` / `log_file`
- `inputs` / `outputs`（路径、存在性、大小）

这部分是后续 SQLite `step_executions` 与 `artifacts` 的主要来源。

---

## 4. 阶段 2：GA 资源如何生成、绑定、继承到 run

对应脚本：`app/services/ga_registry_manager.py`

### 4.1 Workspace 级 GA 版本（源头）

- 存放位置：
  - `{workspace_root}/ga_sets/{ga_set_id}/versions/{ga_version_id}/GA_with_id.csv`

### 4.2 绑定到 run（关键动作）

命令：

```bash
python -m app.services.ga_registry_manager bind-ga-version-to-run \
  --workspace-root <workspace_root> \
  --run-root <run_root> \
  --ga-set-id <ga_set_id> \
  --ga-version-id <ga_version_id>
```

效果：

1. 从 workspace 级版本复制 GA 到 run 快照：
   - `run/tmp/GA_with_id.csv`
2. 写绑定清单：
   - `run/manifests/ga_binding.json`
   - 包含 `ga_set_id/ga_version_id/checksum/source_ga_csv/materialized_path`
3. 写 stale 报告（若下游产物可能受影响）：
   - `run/manifests/ga_stale_report.json`

### 4.3 checksum 的角色

- `checksum` 基于 `GA_with_id.csv` 计算（`sha256:...`）
- 在 SQLite 中落入：
  - `run_ga_bindings.checksum`
  - `models.ga_binding_checksum`
- 用于追踪“模型到底基于哪份 GA 内容训练”

---

## 5. 阶段 3：训练链路（从 step13 到模型）

涉及脚本：

- `app/services/training_manager.py`（in-process 准备）
- `python -m training.data`（生成 `*_records.pkl`）
- `python -m training.train`（输出 ckpt）

### 5.1 准备索引与配置

`training_manager.prepare_training_index(ctx)` 生成：

- `training/index.json`
  - 指向：
    - `tmp/step13_kl_nl_samples.csv`
    - `tmp/m_l3_pairs.csv`
    - `tmp/l3_gac.json`
    - `data/L3_embedding/L3_embedding_ecfp.npz`
    - `data/metal_embedding/element_features_zscore.csv`

`training_manager.prepare_training_config(ctx)` 生成：

- `training/config.yaml`

### 5.2 训练数据与模型

- `training.data` 读 `index.json`，写：
  - `training/split_index.json`
  - `training/train_records.pkl`
  - `training/val_records.pkl`
  - `training/test_records.pkl`

- `training.train` 读 config + pkl，写：
  - `training/ckpts/best.pt`
  - `training/ckpts/last.pt`
  - `training/ckpts/history.json`

---

## 6. 阶段 4：Model-After 任务与结果

脚本：`app/services/model_after_manager.py`

### 6.1 Task 资源

- `workspace/model_after_tasks/<task_id>/task.json`
- 同目录下 `candidates.csv` / `positives.csv` / `negatives.csv`

### 6.2 结果资源

- `workspace/model_after_results/<task_id>/<batch_id>/`
- 常见文件：
  - `evaluate_model(s)_manifest.json`
  - `model_selection_summary.csv`
  - `best_model.json`
  - `models/<model_id>/ranking.csv`
  - `models/<model_id>/metrics.json`
  - `models/<model_id>/report.json`

---

## 7. 阶段 5：SQLite ingest（文件 -> 表）

入口：`python -m app.storage.cli rebuild-index`
核心：`app/storage/ingest.py`

### 7.1 扫描顺序

1. `workspace/projects/*/runs/*` -> `register_run()`
   - 子步骤：
     - `ingest_run`
     - `ingest_run_ga_binding`
     - `ingest_step_executions`
     - `ingest_key_artifacts`
     - `ingest_model`
2. `workspace/ga_sets/*` -> `ingest_ga_set`
3. `workspace/model_after_tasks/*` -> `ingest_model_after_task`
4. `workspace/model_after_results/*/*` -> `ingest_model_after_batch`

### 7.2 关键表与来源

| 表 | 主要来源文件 |
|----|--------------|
| `runs` | `manifests/run.json`, `pipeline_*.json`, 路径解析 |
| `run_ga_bindings` | `manifests/ga_binding.json`（或 legacy/missing 推断） |
| `step_executions` | `manifests/<step>.json` |
| `artifacts` | step manifest 的 `inputs/outputs` + key artifacts |
| `models` | `training/ckpts/*.pt`, `training/index.json`, `config.yaml`, `history.json`, `run_ga_bindings` |
| `model_after_tasks` | `model_after_tasks/<task>/task.json` + csv 行数 |
| `model_after_batches` | `model_after_results/<task>/<batch>/*.json/*.csv` |
| `model_after_model_results` | `metrics.json` + `model_selection_summary.csv` |

### 7.3 `registry_model_id` 如何建立

`ingest._resolve_registry_model_id()` 尝试按：

- `project_id + run_id + checkpoint_path`

映射到 `models.model_id`，并写入：

- `model_after_model_results.registry_model_id`

这是前端结果页正确 join 模型元数据的关键字段。

---

## 8. 阶段 6：Streamlit 页面如何消费 SQLite 与文件

入口：`app/streamlit_ui/Home.py`
页面状态：`app/streamlit_ui/state.py`

### 8.1 DB 连接策略

- `state.get_db_connection(read_only=True)` 使用 `app.storage.db.connect_db`
- DB 默认路径：
  - `CHEMDB_DB_PATH` 或 `workspace/chemdb.sqlite`
- DB 不存在时返回 `None`，页面提示而不崩溃

### 8.2 页面到查询函数映射

| 页面 | 查询函数 |
|------|----------|
| Dashboard | `get_registry_summary`, `list_runs`, `list_models`, `list_all_batches` |
| GA Library | `list_ga_sets`, `list_ga_versions`, `list_runs_by_ga_version` |
| Task Library | `list_model_after_tasks`, `get_model_after_task`, `list_batches` |
| Model Library | `list_runs`, `list_models`, `get_model` |
| Model Selector Results | `list_model_after_tasks`, `list_batches`, `get_batch`, `list_model_results_with_models` |
| Settings | `get_registry_summary` + rebuild 命令 |

### 8.3 文件预览链路

`app/storage/previews.py`：

- `preview_ga_version` -> 读 `ga_versions.ga_csv_path`
- `preview_task_files` -> `task_dir + candidate/positive/negative_file`
- `preview_model_result` -> `ranking_path`, `metrics_path`, `report_path`
- CSV 只读前 100 行（`pandas.read_csv(nrows=...)`）

---

## 9. 关键依赖关系（按“谁依赖谁”）

1. **网页显示依赖 SQLite 与源文件**
   - 页面主表来自 SQLite
   - 预览来自 SQLite 中的路径字段再读文件

2. **SQLite 依赖 manifests 和 run 产物**
   - 如果某步骤未跑，相关表字段会缺失或为 `missing/stub`

3. **模型 GA 元数据依赖 run 绑定**
   - `models.ga_set_id/ga_version_id/ga_binding_checksum` 来源于 `run_ga_bindings`

4. **结果页模型 join 依赖 `registry_model_id`**
   - 必须 `model_after_model_results.registry_model_id = models.model_id`

5. **Rebuild Index 是页面数据刷新唯一写入口（Phase 3A）**
   - Settings 按钮仅调用 `python -m app.storage.cli rebuild-index`

---

## 10. 运维与排障（当前版本）

### 10.1 推荐顺序

```bash
# 1) 先更新 sqlite
python -m app.storage.cli rebuild-index --workspace-root workspace --db-path workspace/chemdb.sqlite

# 2) 启动前端
PYTHONPATH=. streamlit run app/streamlit_ui/Home.py
```

### 10.2 常见问题

1. `DB not found`
   - 先执行 `rebuild-index`

2. 结果页显示 `missing_model`
   - 说明该 result 的 `registry_model_id` 未匹配到 `models`
   - 常见原因：checkpoint 路径不一致或模型未被 ingest

3. 预览文件不存在
   - SQLite 路径记录存在，但源文件被移动/清理
   - 重新 `rebuild-index` 并检查对应 run/task 目录

---

## 11. 与 Phase 3A 边界的一致性

本链路中，Phase 3A 页面严格为只读：

- 不直接调用 orchestrator / training / model-after 执行命令
- 不改 `ChemDB/src/**`
- 唯一允许写操作为 Settings 的 `rebuild-index`（更新 sqlite）

该设计保证前端只做“资源组织与展示”，不会改变科研流程语义。
