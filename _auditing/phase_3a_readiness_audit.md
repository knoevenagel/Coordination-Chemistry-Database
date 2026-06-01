# Phase 3A Readiness Audit — Resource-Centric Read-Only Streamlit Dashboard

| 字段 | 值 |
|------|-----|
| 类型 | 审计 / Phase 3A 就绪评估（**不修改代码**） |
| 日期 | 2026-05-28 |
| 范围 | Phase 3A 只读 Streamlit dashboard（不写 pipeline、不训练、不编辑 GA、不跑 evaluation） |
| 依据 | `app/storage/`、`app/services/`、`workspace/` 真实数据、Phase 2A 验证输出、pytest 62 passed |
| 关联 | [frontend.md](./frontend.md)（Phase 3 总体设计）、[phase_2a_implementation.md](./phase_2a_implementation.md)、[phase_2a_validation_output.md](./phase_2a_validation_output.md) |
| 结论 | **可以进入 Phase 3A**；repository 层与 schema 已覆盖约 80% 只读需求，但缺 join helper、文件 preview helper、Streamlit 骨架与 streamlit 依赖 |

---

## 执行摘要

Phase 2A SQLite registry 已通过 legacy run、bound-GA run、model-after `registry_model_id` join、rebuild-index 幂等等验收。Phase 3A 目标是在此之上提供 **resource-centric 只读 Streamlit dashboard**，展示 GA Library、Task Library、Model Library、Model Selector Results、Dashboard 和 Settings。

**当前缺口（P0）**：

1. 无 `app/streamlit_ui/` 目录
2. scidata 环境未安装 streamlit
3. `workspace/chemdb.sqlite` 本地不存在（需 rebuild-index）
4. `repositories.py` 缺 model join 与 dashboard 聚合查询
5. 无 CSV/JSON preview helper

**明确不做**：启动 pipeline、训练、编辑 GA、创建 task、运行 evaluation、修改 `ChemDB/src`。

---

## A. 当前状态判断

### A.1 已就绪

| 能力 | 状态 | 说明 |
|------|------|------|
| SQLite registry ingest | ✅ | `rebuild-index` 可从 `workspace/` 全量 UPSERT，幂等已验证 |
| Schema 11 表 | ✅ | 覆盖 GA / Task / Model / Result / Run / Artifact / Step |
| 基础 repository 查询 | ✅ | 12 个 list/get 函数，足够支撑第一版列表页 |
| `registry_model_id` ingest | ✅ | ingest 时写入；真实 DB 4 条 result 均可 join `models` |
| 真实 workspace 数据 | ✅ | 1 project、1 run、1 ga_set/v001、2 tasks、3 batches、4 results、2 models |
| pandas | ✅ | scidata 环境 `pandas 2.0.3` |
| DB 连接层 | ✅ | `connect_db(read_only=True)` 支持只读 URI 模式 |

### A.2 未就绪 / 缺口

| 缺口 | 影响 | 优先级 |
|------|------|--------|
| **无 `app/streamlit_ui/`** | 前端从零开始 | P0 |
| **streamlit 未安装** | 无法启动 UI | P0 |
| **`workspace/chemdb.sqlite` 不存在** | 首次打开 UI 需提示 rebuild 或 Settings 触发 | P0 |
| **无 preview helper** | GA/Task/Result 预览需新建 | P0 |
| **无 model join query** | Model Selector Results 页需手写 SQL 或新 helper | P0 |
| **无单条 get helper** | GA version / task / model / batch 详情页缺便捷 API | P1 |
| **无 dashboard 聚合查询** | 需 COUNT 或 `get_registry_summary()` | P1 |
| **ChemDBWebVersion 级 requirements** | 仅有 `ChemDB/requirements.txt`，无 streamlit | P1 |
| **runs / step_executions 无前端页** | Phase 3A 可不做；Dashboard 卡片级展示即可 | P2 |

### A.3 与 Phase 3 总体设计的关系

[frontend.md](./frontend.md) 描述的是完整 Phase 3（含 Model Builder、Run Model Selection、Recommendation）。**Phase 3A 是子集**：只展示已有资源，唯一可选写操作是 **Rebuild SQLite Index**（调用现有 `app.storage.cli rebuild-index`）。

Phase 3A **明确不做**：启动 pipeline、训练、编辑 GA、创建 task、运行 evaluation、修改 `ChemDB/src`。

---

## 1. `repositories.py` 现有查询函数审计

### 1.1 现有函数（12 个）

| 函数 | 签名 | 返回 |
|------|------|------|
| `list_projects` | `(conn)` | 全部 projects |
| `list_runs` | `(conn, project_id=None)` | runs，可按 project 过滤 |
| `get_run` | `(conn, project_id, run_id)` | 单条 run 或 None |
| `list_ga_sets` | `(conn)` | 全部 ga_sets |
| `list_ga_versions` | `(conn, ga_set_id)` | 某 GA set 的全部 versions |
| `get_run_ga_binding` | `(conn, project_id, run_id)` | run 的 GA binding |
| `list_models` | `(conn, project_id=None, run_id=None)` | models，可过滤 |
| `list_model_after_tasks` | `(conn)` | 全部 tasks |
| `list_batches` | `(conn, task_id)` | 某 task 的全部 batches |
| `list_model_results` | `(conn, task_id, batch_id)` | 某 batch 的全部 model results（**无 join**） |
| `get_best_model_for_batch` | `(conn, task_id, batch_id)` | 最佳 result（按 batch.best_model_id 或 selection_score） |

### 1.2 各函数可支撑的页面

| 函数 | 可支撑页面 |
|------|-----------|
| `list_projects` + `list_runs` | Dashboard（概览）、Settings（workspace 信息）、Model Library（run 过滤） |
| `get_run` | Model Library 详情（run 状态）、Dashboard run 卡片 |
| `list_ga_sets` + `list_ga_versions` | **GA Library** 主表 + version 子表 |
| `get_run_ga_binding` | Model Library（模型来源 GA）、Dashboard（bound run 统计） |
| `list_models` | **Model Library** 主表 |
| `list_model_after_tasks` | **Task Library** 主表 |
| `list_batches` | **Model Selector Results**（按 task 分组 batch 列表） |
| `list_model_results` | **Model Selector Results** 排名表（缺 model GA 元数据） |
| `get_best_model_for_batch` | Dashboard 亮点卡片、Results 页默认选中 best |

### 1.3 缺失的只读查询函数

| 建议函数 | 用途 | 页面 |
|----------|------|------|
| `get_ga_set(conn, ga_set_id)` | GA set 详情 | GA Library |
| `get_ga_version(conn, ga_set_id, ga_version_id)` | version 详情 + `ga_csv_path` | GA Library |
| `get_model(conn, model_id)` | 单模型详情 | Model Library |
| `get_model_after_task(conn, task_id)` | task 详情 + 文件路径 | Task Library |
| `get_batch(conn, task_id, batch_id)` | batch 详情 + `summary_path` | Model Selector Results |
| **`list_model_results_with_models(conn, task_id, batch_id)`** | result LEFT JOIN models ON `registry_model_id` | **Model Selector Results**（核心） |
| `list_all_batches(conn)` | 跨 task 的 batch 总览 | Dashboard、Results 首页 |
| `get_registry_summary(conn)` | `{projects, runs, ga_sets, tasks, batches, models, results}` counts | **Dashboard** |
| `list_run_ga_bindings(conn)` | 全部 binding + 可选 join runs | Dashboard、GA Library（usage） |
| `list_models_with_run(conn, ...)` | models JOIN runs（status、ga_binding_status） | Model Library 增强列 |
| `list_runs_by_ga_version(conn, ga_set_id, ga_version_id)` | 哪些 run 绑定了某 GA version | GA Library 详情 |

**Phase 3A 最小必需新增（P0）**：`get_registry_summary`、`get_ga_version`、`get_model_after_task`、`get_batch`、`list_model_results_with_models`。

---

## 2. SQLite Schema 对前端是否够用

### 2.1 表级评估

| 表 | Phase 3A 是否够用 | 说明 |
|----|------------------|------|
| `projects` | ✅ | Dashboard 计数；Settings 展示 workspace root |
| `runs` | ✅ | Model Library 过滤；Dashboard 状态卡片 |
| `ga_sets` | ✅ | GA Library 主表 |
| `ga_versions` | ✅ | GA Library version 列表；`ga_csv_path` 供 preview |
| `run_ga_bindings` | ✅ | 展示 run↔GA 关系；models 已冗余 GA 字段 |
| `models` | ✅ | Model Library 主表；join 目标 |
| `model_after_tasks` | ✅ | Task Library；`task_dir` + 相对文件名 |
| `model_after_batches` | ✅ | Results 分组；`summary_path` 可 preview batch 汇总 |
| `model_after_model_results` | ✅ | Results 排名；`ranking_path`/`metrics_path`/`report_path` |
| `step_executions` | ⚠️ 可选 | Phase 3A 可不展示；Settings/Debug 扩展用 |
| `artifacts` | ⚠️ 可选 | 80 条 artifact 登记；Phase 3A 可不展示 |

**结论：schema 对 Phase 3A 只读 dashboard 足够，无需改 schema。**

### 2.2 页面-表-字段映射（主表 vs 详情）

#### Dashboard

| 数据源 | 主表/卡片字段 | 详情/跳转 |
|--------|--------------|----------|
| 聚合 COUNT | projects、runs、ga_sets、tasks、batches、models | 各 Library 页 |
| `runs` | `status`, `ga_binding_status`, `pipeline_core_ok`, `pipeline_training_ok` | Run 详情（可选） |
| `get_best_model_for_batch` × N | task_id、best MRR | Model Selector Results |

#### GA Library

| 表 | 主表列 | 详情页 |
|----|--------|--------|
| `ga_sets` | `ga_set_id`, `name`, `created_at`, `updated_at` | `description`, `tags_json`, `source_json`, `ga_set_root` |
| `ga_versions` | `ga_version_id`, `num_ga`, `checksum`, `created_at` | `parent_version_id`, `created_by`, `notes`, `source_json`, **`ga_csv_path`** → preview |
| `run_ga_bindings`（可选） | 使用该 version 的 run 数 | `project_id`, `run_id`, `bound_at`, `status` |

#### Task Library

| 表 | 主表列 | 详情页 |
|----|--------|--------|
| `model_after_tasks` | `task_id`, `metal`, `embedding`, `candidate_count`, `positive_count`, `negative_count` | `notes`, `task_dir`, 文件 preview |
| 文件（非 DB） | — | `task_dir/candidate_file` 等 → preview 前 100 行 |

#### Model Library

| 表 | 主表列 | 详情页 |
|----|--------|--------|
| `models` | `model_id`, `project_id`, `run_id`, `checkpoint_stem`, `val_mrr`, `ga_set_id`, `ga_version_id` | `checkpoint_path`, `config_path`, `embedding_backend`, `val_recall_at_*`, `ga_binding_checksum`, `size_bytes`, `history_path` |
| `runs`（join） | `status`, `ga_binding_status` | `run_root`, `last_finished_at` |
| `run_ga_bindings`（join） | binding `status`, `checksum` | `source_ga_csv`, `materialized_path` |

#### Model Selector Results

| 表 | 主表列 | 详情页 |
|----|--------|--------|
| `model_after_batches` | `task_id`, `batch_id`, `command`, `n_models`, `n_success`, `finished_at`, `best_model_id` | `output_dir`, `model_runs_file`, `summary_path` |
| `model_after_model_results` + **join models** | `rank_among_models`, `registry_model_id`, `mrr`, `hit_at_5`, `selection_score`, `status` | `ranking_path`, `metrics_path`, `report_path` + preview |
| join `models` | `ga_set_id`, `ga_version_id`, `val_mrr` | 完整 model 记录 |

#### Settings

| 来源 | 字段 |
|------|------|
| 环境/配置 | `workspace_root`, `db_path` |
| `connect_db` 探测 | DB 是否存在、表行数、`ingest_status` 异常计数 |
| 可选 action | Rebuild Index（调用 CLI，非 schema 写） |

---

## 3. `registry_model_id` Join 逻辑

### 3.1 设计确认

- **正确 join 键**：`model_after_model_results.registry_model_id = models.model_id`
- **不要用**：`model_after_model_results.model_id`（raw id，如 `run_demo`、`run_demo_last`）直接 join `models`
- ingest 已在 `_resolve_registry_model_id()` 写入 `registry_model_id`（合成 id 如 `p001_run_demo_best`）
- Phase 2A 验证：4/4 results 均为 `matched_model`

### 3.2 repositories 现状

| 项 | 状态 |
|----|------|
| join helper | ❌ **不存在** |
| `list_model_results` | 仅 `SELECT * FROM model_after_model_results`，无 join |
| `get_best_model_for_batch` | 用 batch 的 `best_model_id` 匹配 result 的 **`model_id`**（raw），非 `registry_model_id`；当前数据可用，但不返回 model GA 元数据 |

### 3.3 建议新增 helper

```python
def list_model_results_with_models(conn, task_id, batch_id) -> List[Dict]:
    """
    LEFT JOIN models m ON r.registry_model_id = m.model_id
    返回字段：
      - result 全字段（或别名 r_*）
      - m.ga_set_id, m.ga_version_id, m.ga_binding_checksum
      - m.val_mrr, m.checkpoint_stem, m.embedding_backend
      - join_status: 'matched_model' | 'missing_model'
    """
```

可选：`get_model_result_with_model(conn, task_id, batch_id, model_id)` 用于单条详情。

---

## 4. 文件 Preview 需求

### 4.1 现状

**app/ 下无任何 preview / read_csv_head helper**。

### 4.2 各类型 preview 路径解析

| 类型 | DB 字段 | 实际路径 | 读取方式 |
|------|---------|----------|----------|
| **GA** | `ga_versions.ga_csv_path` | ingest 写入绝对路径，如 `.../ga_sets/fixture_run_demo/versions/v001/GA_with_id.csv` | `pandas.read_csv(path, nrows=100)` 或 csv 模块 |
| **Task candidates** | `task_dir` + `candidate_file` | `task_dir/candidates.csv`（DB 存相对名） | 同上 |
| **Task positives/negatives** | `task_dir` + `positive_file` / `negative_file` | 同上；`negative_file` 可能为空 | 空则跳过 |
| **Result ranking** | `model_after_model_results.ranking_path` | 绝对路径；列含 rank, score, candidate_id, smiles 等 | nrows=100 |
| **Result metrics** | `metrics_path` | JSON 文件 | `json.load` 全量（通常小） |
| **Result report** | `report_path` | JSON 文件 | 全量或格式化展示 |
| **Batch summary** | `model_after_batches.summary_path` | `model_selection_summary.csv` | nrows=100 |

### 4.3 建议新增 helper（建议放 `app/storage/previews.py`）

| 函数 | 职责 |
|------|------|
| `resolve_task_file(task_row, which='candidates')` | `Path(task_dir) / candidate_file` |
| `read_csv_preview(path, limit=100)` | 返回 `List[Dict]` 或 DataFrame；文件不存在返回 error dict |
| `read_json_preview(path)` | 读 JSON；过大时可截断 |
| `preview_ga_version(conn, ga_set_id, ga_version_id)` | 查 DB 取 `ga_csv_path` → `read_csv_preview` |
| `preview_task(task_row, limit=100)` | candidates/positives/negatives 三块 preview |
| `preview_model_result(result_row, limit=100)` | ranking CSV + metrics/report JSON |

**约束**：preview 只读、默认 100 行、路径必须来自 DB 已 ingest 的绝对路径（不拼接用户输入防 path traversal）。

---

## 5. Streamlit 目录现状

### 5.1 现状

- **`app/streamlit_ui/` 不存在**
- `app/` 仅有 `services/`、`storage/`、`templates/`
- `initial_assessment.md` 曾规划 `streamlit_ui/`，未实现

### 5.2 建议创建的文件结构

```text
app/streamlit_ui/
├── app.py                      # 入口：st.set_page_config, 导航, DB 连接
├── state.py                    # session_state: db_path, workspace_root, conn cache
├── pages/
│   ├── 01_Dashboard.py
│   ├── 02_GA_Library.py
│   ├── 03_Task_Library.py
│   ├── 04_Model_Library.py
│   ├── 05_Model_Selector_Results.py
│   └── 06_Settings.py
└── components/
    ├── sidebar.py              # 全局 sidebar、DB 状态、Rebuild 入口
    ├── tables.py               # st.dataframe 封装、空态
    ├── cards.py                # metric 卡片、summary 行
    └── previews.py             # 调用 storage/previews.py 渲染
```

可选薄封装：`app/streamlit_ui/db_session.py` — 包装 `connect_db` + repository 调用。

---

## 6. 启动方式与依赖

### 6.1 环境探测（scidata）

| 包 | 状态 |
|----|------|
| **pandas** | ✅ 2.0.3 |
| **streamlit** | ❌ `ModuleNotFoundError` |
| **sqlite3** | ✅ 标准库 |
| **rdkit** | ✅（GA preview 不需要；Phase 3A 只读 CSV 即可） |

### 6.2 当前 DB 状态

- **`workspace/chemdb.sqlite` 不存在**（gitignore，需本地 rebuild）
- 验证用临时 DB：`python -m app.storage.cli rebuild-index --workspace-root workspace --db-path workspace/chemdb.sqlite`

### 6.3 推荐启动命令（Phase 3A 实施后）

```bash
cd ChemDBWebVersion
conda activate scidata
pip install streamlit   # 或 conda install streamlit
python -m app.storage.cli rebuild-index   # 首次或 workspace 变更后
PYTHONPATH=. streamlit run app/streamlit_ui/Home.py
```

环境变量建议：

```bash
export CHEMDB_WORKSPACE_ROOT=/path/to/workspace
export CHEMDB_DB_PATH=/path/to/workspace/chemdb.sqlite
```

### 6.4 Requirements 更新建议

- 新增 **`ChemDBWebVersion/requirements-web.txt`**（或 `requirements.txt`）：

```text
streamlit>=1.28
pandas>=2.0
```

- **不修改** `ChemDB/requirements.txt`（科研 pipeline 依赖隔离）
- 在 `TESTING.md` 增加 Phase 3A smoke test 说明

---

## B. 缺失函数清单（汇总）

### Repository 层（P0–P1）

1. `get_registry_summary(conn)` — P0
2. `get_ga_set(conn, ga_set_id)` — P1
3. `get_ga_version(conn, ga_set_id, ga_version_id)` — P0
4. `get_model(conn, model_id)` — P1
5. `get_model_after_task(conn, task_id)` — P0
6. `get_batch(conn, task_id, batch_id)` — P0
7. **`list_model_results_with_models(conn, task_id, batch_id)`** — P0
8. `list_all_batches(conn)` — P1
9. `list_runs_by_ga_version(conn, ga_set_id, ga_version_id)` — P2

### Preview 层（P0）

1. `read_csv_preview(path, limit=100)`
2. `read_json_preview(path)`
3. `resolve_task_file(task_row, which)`
4. `preview_ga_version(...)` / `preview_task(...)` / `preview_model_result(...)`

### Streamlit 基础设施（P0）

1. 整个 `app/streamlit_ui/` 目录
2. DB 缺失 / stale 检测与 Settings rebuild 按钮

---

## C. 推荐文件结构

```text
ChemDBWebVersion/
├── requirements-web.txt          # 新增：streamlit + pandas pin
├── app/
│   ├── storage/
│   │   ├── repositories.py       # 修改：新增查询
│   │   └── previews.py           # 新增：文件 preview
│   └── streamlit_ui/             # 新增：整目录
│       ├── app.py
│       ├── state.py
│       ├── pages/ (6 pages)
│       └── components/ (4 modules)
├── tests/
│   ├── storage/
│   │   ├── test_phase2a_sqlite_registry.py  # 修改：repository 测试
│   │   └── test_previews.py                 # 新增
│   └── streamlit_ui/
│       └── test_import_smoke.py             # 新增
└── scripts/
    └── run_streamlit.sh                     # 可选：启动包装
```

**不创建 / 不修改**：`ChemDB/src/**`、`app/services/orchestrator.py` pipeline 语义、`schema.sql`。

---

## D. 页面-数据源映射表

| 页面 | 主查询 | 辅助查询 | 文件 Preview | 写操作 |
|------|--------|----------|--------------|--------|
| **01 Dashboard** | `get_registry_summary` | `list_runs`, `list_model_after_tasks`, `get_best_model_for_batch` | 无 | 无 |
| **02 GA Library** | `list_ga_sets`, `list_ga_versions` | `get_ga_version`, `list_runs_by_ga_version` | `ga_csv_path` → GA CSV | 无 |
| **03 Task Library** | `list_model_after_tasks` | `get_model_after_task`, `list_batches(task_id)` | task_dir 下 CSV | 无 |
| **04 Model Library** | `list_models` | `get_model`, `get_run`, `get_run_ga_binding` | 无（或 history.json 可选） | 无 |
| **05 Model Selector Results** | `list_all_batches` 或 task→batch 两级 | `list_model_results_with_models`, `get_batch` | ranking/metrics/report | 无 |
| **06 Settings** | `connect_db` 探测、表 COUNT | — | — | **可选** `rebuild-index` |

---

## E. Phase 3A 实施计划

### E.1 建议新增文件

| 文件 | 内容 |
|------|------|
| `requirements-web.txt` | streamlit, pandas |
| `app/storage/previews.py` | CSV/JSON preview helpers |
| `app/streamlit_ui/Home.py` | 入口、page config |
| `app/streamlit_ui/state.py` | session state、路径默认值 |
| `app/streamlit_ui/pages/01_Dashboard.py` | 计数卡片 + 最近 batch 亮点 |
| `app/streamlit_ui/pages/02_GA_Library.py` | GA set/version 表 + preview |
| `app/streamlit_ui/pages/03_Task_Library.py` | Task 表 + CSV preview |
| `app/streamlit_ui/pages/04_Model_Library.py` | Model 表 + GA/run 列 |
| `app/streamlit_ui/pages/05_Model_Selector_Results.py` | Batch 选择 → 排名表 + preview |
| `app/streamlit_ui/pages/06_Settings.py` | 路径配置、DB 状态、Rebuild |
| `app/streamlit_ui/components/*.py` | sidebar、tables、cards、previews UI |
| `tests/storage/test_previews.py` | preview 单元测试 |
| `tests/storage/test_repositories_phase3a.py` | 新 repository 测试 |
| `tests/streamlit_ui/test_import_smoke.py` | 各 page import smoke |
| `scripts/run_streamlit.sh` | 可选启动脚本 |

### E.2 建议修改文件

| 文件 | 修改 |
|------|------|
| `app/storage/repositories.py` | 新增 §B 中 P0/P1 查询 |
| `TESTING.md` | Phase 3A 测试与启动说明 |
| `_auditing/README.md` | 索引 Phase 3A audit |

### E.3 实施顺序（建议 4 步）

1. **Repository + Preview 层**（可 pytest，无 streamlit）
2. **Streamlit 骨架 + Settings + DB 连接**（处理 DB 不存在）
3. **6 个只读页面**（先列表后 preview）
4. **Smoke test + 真实 workspace 手工验收**

### E.4 明确边界

| 做 | 不做 |
|----|------|
| 只读展示 + preview | 启动 pipeline / training / evaluation |
| Settings 中 Rebuild Index | 编辑 GA / 创建 task / bind GA |
| 调用 `app.storage.cli` | 修改 `ChemDB/src` |
| `connect_db(read_only=True)` 用于展示 | 直接 SQL 写 registry |

---

## F. Phase 3A 验收标准

### F.1 Repository / Preview 测试

| 测试 | 断言 |
|------|------|
| `test_get_registry_summary` | rebuild 后 counts 与 workspace 一致（1/1/1/2/3/2/4） |
| `test_list_model_results_with_models` | join 后含 `ga_set_id`；4 条均为 matched |
| `test_read_csv_preview_limit` | 大文件只返回 100 行 |
| `test_preview_ga_version` | 路径来自 `ga_csv_path` |
| `test_preview_task` | 解析 `task_dir + candidate_file` |

### F.2 Streamlit smoke

| 测试 | 方法 |
|------|------|
| import smoke | `import app.streamlit_ui.pages.02_GA_Library` 等不抛错 |
| 无 streamlit 环境 | CI 可 `@pytest.mark.skipif` 无 streamlit；或仅测 storage 层 |

### F.3 连接真实 DB

```bash
python -m app.storage.cli rebuild-index --workspace-root workspace --db-path workspace/chemdb.sqlite
PYTHONPATH=. streamlit run app/streamlit_ui/Home.py
```

| 检查 | 预期 |
|------|------|
| Settings 显示 DB 路径、表行数 | 与 rebuild summary 一致 |
| DB 不存在时 | 友好提示 + Rebuild 按钮，不 crash |
| 只读连接 | 展示阶段用 `read_only=True` |

### F.4 各页面加载（真实 workspace）

| 页面 | 验收 |
|------|------|
| Dashboard | 显示 1 project、1 run、1 GA set、2 tasks、3 batches、2 models |
| GA Library | `fixture_run_demo` / v001，126 GA，preview 有 SMILES 列 |
| Task Library | `ni_lkb_p_ni`（344 candidates）、`pd_lkb_p_cluster`，preview 可切换 |
| Model Library | `p001_run_demo_best/last`，含 GA metadata |
| Model Selector Results | `batch_dummy_ckpts` 排名 2 模型；join 显示 GA；preview ranking |
| Settings | Rebuild 后行数不变（幂等） |

### F.5 非功能验收

- 全程无写 pipeline / 无修改 workspace 数据文件（Rebuild 仅更新 sqlite）
- 无 `ChemDB/src` diff
- pytest 全绿（storage 新增 + 既有 62 passed）
- 页面加载 < 3s（本地 sqlite + 100 行 preview）

---

## 附录：真实数据快照（rebuild-index 2026-05-28）

| 实体 | 数量 | 示例 ID |
|------|------|---------|
| projects | 1 | p001 |
| runs | 1 | run_demo |
| ga_sets | 1 | fixture_run_demo |
| ga_versions | 1 | v001（126 GA） |
| model_after_tasks | 2 | ni_lkb_p_ni, pd_lkb_p_cluster |
| model_after_batches | 3 | batch_dummy_ckpts, p001_run_demo ×2 |
| models | 2 | p001_run_demo_best, p001_run_demo_last |
| model_after_model_results | 4 | 全部 registry_model_id join 成功 |

---

**总结**：Phase 2A registry 已为 Phase 3A 提供足够的数据底座；主要工作是补 5–7 个 repository 函数、1 个 previews 模块、Streamlit 六页骨架、安装 streamlit、首次 rebuild 本地 sqlite。无需改 schema 或 pipeline，可按上述计划直接进入实现。
