# Phase 3A 实施记录（Resource-Centric Read-Only Streamlit Dashboard）

| 字段 | 值 |
|------|-----|
| 类型 | 实施记录 |
| 日期 | 2026-05-29 |
| 阶段 | Phase 3A |
| 范围 | 只读资源组织与展示（Dashboard / GA Library / Task Library / Model Library / Model Selector Results / Settings） |
| 边界 | 不启动 pipeline；不训练模型；不编辑 GA；不创建 task；不运行 evaluation；不修改 `ChemDB/src/**` |
| 唯一允许写操作 | Settings 页面可选 `rebuild-index`（调用 `python -m app.storage.cli rebuild-index`） |
| 前置审计 | [phase_3a_readiness_audit.md](./phase_3a_readiness_audit.md) |

---

## 1. 实施目标与结果

### 1.1 目标

在现有 Phase 2A SQLite registry 基础上，新增 Streamlit 前端，提供资源中心化只读浏览能力：

- Dashboard
- GA Library
- Task Library
- Model Library
- Model Selector Results
- Settings

### 1.2 结果

Phase 3A 已落地，核心结果如下：

1. `app/storage/repositories.py` 补齐 Phase 3A 所需只读查询函数。
2. 新增 `app/storage/previews.py`，支持 GA/task/result 文件预览（CSV 前 100 行 + JSON）。
3. 新增 `app/streamlit_ui/` 完整目录与 6 个页面。
4. 新增 `requirements-web.txt` 与 `scripts/run_streamlit.sh`。
5. 新增 Phase 3A 测试并通过；全量 pytest 通过。
6. `PYTHONPATH=. streamlit run app/streamlit_ui/Home.py` 已做启动冒烟验证。

---

## 2. 代码变更清单

### 2.1 新增文件

- `app/storage/previews.py`
- `app/streamlit_ui/__init__.py`
- `app/streamlit_ui/Home.py`
- `app/streamlit_ui/state.py`
- `app/streamlit_ui/components/__init__.py`
- `app/streamlit_ui/components/sidebar.py`
- `app/streamlit_ui/components/cards.py`
- `app/streamlit_ui/components/tables.py`
- `app/streamlit_ui/components/previews.py`
- `app/streamlit_ui/pages/__init__.py`
- `app/streamlit_ui/pages/01_Dashboard.py`
- `app/streamlit_ui/pages/02_GA_Library.py`
- `app/streamlit_ui/pages/03_Task_Library.py`
- `app/streamlit_ui/pages/04_Model_Library.py`
- `app/streamlit_ui/pages/05_Model_Selector_Results.py`
- `app/streamlit_ui/pages/06_Settings.py`
- `requirements-web.txt`
- `scripts/run_streamlit.sh`
- `tests/storage/test_repositories_phase3a.py`
- `tests/storage/test_previews.py`
- `tests/streamlit_ui/test_import_smoke.py`

### 2.2 修改文件

- `app/storage/repositories.py`
  - 新增：
    - `get_registry_summary`
    - `get_ga_set`
    - `get_ga_version`
    - `get_model`
    - `get_model_after_task`
    - `get_batch`
    - `list_all_batches`
    - `list_model_results_with_models`（通过 `registry_model_id` join）
    - `list_runs_by_ga_version`
- `tests/storage/test_phase2a_sqlite_registry.py`
  - 调整 `test_no_write_to_real_workspace_or_repo` 断言为“真实 DB 存在状态前后不变”，避免与“先重建 DB 后跑全量”冲突。

---

## 3. 功能说明（页面）

### 3.1 Dashboard

- 展示资源计数卡片：Projects / Runs / GA Sets / GA Versions / Models / Tasks / Batches / Results
- 展示近期 runs/models/batches 表
- DB 不存在时仅提示，不崩溃

### 3.2 GA Library

- 展示 `ga_sets`
- 按 `ga_set_id` 展示 `ga_versions`
- 预览 `ga_csv_path` 指向的 GA CSV 前 100 行
- 展示使用该 GA version 的 run（`list_runs_by_ga_version`）

### 3.3 Task Library

- 展示 task 列表（metal/embedding/counts/notes）
- 展示 task 详情
- 预览 candidates / positives / negatives
- 展示 task 对应 batches

### 3.4 Model Library

- 支持 project/run 过滤
- 展示 model 主表（含 GA 元数据）
- 展示 model 详情（checkpoint/config/index/history 等路径）

### 3.5 Model Selector Results

- 选择 task + batch
- 展示 batch metadata
- 展示结果排序表（必须含 `model_join_status`）
- 使用 `registry_model_id = models.model_id` 做 LEFT JOIN
- 预览 ranking.csv / metrics.json / report.json

### 3.6 Settings

- 展示 workspace_root、db_path、DB 存在状态、Python 信息
- 展示 SQLite 汇总计数
- 提供 Rebuild 按钮（调用 `python -m app.storage.cli rebuild-index ...`）

---

## 4. 测试与验收

### 4.1 依赖安装

```bash
pip install -r requirements-web.txt
conda run -n scidata pip install -r requirements-web.txt
```

### 4.2 DB 重建

```bash
python -m app.storage.cli rebuild-index --workspace-root workspace --db-path workspace/chemdb.sqlite
```

返回 summary：

- projects: 1
- runs: 1
- ga_sets: 1
- tasks: 2
- batches: 3
- errors: []

### 4.3 新增测试

```bash
./scripts/run_pytest.sh tests/storage/test_repositories_phase3a.py -v
./scripts/run_pytest.sh tests/storage/test_previews.py -v
./scripts/run_pytest.sh tests/streamlit_ui/test_import_smoke.py -v
```

结果：

- `test_repositories_phase3a.py`: 4 passed
- `test_previews.py`: 7 passed
- `test_import_smoke.py`: 2 passed

### 4.4 全量测试

```bash
./scripts/run_pytest.sh -q
```

结果：

- 75 passed, 12 deselected

### 4.5 Streamlit 启动冒烟

```bash
timeout 12s PYTHONPATH=. streamlit run app/streamlit_ui/Home.py --server.headless true --server.port 8510
```

结果：

- 服务成功启动并输出 Local URL / Network URL
- `timeout` 触发停止（预期）

---

## 5. 边界符合性复核

本阶段实现满足边界约束：

- 未修改 `ChemDB/src/**`
- 未引入 pipeline/training/evaluation 操作入口
- 页面均为只读展示
- 唯一可写操作仅 Settings 的 `rebuild-index`

---

## 6. 已知限制（Phase 3A 内可接受）

1. 暂无复杂图表，仅表格与基础卡片。
2. 过滤器当前以 project/run 为主（Model Library）。
3. 结果页为既有 batch 浏览，不支持创建新评估（符合 Phase 3A 范围）。

---

## 7. 下一步建议（Phase 3B+）

1. 补充更细粒度筛选与排序（GA/task/model/results）。
2. 增加页面级缓存与大表分页。
3. 在保持边界前提下，逐步规划 Model Builder / Selector 执行面板（后续阶段）。

---

## 8. 2026-05-29 启动修复（Hotfix）

### 问题

`streamlit run app/streamlit_ui/app.py` 时出现：

`ModuleNotFoundError: No module named 'app.streamlit_ui'; 'app' is not a package`

根因是入口文件 `app.py` 与顶层 `app` 包同名，Streamlit 运行时导入发生遮蔽。

### 修复

1. 入口文件重命名：`app/streamlit_ui/app.py` -> `app/streamlit_ui/Home.py`
2. 启动脚本更新为：
   `PYTHONPATH=. streamlit run app/streamlit_ui/Home.py`
3. 测试引用更新：`import app.streamlit_ui.app` -> `import app.streamlit_ui.Home`
4. 文档中的启动命令与路径同步改为 `Home.py`

### 复测

```bash
./scripts/run_pytest.sh tests/streamlit_ui/test_import_smoke.py -v
./scripts/run_pytest.sh -q
timeout 12s env PYTHONPATH=. streamlit run app/streamlit_ui/Home.py --server.headless true --server.port 8510
```

结果：

- import smoke：2 passed
- 全量：75 passed, 12 deselected
- Streamlit 冒烟：服务成功启动并输出 `http://localhost:8510`
