# Phase 3B Implementation Record

## 目标

将 Phase 3A 只读 dashboard 扩展为可执行 MVP，同时保持科研算法语义不变。

## 已落地项

1. **启动与入口收口**
   - `Home.py` 入口稳定，脚本/测试统一引用。

2. **Jobs 基础设施**
   - `schema.sql` 新增 `jobs` 表。
   - `job_manager.py` 提供 create/mark/run 与 async worker。
   - `repositories.py` 增加 jobs 查询函数。
   - `07_Jobs_and_Logs.py` + sidebar recent jobs。

3. **GA / Task 可写化**
   - `ga_registry_manager.create_ga_set_or_version_from_upload`（严格 CSV）。
   - `task_registry_manager.create_task_from_uploads`。
   - 页面 `02_GA_Library.py`、`03_Task_Library.py` 加上传创建入口。

4. **可执行模型流程**
   - `model_build_manager.py` + `08_Model_Builder.py`（异步作业链）。
   - `model_after_manager` 增加 evaluate/recommend job 提交封装。
   - `09_Model_Selector.py` 发起 evaluate-models job。
   - `10_Recommendation.py` 支持 existing/temporary task 与下载。

5. **Hyperparameter 最小生效**
   - `hyperparameter_manager.py`（上传 version + apply-to-run 白名单注入）。
   - SQLite 新增 `hyperparameter_sets` / `hyperparameter_versions`。
   - `ingest.py`、`repositories.py` 增补超参资源 ingest/query。
   - `11_Hyperparameters.py` 页面落地。

6. **UI 收尾**
   - Home/Sidebar 文案升级到 Phase 3B 流程。
   - 表格路径列保持后置展示。

## 新增/修改的核心文件

- 新增：
  - `app/services/job_manager.py`
  - `app/services/task_registry_manager.py`
  - `app/services/model_build_manager.py`
  - `app/services/hyperparameter_manager.py`
  - `app/streamlit_ui/pages/07_Jobs_and_Logs.py`
  - `app/streamlit_ui/pages/08_Model_Builder.py`
  - `app/streamlit_ui/pages/09_Model_Selector.py`
  - `app/streamlit_ui/pages/10_Recommendation.py`
  - `app/streamlit_ui/pages/11_Hyperparameters.py`
  - `README_WEB.md`
  - `docs/phase_3b_mvp_implementation.md`
- 修改：
  - `app/storage/schema.sql`
  - `app/storage/ingest.py`
  - `app/storage/repositories.py`
  - `app/services/ga_registry_manager.py`
  - `app/services/model_after_manager.py`
  - `app/streamlit_ui/pages/02_GA_Library.py`
  - `app/streamlit_ui/pages/03_Task_Library.py`
  - `app/streamlit_ui/Home.py`
  - `app/streamlit_ui/components/sidebar.py`

## 已知限制

- 现阶段 jobs 通过本地后台 worker 执行，不含分布式调度与并发隔离。
- Hyperparameter 仅白名单字段生效，非白名单字段会记录为 ignored。
