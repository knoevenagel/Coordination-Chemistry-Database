# Phase 3B MVP Implementation

## Scope

- 目标：把 Phase 3A 只读前端推进到可执行 MVP。
- 边界：不改 `ChemDB/src/**` 与 `ChemDB/training/**` 算法语义，仅做路径/编排/调用封装。

## Delivered

- Stage 0: `Home.py` 入口与引用统一，保持 `PYTHONPATH=. streamlit run app/streamlit_ui/Home.py`。
- Stage 1: 新增 `jobs` schema、`job_manager`、`07_Jobs_and_Logs.py`、sidebar recent jobs。
- Stage 2: GA 上传创建 version（严格 `GA_SMILES + GA_ID`）并接入 `02_GA_Library.py`。
- Stage 3: task 上传创建服务 `task_registry_manager` 并接入 `03_Task_Library.py`。
- Stage 4: `model_build_manager` + `08_Model_Builder.py` 异步执行 init/bind/core/training/rebuild。
- Stage 5: `09_Model_Selector.py` 可执行 evaluate-models，结果回流 SQLite。
- Stage 6: `10_Recommendation.py` 支持 existing/temporary task 与 ranking 下载。
- Stage 7: `hyperparameter_manager`、超参资源库 schema/ingest/repository、`11_Hyperparameters.py`、Model Builder 白名单注入。
- Stage 8: 首页与侧边栏文案/流程提示统一，路径列持续后置展示。

## Testing

建议顺序：

```bash
./scripts/run_pytest.sh tests/storage -v
./scripts/run_pytest.sh tests/services -v
./scripts/run_pytest.sh tests/streamlit_ui -v
./scripts/run_pytest.sh -q
PYTHONPATH=. timeout 12s streamlit run app/streamlit_ui/Home.py --server.headless true --server.port 8510
```
