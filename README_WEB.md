# ChemDBWebVersion Web MVP (Phase 3B)

本项目的 Web 端已从 Phase 3A 只读升级到 Phase 3B MVP，可在不修改 `ChemDB/src` 算法语义的前提下完成资源注入与异步任务执行。

## 启动

```bash
PYTHONPATH=. streamlit run app/streamlit_ui/Home.py
```

可选：使用脚本

```bash
./scripts/run_streamlit.sh
```

## 主要页面

- `Home`: MVP 总览与推荐操作流
- `01_Dashboard`: SQLite 资源统计
- `02_GA_Library`: GA 上传创建 version（严格 `GA_SMILES + GA_ID`）
- `03_Task_Library`: task 上传创建与预览
- `04_Model_Library`: 模型资源浏览
- `05_Model_Selector_Results`: 选择器结果回看（含 join 状态）
- `06_Settings`: 环境状态与 rebuild-index
- `07_Jobs_and_Logs`: 后台作业状态与日志
- `08_Model_Builder`: 异步模型构建链（init/bind/core/training/rebuild）
- `09_Model_Selector`: 异步 evaluate-models
- `10_Recommendation`: 推荐任务（existing/temporary task）+ ranking 下载
- `11_Hyperparameters`: 超参集上传与版本浏览

## Phase 3B 写入范围

- 仅写 `workspace/`：
  - `workspace/ga_sets/**`
  - `workspace/model_after_tasks/**`
  - `workspace/model_after_results/**`
  - `workspace/recommendation_results/**`
  - `workspace/hyperparameter_sets/**`
  - `workspace/jobs/**`
  - `workspace/chemdb.sqlite`

## 关键机制

- 长任务通过 `jobs` 表 + `workspace/jobs/*.log` 追踪，页面不阻塞等待。
- 每次资源注入或任务完成后都通过 `app.storage.cli rebuild-index` 回流 SQLite。
- `Model Selector Results` 继续通过 `registry_model_id -> models.model_id` 左连接展示模型元信息。
