# Phase 3B Readiness Audit

## 审计结论（实施前）

基于 Phase 3A 状态，进入 Phase 3B 具备以下前提：

- 已有稳定 Streamlit 基线：`Home.py + pages + components + state`。
- 已有 SQLite 读模型与文件预览能力：`repositories.py`、`previews.py`。
- 已有核心 run-based 编排能力：`orchestrator.py`、`ga_registry_manager.py`、`model_after_manager.py`。
- 已验证 `workspace-root` 隔离语义（GA 绑定与 SQLite ingest）。

## 关键缺口（实施前）

- 无 jobs 表与标准异步状态机。
- GA/Task 仅可读，缺少上传创建入口。
- 无可执行 Model Builder / Model Selector / Recommendation 页面。
- 无 Hyperparameter 资源库与最小注入机制。

## Phase 3B 执行策略

- 写入统一限定在 `workspace/**`。
- 长任务统一走 jobs（subprocess + log file + SQLite）。
- 结果统一通过 `rebuild-index` 回流到 registry。
- UI 页面维持“读 SQLite + 读文件预览”风格，避免前端直接耦合算法实现。
