# ChemDBWebVersion 审计与变更记录

本目录存放 Run-based 改造相关的**评估、实施记录与后续变更日志**。每次阶段性落地或重要补丁，应在此新增或更新文档，便于回溯「做了什么、为什么、验收如何」。

## 文档索引

| 文档 | 类型 | 说明 |
|------|------|------|
| [initial_assessment.md](./initial_assessment.md) | 评估 | 改造前源码路径扫描、风险与 Phase 规划（基线） |
| [phase_1a_implementation.md](./phase_1a_implementation.md) | 实施记录 | Phase 1A 落地过程、文件清单、约束与验收 |
| [phase_1b_implementation.md](./phase_1b_implementation.md) | 实施记录 | Phase 1B：1% PubChem 子集、integration pytest、Tier A/B |
| [phase_1c_plan.md](./phase_1c_plan.md) | 实施计划 | Phase 1C：ECFP embedding + per-run training 路径隔离 |
| [phase_1c_implementation.md](./phase_1c_implementation.md) | 实施记录 | Phase 1C 落地、测试与已知限制 |
| [phase_1abc_overall_assessment.md](./phase_1abc_overall_assessment.md) | **整体评估** | Phase 1A–1C 审计、run_demo 产物、模型后现状与 Phase 1D/2/3 建议 |
| [phase_1d_plan.md](./phase_1d_plan.md) | 实施计划 | Phase 1D：model-after run-based wrapper |
| [phase_1d_implementation.md](./phase_1d_implementation.md) | **实施记录** | Phase 1D：加载/打分/评估/多模型选择、CLI 与测试 |
| [phase_2_readiness_after_ga_workspace.md](./phase_2_readiness_after_ga_workspace.md) | **评估** | Phase 2 就绪审计（post Phase 1E） |
| [phase_2a_implementation.md](./phase_2a_implementation.md) | **实施记录** | Phase 2A：SQLite read-only registry + ingest CLI |
| [phase_2a_validation_output.md](./phase_2a_validation_output.md) | 验证输出 | Phase 2A `rebuild-index` 对真实 workspace 的完整终端记录 |
| [phase_2a_bound_ga_validation_output.md](./phase_2a_bound_ga_validation_output.md) | 验证输出 | Phase 2A bound-GA run + 复制 workspace 验证 |
| [frontend.md](./frontend.md) | 设计草案 | Phase 3 resource-centric 前端总体设计 |
| [phase_3a_readiness_audit.md](./phase_3a_readiness_audit.md) | **评估** | Phase 3A：只读 Streamlit dashboard 就绪审计 |
| [phase_3a_implementation.md](./phase_3a_implementation.md) | **实施记录** | Phase 3A：只读 Streamlit dashboard 落地、测试与验收 |
| [phase_3a_end_to_end_dataflow.md](./phase_3a_end_to_end_dataflow.md) | 架构说明 | 从原始数据到 SQLite ingest 再到网页展示的完整脚本/文件依赖链路 |
| [phase_3b_readiness_audit.md](./phase_3b_readiness_audit.md) | **评估** | Phase 3B：MVP 可写化与异步执行就绪审计 |
| [phase_3b_implementation.md](./phase_3b_implementation.md) | **实施记录** | Phase 3B：GA/Task 上传、Jobs、Model Builder、Selector、Recommendation、Hyperparameters |

## 阶段路线图（简）

| 阶段 | 内容 |
|------|------|
| Phase 1A | Run 隔离 + 编排层 + step 路径补丁 |
| Phase 1B | integration 数据集 + 真实 step 测试 |
| Phase 1C | ECFP embedding + per-run training/index/config/data/train（见 phase_1c_implementation.md） |
| Phase 1D | model-after：加载/打分/评估/多 model_run 对比（见 phase_1d_implementation.md） |
| Phase 1E | workspace GA sets + run binding（见 phase_1e_implementation.md） |
| Phase 2A | SQLite read-only registry + ingest CLI（见 phase_2a_implementation.md） |
| Phase 3A | 只读 Streamlit dashboard（见 phase_3a_readiness_audit.md、phase_3a_implementation.md、frontend.md） |
| Phase 3B | MVP 可写化 + 异步作业 + 推荐 + 超参数资源库 |
| 之后 | 可选并发隔离、任务调度增强、`--write-db` hook 等 |

## 后续步骤如何记档（约定）

1. **默认流程（从现在开始）**：每个阶段至少维护两份文档：`phase_<id>_readiness_audit.md`（阶段开始前审计）+ `phase_<id>_implementation.md`（阶段完成后实现记录），并在本文「文档索引」补齐链接。
2. **新阶段**：新增实施记录（`phase_<id>_implementation.md`，或 `changelog_YYYY-MM-DD_<topic>.md`）。
3. **每条记录建议包含**：
   - 日期、阶段 ID、关联需求/对话要点
   - **目标**与**不在范围内**的事项
   - **目录/命名**变更（若有）
   - **新增文件**与**修改文件**列表（路径 + 一行说明）
   - **行为变更**（CLI、默认路径、测试策略）
   - **验收**：命令、通过条件、已知限制
4. **评估文档**：`initial_assessment.md` 仅在基线事实变化时修订；阶段结论写在对应 `phase_*` 文档，避免单文件无限膨胀。
5. **与代码同步**：`ChemDB`、`app/services`、`tests/`、`pytest.ini` 中的路径名以 **`ChemDB`** 为准；勿再引入 `ChemDB_run_based` 作为默认名。

## 测试

所有 pytest **默认须在 `scidata` conda 环境**中运行，见 [TESTING.md](../TESTING.md) 与 `./scripts/run_pytest.sh`。

## 仓库布局（Phase 1A 后）

```
ChemDBWebVersion/
├── ChemDB/                 # 流水线源码（原 ChemDB_run_based 已重命名）
├── app/services/           # RunContext、StepRegistry、Orchestrator
├── tests/run_isolation/    # 默认轻量测试；integration 需 -m integration
├── workspace/              # 真实 run 工作区（git 仅保留占位）
└── _auditing/              # 本目录
```

## 环境变量

| 变量 | 默认 |
|------|------|
| `CHEMDB_REPO_ROOT` | `ChemDBWebVersion/ChemDB` |
| `CHEMDB_RUN_ROOT` | 由 orchestrator 在子进程中设为当前 `--run-root` |
