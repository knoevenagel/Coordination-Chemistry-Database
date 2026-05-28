# Phase 1E-MVP 实施记录：Workspace GA + Run Binding

**日期**：2026-05-21  
**状态**：已完成（MVP）

## 目标

GA 从 core pipeline 自动中间产物改为 **workspace 级可版本化输入**；每个 run **绑定** GA version 并 **物化** 到 `run/tmp/GA_with_id.csv` 后，core 才继续 `apply-ga` 与 step6_7+。

## 实现摘要

| 组件 | 变更 |
|------|------|
| `ChemDB/src/step4_5.py` | `--mode full-auto \| generate-ga \| apply-ga` |
| `app/services/ga_registry_manager.py` | GA workspace CRUD + bind + apply CLI |
| `app/services/step_registry.py` | core 顺序：`require_bound_ga` → `apply_ga_to_run`；`step4_5` 保留为 legacy |
| `app/services/orchestrator.py` | in-process `require_bound_ga`；`--through step4_5` 明确报错 |
| 测试 fixture | `tests/run_isolation/fixtures/ga_sets/{fixture_run_demo,fixture_restructured}/versions/v001/` |
| 1% integration | Tier B/C：step1–2 → bind fixture → apply-ga → downstream |

## 新 core pipeline

```text
step1 → step2 → require_bound_ga → apply_ga_to_run → step6_7 → step8 → step12 → step13
```

## GA workspace 布局

```text
workspace/ga_sets/{ga_set_id}/ga_set.json
workspace/ga_sets/{ga_set_id}/versions/{ga_version_id}/GA_with_id.csv
workspace/ga_sets/{ga_set_id}/versions/{ga_version_id}/ga_version.json
```

Run 绑定：`run/manifests/ga_binding.json` + 物化 `run/tmp/GA_with_id.csv`。

## CLI 示例

```bash
# 从 run 配体库生成 workspace GA version（authoring）
python -m app.services.ga_registry_manager generate-ga-from-run \
  --run-root workspace/projects/p001/runs/run_demo \
  --ga-set-id p001_run_demo_derived

# 绑定已有 version 到 run
python -m app.services.ga_registry_manager bind-ga-version-to-run \
  --run-root workspace/projects/p001/runs/run_demo \
  --ga-set-id fixture_run_demo --ga-version-id v001

# 对已绑定 run 计算 GAC/IRL
python -m app.services.ga_registry_manager apply-ga-to-run \
  --run-root workspace/projects/p001/runs/run_demo

# Core pipeline（需先 bind）
python -m app.services.orchestrator run-pipeline \
  --run-root workspace/projects/p001/runs/run_demo \
  --pipeline core --through step6_7
```

`run-pipeline --through step4_5` **已拒绝**，提示使用 GA registry。

## 测试 fixture GA 版本

| ga_set_id | 来源 | 约 num_ga |
|-----------|------|-----------|
| `fixture_run_demo` | `workspace/.../run_demo/tmp/GA_with_id.csv` | 126 |
| `fixture_restructured` | `ChemDB_restructured/ChemDB/tmp/GA_with_id.csv` | 727 |

## 未做（MVP 外）

- SQLite / Streamlit
- `create-ga-version-from-csv` 已在代码中实现（含校验），集成测试覆盖
- training / model_after 逻辑未改；stale 报告含 embedding、ckpt、model_after_results 路径模式
- orchestrator 不因 stale 硬阻断（仅 `ga_stale_report.json`）

## pytest

```bash
./scripts/run_pytest.sh tests/run_isolation/test_phase1e_ga_workspace.py -v
./scripts/run_pytest.sh -q   # 49 passed + 12 deselected integration
```
