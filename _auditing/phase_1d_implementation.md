# Phase 1D 实施记录：model-after run-based wrapper

**日期**：2026-05-20（初版）；2026-05-21（workspace 共享 task、Ni/Pd 任务、路径约定）  
**状态**：已完成

## 目标

在 Phase 1C 的 `model_run_root` 上：加载 RankModel、读 task、positives mean 打分、任务适配 metrics、多 run 对比。产物只写显式 `output_dir`；**task 目录仅输入**。

**不做**：SQLite、Streamlit、重训、改 RankModel/loss、整包 legacy evaluation、MolCLR/ChemBERT、写 `ChemDB/training/evaluation/results` 与 `result_analysis`。

---

## 三层职责（避免混淆）

| 层级 | 位置 | 内容 |
|------|------|------|
| 任务 | `workspace/model_after_tasks/{task_id}/` | task.json + CSV；**workspace 级**，跨 project |
| 模型 | `workspace/projects/{pid}/runs/{run_id}/` | best.pt、index、npz、metal zscore |
| 结果 | `workspace/model_after_results/{task_id}/{batch_id}/` | ranking、metrics、多模型 summary |

CLI 用显式参数连接：`--task-dir`、`--model-run-root`、`--output-dir`。评估的仍是 **run 级模型**；结果路径按 **task + 批次** 归档，便于同一 task 对比多个 run。

---

## Workspace 约定

```text
workspace/
├── projects/.../runs/{run_id}/          # 模型
├── model_after_tasks/{task_id}/         # 任务（仅输入）
└── model_after_results/{task_id}/{batch_id}/   # 产物
```

- 核心库 **无** `projects/.../model_after_tasks` 硬编码；`load_task_bundle` 接受任意 `--task-dir`。
- CLI 默认输出：`model_after_results/{task_id}/{project}_{run_id}/`（evaluate-model 加 `_eval`）；多模型默认 `batch_default`（`--batch-id` 可覆盖）。
- `model_runs.csv` 可含不同 project 的 `model_run_root`。
- `model_after_results/**` 已 gitignore；task 目录已提交内容可进 git。

---

## 代码结构（`ChemDB/training/model_after/`）

| 模块 | 职责 |
|------|------|
| `bundle.py` | 从 run 加载 RankModel、npz、metal（仅 ecfp，路径须在 run 内） |
| `io.py` | 读 task.json + CSV |
| `scoring.py` | ECFP、positives mean context、打分 → ranking.csv |
| `eval_task.py` | held-out MRR/Hit@k、margin → metrics.json |
| `batch.py` | `model_runs.csv` 多 run、`best_model.json` |
| `README.md` | **科学逻辑改哪里**（非数据流） |

CLI：`app/services/model_after_manager.py` + `app/services/model_after_paths.py`。

**改评价科学**：`scoring.py`、`eval_task.py`、选模规则 `batch._selection_score`；I/O/加载见 `io.py` / `bundle.py`。

---

## 主要文件清单

**包**：`ChemDB/training/model_after/{bundle,io,scoring,eval_task,batch,__init__,README}.py`  
**服务**：`app/services/model_after_manager.py`、`model_after_paths.py`  
**模板**：`app/templates/model_after_task_template.json`、`model_runs_template.csv`  
**测试**：`tests/run_isolation/test_phase1d_model_after_run_based.py`（15 项默认 + 3 integration）  
**任务数据**：`workspace/model_after_tasks/ni_lkb_p_ni/`、`pd_lkb_p_cluster/`（自 legacy quest 转换）

---

## Legacy Ni/Pd → workspace 任务

源：`ChemDB_restructured/.../evaluation/quests/{ni_lkb_p_ni,pd_lkb_p_cluster}.csv`  
规则：`No`+`SMILES` → candidates；`Label=positive/negative` → 对应表；其余 label 仅在 candidates。

| task_id | 金属 | 候选 | 正例 | 负例 |
|---------|------|------|------|------|
| `ni_lkb_p_ni` | Ni | 344 | 11 | 0 |
| `pd_lkb_p_cluster` | Pd | 343 | 10 | 1 |

---

## 实测：`p001/runs/run_demo`（CPU）

| 任务 | output_dir | status | mrr | hit@50 | pos_neg_margin |
|------|------------|--------|-----|--------|----------------|
| Ni | `.../ni_lkb_p_ni/p001_run_demo` | success | 0.016 | 0.26 | null（无负例） |
| Pd | `.../pd_lkb_p_cluster/p001_run_demo` | success | 0.016 | 0.30 | 0.63 |

产物含 `ranking.csv`、`metrics.json`、`heldout_ranking.csv` 等；**未写入** task 目录。数值为跨任务适配，不与 legacy 全局 `training/models` 结果直接对比。

---

## CLI 示例

```bash
python -m app.services.model_after_manager evaluate-model \
  --model-run-root workspace/projects/p001/runs/run_demo \
  --task-dir workspace/model_after_tasks/ni_lkb_p_ni \
  --output-dir workspace/model_after_results/ni_lkb_p_ni/p001_run_demo \
  --device cpu
```

`recommend` / `evaluate-models` 同理；`evaluate-models` 需 `model_runs.csv`（宜放在 task 外或 results 批次目录，勿污染 task 输入目录）。

---

## Schema（摘要）

- **task.json**：`task_id`、`metal`、`embedding=ecfp`、CSV 文件名  
- **CSV**：`candidate_id|positive_id|negative_id, smiles, did, name, source`；有 smiles 优先现场 ECFP  
- **metrics.json**：mrr、hit@5/10/20/50、margin、status（positives&lt;2 → partial）  
- **model_runs.csv**：`model_id, model_run_root, checkpoint_path, ...`

---

## 测试

```bash
./scripts/run_pytest.sh -q                                    # 34 passed, 12 deselected
./scripts/run_pytest.sh tests/run_isolation/test_phase1d_model_after_run_based.py -v
# integration: -m integration_model_after_rank|eval|models
```

要点：共享 task 跨 project、默认不写 task_dir、no-write 源码目录、无 MolCLR import。

---

## 已知限制

1. 仅 ecfp；不从全局 `ChemDB/training/index.json` fallback。  
2. recommend 用 **全量 positives 均值** context；held-out 用 **正例子集组合**（与 legacy 略异）。  
3. 无 NDCG/AUC。  
4. 库直调且未传 `output_dir` 时，`scoring` 仍默认 `run/reports/model_after/`（CLI 已指向 `model_after_results`）。  
5. `orchestrator pipeline=model_after` 未接；验收以 `model_after_manager` CLI 为准。

---

## Phase 2 草案（SQLite，未做）

`model_after_tasks`、`model_after_runs`、`model_after_metrics`、`model_selection_batches` 表字段见原计划；`task_dir_path` 指向 workspace 绝对路径即可。

---

## 关联

- [phase_1d_plan.md](./phase_1d_plan.md)
- [phase_1abc_overall_assessment.md](./phase_1abc_overall_assessment.md)
- 包内说明：`ChemDB/training/model_after/README.md`
