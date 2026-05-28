# Phase 1D 实施计划：Model-After Run-Based Wrapper

| 字段 | 值 |
|------|-----|
| 状态 | **待确认**（本文件仅为计划，不修改代码） |
| 前置 | [phase_1abc_overall_assessment.md](./phase_1abc_overall_assessment.md)、Phase 1A–1C 已完成 |
| 目标 | 在 `{run_root}` 内完成模型加载、候选打分、任务适配评估、同 run 多 ckpt 批量对比；标准化 I/O；manifest/log；no-write 测试 |
| 不在范围 | SQLite、Streamlit、UI、DB 表、重训、改 RankModel/损失、MolCLR/ChemBERT、迁入整个 `evaluation/legacy`、写 `ChemDB/training/evaluation/results` |

---

## 0. 现状依据（代码事实）

| 事实 | 来源 |
|------|------|
| WebVersion `ChemDB/training/` 仅有 `data.py`, `train.py`, `model.py`, `__init__.py` | 仓库目录 |
| 可复用逻辑在 `ChemDB_restructured/ChemDB/training/evaluation/core.py` | `load_task_data`, `score_and_rank`, `run_eval_for_context`, `evaluate_checkpoint_for_task`, `load_rank_model` |
| 批量入口 `run_all_quests.py` 写 `evaluation/results/`，默认 `training/index.json`、`training/models/*.pt` | restructured 脚本 |
| Phase 1C 产物契约 | `{run_root}/training/index.json`, `config.yaml`, `ckpts/best.pt`, `data/L3_embedding/L3_embedding_ecfp.npz`, `data/metal_embedding/element_features_zscore.csv` |
| 训练时 context = **KL DID 向量 mean**；candidates = T+NL+UL（`RankDataset`，161 维） | `ChemDB/training/data.py` |
| 评估时 context = **正例整数 id 组合** mean embedding；候选 = quest 表内全部 valid id | `evaluation/core.py` |
| `RunContext` 已有 `report_dir` | `app/services/run_context.py` |
| Orchestrator 已支持 `pipeline=core|training`、in-process / subprocess、`manifest` 增强结构 | `app/services/orchestrator.py` |

**Phase 1D 只支持 `embedding=ecfp`**（与 Phase 1C 一致）；`get_embeddings_for_ids` 中 molclr 分支**不迁入**。

---

## 1. 目标 run 目录扩展

在 Phase 1C 的 `{run_root}` 上增加（不改动既有 `training/` 语义）：

```text
{run_root}/
├── reports/                          # Phase 1D 主输出根（新建，create_dirs 可选补齐）
│   ├── ranking.csv                   # 单模型候选排名
│   ├── metrics.json                  # 任务适配指标摘要
│   ├── report.json                   # 运行摘要（模型、输入、输出路径）
│   └── eval/                         # 详细评估（可选深度输出）
│       └── {checkpoint_stem}/
│           ├── {stem}_context.csv
│           ├── {stem}_heldout.csv
│           └── {stem}_summary_by_heldout.csv
├── reports/model_batch/              # 多模型批量对比（可选）
│   ├── summary_by_model.csv
│   └── run_manifest.json
├── input/                            # 可选：用户放入本次任务输入
│   ├── task.json                     # 推荐：统一任务描述
│   ├── candidates.csv
│   ├── positives.csv
│   └── negatives.csv                 # 可选
├── logs/
│   └── {step_id}.log
└── manifests/
    └── {step_id}.json
```

**原则**：所有 Phase 1D 写入路径必须在 `run_root` 下；禁止默认写 `ChemDB/training/evaluation/results` 或 `result_analysis`。

---

## 2. 输入 / 输出格式（标准化）

### 2.1 推荐：`task.json`（主契约）

放置：`{run_root}/input/task.json` 或通过 CLI `--task {path}` 传入（须在 run 内或复制进 run）。

```json
{
  "task_id": "demo_rank_001",
  "metal": "Pd",
  "embedding": "ecfp",
  "checkpoint": "training/ckpts/best.pt",
  "positives": [
    {"id": "1", "smiles": "CCO"},
    {"id": "2", "smiles": "CC(C)O"}
  ],
  "negatives": [
    {"id": "99", "smiles": "c1ccccc1"}
  ],
  "candidates": [
    {"id": "1", "smiles": "CCO"},
    {"id": "2", "smiles": "CC(C)O"},
    {"id": "3", "smiles": "CCC"}
  ],
  "eval": {
    "max_k": null,
    "mode": "heldout"
  }
}
```

| 字段 | 用途 |
|------|------|
| `metal` | 查 `index.json` → `metal_embedding.path` |
| `embedding` | 必须为 `ecfp`；与 `index` / ckpt 内 config 一致 |
| `checkpoint` | 相对 `run_root` 或绝对路径（须在 run 内） |
| `positives` | **打分**：KL/context 池化（见 §4.1）；**评估**：held-out 正例集 |
| `negatives` | Phase 1D **仅记录/透传**，不参与 `evaluate_checkpoint_for_task`（与旧 quest 逻辑一致）；供后续推荐过滤扩展 |
| `candidates` | 打分排序列表 |

### 2.2 CSV 兼容（quest 与简化）

| 文件 | 列 | 兼容 |
|------|-----|------|
| `candidates.csv` | `id`, `smiles`（或 `No`, `SMILES`） | 打分 |
| `positives.csv` | `id`, `smiles`, `label`（positive） | 评估 + 可选作 context 来源 |
| `negatives.csv` | 同上，`label=negative` | 元数据；v1 不进入 held-out 公式 |

**Quest 兼容评估**：若提供单个 `quest.csv`（列 `No`, `SMILES`, `Label`），可映射为 `evaluate_checkpoint_for_task` 所需格式（与 restructured `load_task_data` 一致）。金属符号由 `task.json` 的 `metal` 或文件名 `pd_*` / `ni_*` 推断规则**仅作兼容**，主路径用 `task.json`。

### 2.3 输出

| 文件 | 内容 |
|------|------|
| `reports/ranking.csv` | `id`, `smiles`, `rank`, `score`（按 score 降序） |
| `reports/metrics.json` | `mrr`, `hit_at_5`, `hit_at_10`, `hit_at_20`, `hit_at_50`, `n_positives`, `n_candidates`, `model`, `metal` |
| `reports/report.json` | 输入路径、输出路径、checkpoint、embedding、耗时、status |
| `reports/eval/{stem}_*.csv` | 移植自 `core.py` 的详细表（评估步骤） |
| `reports/model_batch/summary_by_model.csv` | 多 ckpt 对比（`mrr_mean`, `hit_at_*_mean`） |

---

## 3. 代码结构（新增，最小迁入）

### 3.1 新增文件

| 路径 | 作用 |
|------|------|
| `ChemDB/training/model_after/__init__.py` | 包入口 |
| `ChemDB/training/model_after/bundle.py` | **模型加载 wrapper**：`load_run_bundle(run_root, ckpt)` |
| `ChemDB/training/model_after/scoring.py` | **单模型打分**：`score_and_rank`（从 core 抽取，ECFP-only） |
| `ChemDB/training/model_after/eval_task.py` | **任务适配评估**：`evaluate_task_fit`（从 `evaluate_checkpoint_for_task` 抽取） |
| `ChemDB/training/model_after/batch.py` | **同 run 多 ckpt 批量**（替代 `run_all_quests` 的全局目录语义） |
| `ChemDB/training/model_after/io.py` | 读写 `task.json` / CSV → 内部结构 |
| `app/services/model_after_manager.py` | 编排层调用：路径校验、写 reports、返回 manifest payload |
| `app/templates/model_after_task_template.json` | `task.json` 模板 |

**不新增**：`evaluation/legacy/`、`run_all_quests.py` 整文件复制、`results_*` 脚本。

### 3.2 `load_run_bundle`（核心 wrapper）

```python
# 逻辑签名（计划）
load_run_bundle(
    run_root: Path,
    checkpoint: Path | None = None,  # 默认 training/ckpts/best.pt
    device: str | None = None,       # 默认 cuda，复用 resolve_training_device
) -> RunModelBundle
```

`RunModelBundle` 字段（dataclass）：

| 字段 | 来源 |
|------|------|
| `run_root` | 参数 |
| `index` | `training/index.json` |
| `config` | `training/config.yaml` 或 ckpt 内嵌 |
| `ckpt_path` | 解析后绝对路径 |
| `model` | `RankModel` + `load_state_dict`（复用 `core.load_rank_model` 逻辑） |
| `device` | `torch.device` |
| `ligand_lookup`, `d_l` | `training.data.load_ligand_embedding(index ecfp npz)` |
| `metal_lookup`, `d_m` | `training.data.load_metal_embedding(index metal path)` |
| `embedding_backend` | 固定 `"ecfp"` |

**校验**：

- `index.json` 路径均在 `run_root` 下；
- `ligand_embedding.backend == "ecfp"`；
- npz / metal csv / ckpt 存在；
- `d_m`, `d_l` 与 ckpt 中 `model_state_dict` 兼容（维度不一致则失败）。

### 3.3 打分逻辑（对齐训练语义）

**训练**（`RankDataset`）：`context = mean(KL vectors)`，`candidates = stack(T, NL, UL)`。

**Phase 1D 推理**（新 API `rank_candidates`，避免 161 固定候选数）：

| 输入 | 映射 |
|------|------|
| `positives[]` | 当作 **KL 池化集合**（与 eval 的 context 一致：对 id/smiles 取 embedding 后 mean） |
| `candidates[]` | 候选向量矩阵，逐条打分 |
| `metal` | metal 向量 |

实现：复用 `core.score_and_rank(model, metal_vec, context_vec, candidate_embeddings, candidate_ids, device)`。

- 候选 id 可为整数或字符串；**smiles** 用 RDKit ECFP（与 `core.compute_ecfp_embeddings` 一致，2048 bit）。
- 若 id 可映射到 `ligand_lookup`（归一化 DID），**优先用 npz**（与训练一致），否则 fallback ECFP(smiles)。

输出 `reports/ranking.csv`。

### 3.4 任务适配评估（移植 `evaluate_checkpoint_for_task`）

| 项 | 说明 |
|----|------|
| 输入 | `quest.csv` 或 `positives.csv`+`candidates.csv` 合成 quest 视图 |
| 逻辑 | 正例 `label=positive`；held-out 组合评估；指标 **MRR, hit@5/10/20/50** |
| 输出 | `reports/metrics.json` + `reports/eval/{stem}_*.csv` |
| 限制 | 仍需 **≥2 个 positive**（与 core 一致） |
| ECFP | 仅 `model_type=ecfp`；禁止走 molclr 分支 |

`negatives`：写入 `report.json` 统计，**不改变** held-out 公式（与现网评估一致）。

### 3.5 多模型批量（最小 wrapper）

**范围**：仅同一 `run_root` 内多个 `.pt`（如 `best.pt`, `last.pt`，或 `reports/model_batch/checkpoints/*.pt` 列表文件）。

| 项 | 说明 |
|----|------|
| 输入 | `task.json` + `--checkpoints rel1,rel2` 或目录扫描 |
| 行为 | 对每个 ckpt 调 `evaluate_task_fit`（或仅汇总 ranking 指标若任务为 rank 模式） |
| 输出 | `reports/model_batch/summary_by_model.csv` |
| 不做 | 跨 run、跨 index 的模型库（留给 Phase 2 SQLite） |

---

## 4. Orchestrator 扩展

### 4.1 新 pipeline：`model_after`

```text
MODEL_AFTER_PIPELINE_ORDER =
  validate_model_after_inputs   # in-process
  recommend_rank                # in-process（或 subprocess，优先 in-process 便于调试）
  evaluate_task_fit             # in-process
  batch_compare_models          # in-process，optional：仅当 task 指定多 ckpt
```

**简化默认三步**（batch 合并进 `evaluate_task_fit` 的可选参数）：

| step_id | runner | 说明 |
|---------|--------|------|
| `validate_model_after_inputs` | in-process | 检查 index/ckpt/npz/task 输入 |
| `recommend_rank` | in-process | 写 `ranking.csv` |
| `evaluate_task_fit` | in-process | 写 `metrics.json` + eval csv |

`batch_compare_models`：当 `task.json` 内 `checkpoints: [...]` 长度 >1 时由 orchestrator 追加，或单独 `--through batch_compare_models`。

### 4.2 CLI 示例

```bash
# 全流程（需已有 Phase 1C 产物 + input/task.json）
python -m app.services.orchestrator run-pipeline \
  --run-root {run_root} \
  --pipeline model_after \
  --through evaluate_task_fit \
  --stream
```

```bash
# 单步
python -m app.services.orchestrator run-step \
  --run-root {run_root} \
  --step recommend_rank
```

### 4.3 独立模块 CLI（便于调试，与 training.data 并列）

```bash
cd {ChemDB_repo_root}
python -m training.model_after.rank \
  --run-root {run_root} \
  --task {run_root}/input/task.json \
  --output {run_root}/reports/ranking.csv

python -m training.model_after.eval_task \
  --run-root {run_root} \
  --task {run_root}/input/task.json \
  --output-dir {run_root}/reports/eval
```

`cwd`、`PYTHONPATH` 与 Phase 1C 相同：`cwd=ChemDB`，`PYTHONPATH=ChemDB:ChemDB/src`。

---

## 5. 依赖关系（DAG）

```text
Phase 1C 已完成
  training/index.json
  training/ckpts/best.pt
  data/L3_embedding/L3_embedding_ecfp.npz
  data/metal_embedding/element_features_zscore.csv
        │
        ▼
用户准备 input/task.json（或 CSV）
        │
        ├─► recommend_rank ──► reports/ranking.csv
        │
        └─► evaluate_task_fit ──► reports/metrics.json
                              └─► reports/eval/*.csv
```

`recommend_rank` 与 `evaluate_task_fit` 可独立运行；evaluate 不依赖 ranking 输出。

---

## 6. Manifest / Log

复用 Phase 1C 结构，每 logical step 写入：

- `{run_root}/logs/{step_id}.log`
- `{run_root}/manifests/{step_id}.json`

**`inputs[]` / `outputs[]`**：绝对路径、`exists`、`size_bytes`；`cwd=ChemDB`；`run_root` 必填。

**`recommend_rank` outputs**：`reports/ranking.csv`（mandatory）

**`evaluate_task_fit` outputs**：`reports/metrics.json`（mandatory），`reports/eval/...`（mandatory 至少 summary 或 metrics 内嵌）

---

## 7. 测试计划

新文件：`tests/run_isolation/test_phase1d_model_after_run_based.py`

| # | 测试 | 类型 |
|---|------|------|
| 1 | `load_run_bundle` 读取 `run_demo` 的 best.pt + index | 真实（需 run_demo 或 fixture run） |
| 2 | `rank_candidates` 小 task → ranking.csv 行数正确 | 真实 RDKit + 小 npz |
| 3 | `evaluate_task_fit` mini quest（≥2 positive）→ metrics.json 字段存在 | 真实或 skip 若无 GPU |
| 4 | 所有输出路径在 run_root 内 | 单元 |
| 5 | `test_no_write_repo_model_after` 快照 `ChemDB/training/evaluation`、`result_analysis` | 与 Phase 1C guard 扩展 |
| 6 | orchestrator `recommend_rank` manifest 结构 | mock in-process |
| 7 | 禁止 molclr import 于 model_after 包 | 静态检查 / import 测试 |
| 8 | 两 run 的 reports 互不覆盖 | 两 run_root |

**Integration marker**（可选，默认 deselected）：

- `integration_model_after_rank`：对 `run_demo` 或 1% run 做 rank
- `integration_model_after_eval`：对 mini quest 做 eval

**默认 pytest**：mock 为主；不强制 GPU。

**运行**：`./scripts/run_pytest.sh`（scidata）。

---

## 8. no-write 保障

### 8.1 测试 guard 扩展

在 Phase 1C 模式基础上增加（若目录存在则快照）：

- `ChemDB/training/evaluation/`（整个目录无新增文件）
- `ChemDB/result_analysis/`（若存在）
- 禁止 `ChemDB/training/models/`（restructured 惯例）

允许：`ChemDB/training/model_after/*.py` 源码新增。

### 8.2 代码层约束

- `model_after_manager` 内 `ctx.resolve_under_run()` 校验所有读写路径；
- 禁止默认 `results_dir=Path(__file__).parent / "results"`；
- `batch.py` 不得 import `run_all_quests` 的全局默认路径。

---

## 9. 与 Phase 2 / Streamlit 的接口

Phase 1D 结束后，SQLite / Streamlit 只需调用：

| 能力 | 调用方式 |
|------|----------|
| 加载模型 | `load_run_bundle(run_root, model_id → ckpt 路径)` |
| 推荐 | orchestrator `recommend_rank` 或 Python API `rank_candidates` |
| 评估 | `evaluate_task_fit` |
| 元数据 | 读 `manifests/*.json` + `reports/report.json` |

**Model 表字段** 可从 `training/index.json` + `manifests/training_train.json` + `reports/metrics.json` 汇总（见 overall assessment §9）。

---

## 10. 实施顺序（确认后执行）

1. 新增 `ChemDB/training/model_after/` 包（从 `evaluation/core.py` **抽取** ECFP 函数，删除 molclr 路径）。
2. 实现 `bundle.py` + `scoring.py`（`rank_candidates`）。
3. 实现 `eval_task.py`（`evaluate_task_fit`）。
4. 实现 `batch.py`（同 run 多 ckpt）。
5. 新增 `model_after_manager.py` + `task` 模板。
6. 扩展 `step_registry` + `orchestrator`（`pipeline=model_after`）。
7. `run_context.create_dirs` 可选增加 `reports/`, `input/`。
8. 新增 pytest + pytest markers。
9. 使用 `run_demo` 做手动冒烟（rank + eval）。
10. 撰写 `_auditing/phase_1d_implementation.md`。

---

## 11. 风险与明确不做

| 风险 | 缓解 |
|------|------|
| 训练 DID（L3_*）与 quest 整数 id 不一致 | 文档标明；rank API 同时支持 smiles；eval 沿用 quest 整数 id |
| 推理 context 用 positives 而非训练 KL 列表 | 与 `evaluate` 一致；在 `report.json` 注明 `context_source=positives_mean` |
| 无 GPU 环境失败 | 与 training 一致，显式报错 |
| 负例未参与 loss/评估 | Phase 1D 文档化；后续 Phase 可扩展 contrastive |

| 不做 | 原因 |
|------|------|
| 跨 run 模型库 | Phase 2 |
| NDCG/AUC/weaklift | legacy 部分有；非 1D 最小集 |
| 161 固定候选推理 | 仅训练 DataLoader 需要；推理按实际候选数动态 K |
| Streamlit 表单 | Phase 3 |

---

## 12. 验收标准（Phase 1D Done）

1. 给定含 Phase 1C 产物的 `{run_root}` + `input/task.json`，`recommend_rank` 生成 `reports/ranking.csv`。
2. `evaluate_task_fit` 生成 `reports/metrics.json`（含 MRR、hit@k）。
3. 同 run 内 2 个 ckpt 批量对比生成 `summary_by_model.csv`。
4. 默认 pytest 通过，无 repo evaluation/results 污染。
5. `run_demo` 手动跑通 rank 或 eval 至少各一项并记录于 `phase_1d_implementation.md`。

---

## 13. 待用户确认项

1. 输出根目录用 `reports/` 还是 `training/eval/`？（计划默认 **`reports/`**，与 `RunContext.report_dir` 一致。）
2. `recommend_rank` 的 context 是否必须来自 `positives` 列表？（计划：**是**，与评估 context 一致。）
3. 是否必须同时交付 `batch_compare_models` 为独立 step？（计划：**可选**，多 ckpt 时启用。）
4. 输入是否强制 `task.json`？（计划：**推荐**；CSV-only 作为兼容 reader。）

确认后状态改为「可实施」，并另写 `phase_1d_implementation.md`（完成版）。
