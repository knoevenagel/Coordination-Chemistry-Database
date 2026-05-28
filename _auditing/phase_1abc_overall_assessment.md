# Phase 1A–1C 整体评估报告

| 字段 | 值 |
|------|-----|
| 类型 | 审计 / 总结（**不修改代码**） |
| 日期 | 2026-05-21 |
| 依据 | 当前 `ChemDBWebVersion` 仓库源码、`run_demo` 实测产物、pytest 与 `_auditing/phase_*` 记录 |
| 关联 | [phase_1a_implementation.md](./phase_1a_implementation.md)、[phase_1b_implementation.md](./phase_1b_implementation.md)、[phase_1c_implementation.md](./phase_1c_implementation.md) |

---

## 一、当前仓库结构

```
ChemDBWebVersion/
├── ChemDB/
│   ├── data/
│   │   ├── pubchem/              # 全量 PubChem（源码种子，非 run 写入）
│   │   ├── metal_list.txt
│   │   ├── p_elements_list.txt
│   │   └── metal_embedding/
│   │       ├── element_features.csv      # init-run 只读种子
│   │       ├── element_features_zscore.csv  # 仓库内可能存在；run 应写 run_root
│   │       └── zscore_element_features.py
│   ├── src/
│   │   ├── step1.py … step12.py, step_13_complete.py   # 核心链（Phase 1A CLI）
│   │   ├── main.py, server.py, proxy.py, molclr_api.py, affinity_api.py
│   │   └── tools/build_L3_embedding_index.py            # Phase 1C ECFP
│   └── training/                 # Phase 1C 迁入（仅 4 个源码文件，无 evaluation/）
│       ├── __init__.py, data.py, train.py, model.py
├── app/
│   ├── services/
│   │   ├── run_context.py
│   │   ├── step_registry.py
│   │   ├── orchestrator.py
│   │   └── training_manager.py
│   └── templates/
│       ├── training_index_template.json
│       └── training_config_template.yaml
├── tests/run_isolation/
│   ├── conftest.py, integration_helpers.py
│   ├── fixtures/{metal_list,p_elements_list,minimal_pubchem}.csv
│   ├── test_run_context.py, test_step_registry.py, test_two_runs.py
│   ├── test_orchestrator.py, test_no_write_to_repo_tmp.py
│   ├── test_integration.py, test_integration_1pct.py
│   └── test_phase1c_training_run_based.py
├── integration_data/pubchem_1pct/   # 1% 子集（~79 CSV + MANIFEST.json，~18M）
├── scripts/
│   ├── generate_pubchem_1pct_subset.py
│   └── run_pytest.sh
├── workspace/projects/p001/runs/run_demo/   # 实测成功 run 示例
├── TESTING.md
└── _auditing/
```

**说明**：`ChemDB_restructured/ChemDB/training/evaluation/` 含完整评估/推荐脚本，**未迁入** `ChemDBWebVersion`；模型后流程审计需对照该目录（见第八节）。

---

## 二、Phase 1A–1C 已完成内容总结

### Phase 1A

| 项 | 内容 |
|----|------|
| 新增 `app/services/` | `run_context.py`、`step_registry.py`、`orchestrator.py`（`__init__.py`、`__main__.py`） |
| `RunContext` 字段 | `run_root`, `data_dir`, `tmp_dir`, `training_dir`, `log_dir`, `report_dir`, `manifest_dir`；`create_dirs()` 含 `pubchem`, `metal_embedding`, `L3_embedding`, `training/ckpts` |
| `StepSpec` | `step_id`, `script`, `command_argv`, `required_inputs`, `expected_outputs`, `optional_outputs`, `script_kind`, `runner` |
| 核心 `PIPELINE_ORDER` | `step1` → `step2` → `step4_5` → `step6_7` → `step8` → `step12` → `step13` |
| Orchestrator CLI | `init-run`, `run-step`, `run-pipeline`（`--pipeline core` \| `training`） |
| 路径 CLI 补丁的 step | `step1`（`--pubchem-dir`, `--metal-list`, `--tmp-dir`）；`step2/4_5/6_7/8/12/13`（`--input-dir`/`--output-dir`）；`step8` 另需 `--metal-list` |
| Mock 测试 | `test_orchestrator`（patch `_run_subprocess` 模拟 step2）；`test_no_write`（同上）；非 subprocess 的 registry/context 纯单元测试 |

编排子进程：`cwd=ChemDB`，`PYTHONPATH=ChemDB:ChemDB/src`，`CHEMDB_RUN_ROOT={run_root}`。

### Phase 1B

| 项 | 内容 |
|----|------|
| 1% 生成 | `python scripts/generate_pubchem_1pct_subset.py --source ChemDB/data/pubchem --output integration_data/pubchem_1pct` |
| 子集路径 | `integration_data/pubchem_1pct/`（`MANIFEST.json`：`sample_rate=0.01`, `stratified_per_file`, `seed=42`） |
| 规模 | 源 1,946,214 行 → 输出 **19,460** 行，**79** 个金属 CSV，目录约 **18M** |
| pytest 真实跑通 | Tier A：`init-run`+`step1`+`step2`；Tier B：`run-pipeline` through **`step4_5`**；Tier C（可选）：`step6_7`；Enhanced（可选）：`step8` |
| 实测耗时 | Tier A+B：**~21 min**（`3 passed` in 1269.72s，见 phase_1b_implementation.md） |
| pytest 输出位置 | `tmp_path` 下临时 run（如 `run_1pct_tier_b`），**不**写入 `workspace/.../run_demo` |
| 未纳入最低验收 | `step8`–`step13`、training 链；Enhanced marker 单独 |

### Phase 1C

| 步骤 | 实现方式 |
|------|----------|
| ECFP | `ChemDB/src/tools/build_L3_embedding_index.py --tmp-dir --out-dir --backends ecfp` → `{run_root}/data/L3_embedding/L3_embedding_ecfp.npz` |
| metal zscore | `data/metal_embedding/zscore_element_features.py --input/--output` → run 内 zscore CSV |
| `training/index.json` | `training_manager.prepare_training_index` + 模板 |
| `training/config.yaml` | `training_manager.prepare_training_config`；`CHEMDB_TRAINING_INTEGRATION=1` → 3 epoch、cuda、early_stop=0 |
| `training.data` | `python -m training.data --training-dir --index --seed` |
| `training.train` | `python -m training.train --config` → `{run_root}/training/ckpts/best.pt`, `history.json` |
| 模型输出 | **`{run_root}/training/ckpts/best.pt`**（必需），`last.pt`（可选） |

| 测试类型 | 项 |
|----------|-----|
| Mock | `prepare_*` 进程内步骤；`test_training_manifest_paths_under_run_mock` 整条 training 链 |
| 真实 subprocess | `test_zscore_cli_*`；`test_build_l3_embedding_ecfp_*`（需 RDKit）；integration markers（默认不跑） |
| 真实训练模型 | **`run_demo` 手动 E2E**（2026-05-21 成功）；pytest `integration_training_*` 默认 **deselected** |

---

## 三、当前可用命令清单

**说明**：预处理 pipeline 在代码中名为 **`core`**，不是 `preprocessing`。

| 用途 | 命令示例 |
|------|----------|
| init-run | `python -m app.services.orchestrator init-run --run-root {run_root} [--pubchem-source {pubchem_source}]` |
| 单步 | `python -m app.services.orchestrator run-step --run-root {run_root} --step step1 [--stream]` |
| 核心链 | `python -m app.services.orchestrator run-pipeline --run-root {run_root} --pipeline core --through step13 [--stream]` |
| 训练链 | `export CHEMDB_TRAINING_INTEGRATION=1` 后 `… run-pipeline --pipeline training --through training_train` |
| 生成 1% | `python scripts/generate_pubchem_1pct_subset.py --source ChemDB/data/pubchem --output integration_data/pubchem_1pct` |
| 默认 pytest | `./scripts/run_pytest.sh` 或 `cd ChemDBWebVersion && ./scripts/run_pytest.sh -q` |
| 1% 集成 pytest | `export CHEMDB_INTEGRATION_PUBCHEM_SOURCE=$(pwd)/integration_data/pubchem_1pct`；`./scripts/run_pytest.sh -s tests/run_isolation/test_integration_1pct.py -m "integration_1pct_tier_a or integration_1pct_tier_b" -v` |
| training.data | `cd ChemDB && python -m training.data --training-dir {run_root}/training --index {run_root}/training/index.json --seed 42` |
| training.train | `cd ChemDB && python -m training.train --config {run_root}/training/config.yaml` |

**环境约定**

| 项 | 值 |
|----|-----|
| 推荐环境 | conda **`scidata`**（`./scripts/run_pytest.sh` 强制） |
| orchestrator 子进程 `cwd` | `{workspace_root}/ChemDB` |
| `PYTHONPATH` | `{ChemDB}` + `{ChemDB}/src`（编排器设置） |
| `CHEMDB_REPO_ROOT` | 默认 `ChemDBWebVersion/ChemDB`；可被环境变量覆盖 |
| 写 `workspace/` | `init-run`、`run-step`、`run-pipeline`、手动训练命令 |
| 只读 | `init-run` 从 `ChemDB/data/*` 复制种子；step 脚本在显式 CLI 下不写默认 `./tmp` |

占位符：`{workspace_root}` = 仓库根；`{run_root}` = 某次 run 绝对路径；`{pubchem_source}` = 目录或单个 CSV。

---

## 四、Run 目录规范与 `run_demo` 实测产物

示例：`workspace/projects/p001/runs/run_demo`（1% PubChem + core@step13 + training@train，2026-05-21）。

```
{run_root}/
├── data/
│   ├── pubchem/*.csv          # 79 个 1% 文件（init-run 复制）
│   ├── metal_list.txt, p_elements_list.txt
│   ├── metal_embedding/
│   │   ├── element_features.csv           [mandatory seed]
│   │   └── element_features_zscore.csv  [mandatory for training]
│   └── L3_embedding/
│       ├── L3_embedding_ecfp.npz          [mandatory]
│       └── build.log                      [optional]
├── tmp/                                   # core 链 ~134M
│   ├── step13_kl_nl_samples.csv           [mandatory]
│   ├── m_l3_pairs.csv, l3_gac.json      [mandatory]
│   ├── repaired_ligand_data.csv, metal_l3_index.csv  [mandatory for ECFP]
│   └── neo4j_*.csv, step*_stats.json, …   [optional / 中间产物]
├── training/
│   ├── index.json, config.yaml            [mandatory]
│   ├── split_index.json                   [mandatory after training_data]
│   ├── train/val/test_records.pkl         [mandatory after training_data]
│   └── ckpts/
│       ├── best.pt                        [mandatory]
│       ├── history.json                   [mandatory]
│       └── last.pt                        [optional]
├── logs/{step_id|training_*}.log          [mandatory 每步应有]
├── manifests/*.json                       [mandatory 每步应有]
└── reports/                               [optional，当前空]
```

**耗时（run_demo，manifest）**：core 墙钟约 **16 min**（step1 ≈864s）；training 六步合计约 **16 s**。

---

## 五、`training/index.json` 与 `config.yaml` schema

### index.json（真实结构，路径为绝对路径且在 `run_root` 下）

```json
{
  "run_root": "{abs_run_root}",
  "samples_csv": "{run_root}/tmp/step13_kl_nl_samples.csv",
  "m_l3_pairs": "{run_root}/tmp/m_l3_pairs.csv",
  "l3_gac": "{run_root}/tmp/l3_gac.json",
  "m_l3_pairs_path": "{run_root}/tmp/m_l3_pairs.csv",
  "kl_nl_ul_index": { "path": ".../step13_kl_nl_samples.csv", "format": "csv" },
  "ligand_embedding": {
    "backend": "ecfp", "mode": null,
    "path": ".../L3_embedding_ecfp.npz", "dir": ".../L3_embedding",
    "variants": { "ecfp": { "path": "...", "file": "L3_embedding_ecfp.npz" } },
    "npz_keys": ["dids", "smiles", "embeddings"]
  },
  "metal_embedding": {
    "path": ".../element_features_zscore.csv",
    "format": "csv", "lookup_column": "element"
  }
}
```

| 字段 | `training.data` 是否读取 | 说明 |
|------|-------------------------|------|
| `kl_nl_ul_index.path` | **是** | 作为 `kl_nl_ul_csv` |
| `m_l3_pairs_path` / `m_l3_pairs` | **是** | 划分唯一 T |
| `ligand_embedding.variants.*` | **是** | 经 `load_index_paths` → `ligand_embedding_ecfp` |
| `metal_embedding.path` | **是** | metal 特征表 |
| `run_root`, `samples_csv`, `l3_gac` | **否**（当前） | 语义/未来 SQLite、Registry |
| `ligand_embedding.backend/mode` | **否**（当前） | 由 `config.yaml` 的 `data.embedding` 决定训练用哪条 variant |

路径：**全部为绝对路径**；`training_manager` 写入前校验 `resolve_under_run`。

### config.yaml（run_demo 冒烟实例）

| 段 | 关键字段 |
|----|----------|
| `data` | `training_dir`, `index_path`, `embedding: ecfp`, `batch_size`, `max_candidates`, `seed` |
| `model` | `d_m: null`, `d_l: null`（训练时从 Dataset 推断） |
| `train` | `epochs`, `lr`, `ckpt_dir`, `device: cuda`, `early_stop_patience`, `scheduler` |

| 模式 | 差异 |
|------|------|
| 正式（`integration=False`） | `epochs=50`, `batch_size=32`, `device=cuda`, `early_stop=8`, `scheduler=reduce_on_plateau` |
| 冒烟（`CHEMDB_TRAINING_INTEGRATION=1`） | `epochs=3`（可 `CHEMDB_TRAINING_INTEGRATION_EPOCHS` 覆盖）, `batch_size=4`, `early_stop=0`, `scheduler=null` |

`d_m`/`d_l`：`train.py` 若 config 为 null，从 `RankDataset` 加载的 embedding 维度推断。

---

## 六、测试结果汇总

**收集命令**（scidata，2026-05-21）：`./scripts/run_pytest.sh --collect-only -q` → **19 collected / 9 deselected**（共 28 项）。

**默认**：`./scripts/run_pytest.sh -q` → **`19 passed, 9 deselected` in ~1.0s**。

| 类别 | 文件 | 性质 |
|------|------|------|
| Phase 1A | `test_run_context`, `test_step_registry`, `test_two_runs`, `test_orchestrator`, `test_no_write` | 多数单元；step2 **mock** `_run_subprocess` |
| Phase 1B | `test_integration_1pct`（3+2 用例） | **真实** subprocess + RDKit；**deselected** 默认 |
| Phase 1C | `test_phase1c_training_run_based`（~12 默认） | 混合：zscore/ECFP **真实**；manifest 链 **mock** |

| Integration marker | 默认 | 真实训练 |
|------------------|------|----------|
| `integration_1pct_tier_a/b` | deselected | 否 |
| `integration_1pct_tier_c`, `integration_1pct_enhanced` | deselected | 否 |
| `integration_embedding` | deselected | 否 |
| `integration_training_data/train` | deselected | 否（**run_demo 已真实训练**） |

**耗时（文档/实测）**：Tier A+B ~21 min；run_demo core ~16 min；run_demo training ~16 s。

**Flaky**：step1 末段 tqdm 长时间停在 99%（数据特性，非 pytest 死锁）；无自动 retry。

**no-write guard**

| 机制 | 范围 |
|------|------|
| `conftest` autouse | `ChemDB/tmp/` 文件名集合不变；`ChemDB/src/tmp` 不得存在 |
| `test_phase1c` `REPO_ARTIFACT_PATTERNS` | `ChemDB/data/L3_embedding/*.npz`；`element_features_zscore.csv`；`ChemDB/training/{index,config,split_index,*_records.pkl,ckpts/**}` |

允许 `ChemDB/training/{data,train,model,__init__}.py` 源码变更。

---

## 七、no-write 保障与残留风险

### 已保障（测试 + 编排设计）

- `ChemDB/tmp/`、`ChemDB/src/tmp/`
- `ChemDB/data/L3_embedding/` 运行 npz
- `ChemDB/data/metal_embedding/element_features_zscore.csv`（运行生成）
- `ChemDB/training/` 下运行产物：`index.json`, `config.yaml`, `*.pkl`, `split_index.json`, `ckpts/*`

### 仍可能污染源码目录（**ChemDBWebVersion 内**）

| 路径 | 风险类型 |
|------|----------|
| `ChemDB/src/step*.py` 默认 `OUTPUT_DIR=./tmp` | 无 CLI 直接 `python step1.py` 写 `ChemDB/tmp` |
| `ChemDB/src/tools/build_L3_embedding_index.py` | 默认 `REPO_ROOT/data/L3_embedding` |
| `ChemDB/training/data.py` / `train.py` | 无 `--training-dir` 时默认 `Path(__file__).parent` |
| `ChemDB/src/main.py` | 旧编排，硬编码相对路径 |
| `ChemDB/src/server.py`, `proxy.py`, `molclr_api.py`, `affinity_api.py` | 服务/API，未纳入 run-based 测试 |
| `ChemDB/src/step_13_complete.py` | 注释 fallback `result_analysis`（CLI 已可覆盖） |
| `ChemDB_restructured/.../training/evaluation/*` | 默认 `training/index.json`、`evaluation/results/`（**不在 WebVersion 包内**，但可能被手动执行） |

---

## 八、模型后流程现状审计

**结论**：`ChemDBWebVersion/ChemDB/training/` **无** `evaluation/` 包；模型后能力在 **`ChemDB_restructured/ChemDB/training/evaluation/`**，尚未 run-based 化。

| file_path | 用途 | 输入 | 输出 | 全局路径 | 可 run-based | 最小改造 |
|-----------|------|------|------|----------|-------------|----------|
| `ChemDB_restructured/.../evaluation/core.py` | 单 ckpt 对 quest CSV 评估 | quest CSV（No,SMILES,Label）；`index.json`；`*.pt` | `*_context.csv`, `*_heldout.csv`, `*_summary_by_heldout.csv` | 默认 `training/index.json` | 需 `--index`/`--output-dir` | 迁入 WebVersion + CLI 显式路径 |
| `ChemDB_restructured/.../evaluation/run_all_quests.py` | 批量模型×quest | `models/`, `quests/`, `index` | `evaluation/results/` | **是** | 否 | 改为 `{run_root}` 下 models 与 results |
| `ChemDB_restructured/.../evaluation/quests/*.csv` | 任务定义 | No, SMILES, **Label**（positive/其他） | — | — | — | 与训练 records 格式不同 |
| `ChemDBWebVersion/ChemDB/training/train.py` | 训练 + val **Recall@1/5, MRR** | records pkl + npz + metal csv | `ckpts/best.pt`, `history.json` | fallback `training/` | 已支持 `--config` 绝对路径 | 无推理 CLI |
| `ChemDBWebVersion/ChemDB/training/data.py` | 划分与 Dataset | index.json 路径字段 | `*_records.pkl` | fallback 包目录 | 已 `--training-dir` | — |
| `ChemDBWebVersion/tmp/step13_kl_nl_samples.csv` | KL/NL 样本 | T,M,label,candidate_did,… | — | — | 已在 run | 非 quest 格式 |

### 关键问答（基于现状）

1. **给定 checkpoint + 候选列表 → 排名？**  
   **WebVersion 无独立脚本**。`evaluation/core.py` 的 `score_and_rank` 可复用，但输入是 **整数 id + SMILES 表**，需自写 wrapper；非 `candidate_did` 列表 CLI。

2. **正例/负例/候选 → 模型适配评估？**  
   **有（restructured）**：quest CSV 中 `Label=positive` 为已知正例，其余 valid 为候选；held-out 评估 **MRR、hit@5/10/20/50**。**无**显式 negative 列；NL 在训练 records 中，不在 quest 评估流。

3. **批量评估多模型？**  
   **有**：`run_all_quests.py`（`models_dir` × `quests_dir`）。

4. **Quest CSV 格式**：`No, SMILES, Label`（`pd_lkb_p_cluster.csv` 等列很多，评估只用前三列逻辑）。

5. **训练输入格式**：`train_records.pkl` 内 `T,M,kl_dids,nl_dids,ul_dids`（L3 DID）；与 quest 的整数 id **不一致**。

6. **加载模型需要**：`best.pt`（含 `model_state_dict` + 嵌套 `config`）；评估时另需 **`index.json` 的 metal_embedding.path**；配体向量对 ECFP **现场算 SMILES** 或需 L3 npz（训练用 DID 查表）。

7. **除 best.pt 外**：**需要** config（ckpt 内嵌 + 文件）、index（metal）、embedding 源；**不需要** split_index/pkl 做推理。

8. **输出**：context/heldout/summary CSV；非 `ranking.csv`。

9. **指标**：**MRR、hit@k**（evaluation/core）；训练 val 为 **Recall@1/5、MRR**；legacy 脚本有 **NDCG@10/20**；**未见 weaklift/AUC** 于主路径。

10. **仍写全局目录**：`run_all_quests.py` → `training/evaluation/results/`；`generate_sweep_configs` → `training/results_*`；**不在** ChemDBWebVersion 默认包内。

---

## 九、Model Registry 准备度

基于 `run_demo` 成功训练，建议登记字段：

| 字段 | 来源示例 |
|------|----------|
| `model_id` | `{project_id}_{run_id}_ecfp_smoke` 或 uuid |
| `project_id` / `run_id` | `p001` / `run_demo` |
| `run_root` | manifest `run_root` |
| `checkpoint_path` | `training/ckpts/best.pt` |
| `config_path` | `training/config.yaml` |
| `index_path` | `training/index.json` |
| `embedding_backend` | `ecfp`（index + config） |
| `ligand_embedding_path` | `data/L3_embedding/L3_embedding_ecfp.npz` |
| `metal_embedding_path` | `data/metal_embedding/element_features_zscore.csv` |
| `samples_csv` | `tmp/step13_kl_nl_samples.csv` |
| `m_l3_pairs_path`, `l3_gac_path` | tmp 内路径 |
| `split_index_path`, `*_records.pkl` | training/ |
| `history_path` | `training/ckpts/history.json` |
| `metrics` | history.json 内 recall/MRR/loss |
| `train.device`, `epochs`, `seed` | config.yaml |
| `created_at` | manifest `finished_at` |
| `git_commit` | **当前未写入 manifest** |
| `step_manifests` | `manifests/training_train.json` 等 |

**Manifest 缺口**：无统一 `artifact` 类型字段；无 `model_id`；无 git/version；`pipeline_*.json` 无逐步 artifact 哈希；prepare 步骤 manifest 较简。够做 **文件系统 Registry**；**不够**直接映射 SQLite `Model` 表 without 扩展 schema。

---

## 十、下一阶段建议

### 是否先做 Phase 1D 再做 SQLite？

**同意。** 理由：模型后脚本仍在 **restructured + 全局路径**；无 run-based 推理/评估则 SQLite 只能存路径却无法端到端复现。Phase 1C 已证明 **per-run 训练产物闭环**；1D 应复用同一 `run_root` 契约。

### Phase 1D 范围（建议）

- 不做 SQLite / Streamlit。
- 从 `ChemDB_restructured/.../evaluation/` **最小迁入** `core.py` + 薄 CLI（或 `app/services/recommendation_manager.py`）。
- 输入契约（需新设计，与 step13/quest 对齐）：
  - `candidates.csv`：ligand DID 或 SMILES；
  - `positives.csv` / `negatives.csv`（可选）；
  - `metal`，`checkpoint`，`index.json`，`run_root`。
- 输出：`ranking.csv`、`metrics.json`（可选 `report.json`）。
- 模型加载接口：`load_run_model(run_root, ckpt=best.pt)` → RankModel + paths。

**优先改动文件（建议顺序）**

1. 迁入/封装 `evaluation/core.py` → `ChemDB/training/inference.py`（或 `app/services/`）
2. 新 CLI `python -m training.recommend` 或 orchestrator step `recommend` / `evaluate_models`
3. `index.json` 扩展可选 `inference` 段（仍绝对路径）
4. `tests/run_isolation/test_phase1d_*.py`：mock + run_demo `best.pt` 冒烟
5. `_auditing/phase_1d_plan.md`

**Phase 1D 最低测试**

- 给定 `run_demo` 的 `best.pt` + mini candidate CSV → `ranking.csv` 存在且行数>0；
- 输出仅在 `{run_root}/reports/` 或 `{run_root}/training/eval/`；
- no-write guard 扩展 eval 输出目录；
- 不调用 `run_all_quests` 全局 `results/`。

### Phase 2 SQLite（1D 后）

表建议：`Project`, `Run`, `StepExecution`, `Artifact`, `Model`, `RecommendTask`, `EvaluationResult`。大文件只存路径 + `manifests/*.json` 摘要。字段可由第九节映射。

### Phase 3 Streamlit

依赖 Phase 2 元数据 + Phase 1D 可触发评估/推荐的 API 或 CLI 封装。

---

## 附录：与「模型后六步」的映射

| 目标能力 | 当前状态 |
|----------|----------|
| 1. 模型库管理 | 仅文件系统（`run_root/training/ckpts`）；无 Registry API |
| 2. 用户输入 pos/neg/candidate | 无 Web UI；quest CSV（restructured）仅部分覆盖 |
| 3. 自动评估多模型吻合度 | `run_all_quests.py`（全局路径，未迁入） |
| 4. 选择最合适模型 | 人工看 summary CSV |
| 5. 推荐结果 | `score_and_rank` 存在，无 run-based CLI |
| 6. SQLite + Streamlit | 未开始 |

**Phase 1A–1C 已完成**：run_root 下 **数据准备 → 训练 → `best.pt`**。**未完成**：run-based **推理 / 评估 / 推荐 / 模型库元数据**。
