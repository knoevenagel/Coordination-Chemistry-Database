# Phase 1C 实施计划（Embedding + Training Run-based）

| 字段 | 内容 |
|------|------|
| 阶段 | Phase 1C |
| 状态 | **待用户确认**（本文档为实施计划，确认后再改代码） |
| 完成日期 | — |
| 前置 | [phase_1a_implementation.md](./phase_1a_implementation.md)、[phase_1b_implementation.md](./phase_1b_implementation.md) |
| 范围 | ECFP-only L3 embedding、metal zscore、per-run `training/index.json` + `config.yaml`、`training.data`、`training.train` |
| 不在范围 | Streamlit、SQLite、API、evaluation、recommendation、ChemBERT/MolCLR/GIN/GCN、RankModel/loss/sampling 算法改动 |

确认实施后，落地记录写入 [phase_1c_implementation.md](./phase_1c_implementation.md)（实施完成版）。

---

## 1. Run 目录目标结构

每个 `run_root` 在 Phase 1C 完成后应满足（与 Phase 1A 目录约定一致，并补齐 embedding/training 子树）：

```text
{run_root}/
├── data/
│   ├── pubchem/                          # init-run 复制（Phase 1A/1B）
│   ├── metal_list.txt
│   ├── p_elements_list.txt
│   ├── metal_embedding/
│   │   ├── element_features.csv          # init-run 从 ChemDB 种子复制
│   │   └── element_features_zscore.csv   # zscore_metal_embedding
│   └── L3_embedding/
│       └── L3_embedding_ecfp.npz         # build_l3_embedding_ecfp（Phase 1C 仅 ECFP）
├── tmp/
│   ├── repaired_ligand_data.csv          # step2
│   ├── metal_l3_index.csv                # step13（build_l3 依赖）
│   ├── step13_kl_nl_samples.csv          # step13
│   ├── m_l3_pairs.csv                    # step12
│   ├── l3_gac.json                       # step12
│   └── ...                               # 核心链其它产物
├── training/
│   ├── index.json                        # prepare_training_index
│   ├── config.yaml                       # prepare_training_config
│   ├── split_index.json                  # training_data
│   ├── train_records.pkl
│   ├── val_records.pkl
│   ├── test_records.pkl
│   └── ckpts/
│       ├── best.pt                       # training_train
│       ├── last.pt                       # optional
│       └── history.json
├── logs/
│   └── {step_id}.log                     # 每 logical step
├── reports/
└── manifests/
    └── {step_id}.json                    # Phase 1C 增强 manifest
```

### 1.1 禁止写入的仓库路径（硬约束）

测试与编排层均需保证**不产生新的运行产物**于：

| 禁止路径 | 说明 |
|----------|------|
| `ChemDB/data/L3_embedding/` | 全局 L3 npz |
| `ChemDB/data/metal_embedding/`（运行生成物） | 仅允许只读种子；`element_features_zscore.csv` 不得被测试覆盖 |
| `ChemDB/training/` | 全局 training 目录 |
| `ChemDB/training/ckpts/` | 全局 ckpt |
| `ChemDB/training/evaluation/results/` | 评估产物 |
| `ChemDB/tmp/` | 仓库 tmp |
| `ChemDB/src/tmp/` | 历史陷阱目录 |

**允许**：读取 `ChemDB/data/metal_embedding/element_features.csv` 作为 **init-run 复制源**（只读，不写入）。

---

## 2. Logical Steps（StepRegistry + Orchestrator）

在 `app/services/step_registry.py`（或 `training_registry.py` 再合并）中新增 **`TRAINING_PIPELINE_ORDER`**，与现有 `PIPELINE_ORDER`（至 step13）分离。

完整训练 run 前置条件：**核心链已执行至 `step13`**（Tier 集成可复用 1% run 或 fixture run）。

### 2.1 步骤一览

| step_id | 执行方式 | 说明 |
|---------|----------|------|
| `build_l3_embedding_ecfp` | subprocess | 仅 ECFP，禁止 MolCLR 后端 |
| `zscore_metal_embedding` | subprocess | 显式 `--input` / `--output` |
| `prepare_training_index` | `training_manager` 内 Python | 写 `training/index.json` |
| `prepare_training_config` | `training_manager` 内 Python | 写 `training/config.yaml` |
| `training_data` | subprocess `python -m training.data` | 显式 `--training-dir`、`--index` |
| `training_train` | subprocess `python -m training.train` | 显式 `--config` |

编排入口建议：

```bash
# 仅 Phase 1C（假定 run 已有 step13 产物）
python -m app.services.orchestrator run-pipeline --run-root RUN --through training_train \
  --pipeline training

# 或从 step13 一路到 training（两段 pipeline 串联，实施时定 CLI 形态）
```

### 2.2 各步 required_inputs / expected_outputs

路径均相对于 `run_root`（registry 内用 `tmp/`、`data/`、`training/` 前缀，编排层 `resolve_run_path` 展开）。

#### （1）`build_l3_embedding_ecfp`

| 类型 | 路径 |
|------|------|
| **required_inputs** | `tmp/repaired_ligand_data.csv`、`tmp/metal_l3_index.csv`（与现网 `build_L3_embedding_index.py::load_l3_did_smiles` 一致，**不可省略**） |
| **expected_outputs** | `data/L3_embedding/L3_embedding_ecfp.npz` |
| **optional_outputs** | `data/L3_embedding/build.log` |

**命令草案**（`cwd=ChemDB`，`PYTHONPATH=ChemDB/src`）：

```text
python src/tools/build_L3_embedding_index.py \
  --tmp-dir {tmp_dir} \
  --out-dir {data_dir}/L3_embedding \
  --backends ecfp
```

**要求**：

- `--backends ecfp` 唯一；不调用 GIN/GCN。
- **代码改动**：`build_L3_embedding_index.py` 当前在模块顶层 `from molclr_api import ...`，会在 import 时拉 MolCLR。为满足 `test_no_molclr_required_for_phase1c`，需将 MolCLR **延迟导入**到 `embed_molclr` / `gin|gcn` 分支内；ECFP 路径仅依赖 RDKit。
- 输入缺失 → `status=failed`，manifest 记录 `error`，**禁止** fallback 到 `ChemDB/tmp` 或 `result_analysis`。

#### （2）`zscore_metal_embedding`

| 类型 | 路径 |
|------|------|
| **required_inputs** | `data/metal_embedding/element_features.csv` |
| **expected_outputs** | `data/metal_embedding/element_features_zscore.csv` |

**命令草案**：

```text
python data/metal_embedding/zscore_element_features.py \
  --input {data_dir}/metal_embedding/element_features.csv \
  --output {data_dir}/metal_embedding/element_features_zscore.csv
```

**要求**：保留无参数时的旧默认（`__file__.parent` 下读写）；编排**必须**显式传 run 内路径。

#### （3）`prepare_training_index`

| 类型 | 路径 |
|------|------|
| **required_inputs** | `tmp/step13_kl_nl_samples.csv`、`tmp/m_l3_pairs.csv`、`tmp/l3_gac.json`、`data/L3_embedding/L3_embedding_ecfp.npz`、`data/metal_embedding/element_features_zscore.csv` |
| **expected_outputs** | `training/index.json` |

**要求**：

- 全部由 `training_manager.prepare_training_index(ctx)` 生成；不读 `ChemDB/training/index.json`。
- 所有写入路径为 **`run_root` 绝对路径**。
- **embedding 固定 ECFP**；`ligand_embedding.backend = "ecfp"`，`mode = null`（JSON 中为 `null`）。

**与 `training/data.py::load_index_paths` 的兼容**（实施时二选一，推荐 A）：

- **方案 A（推荐）**：`index.json` 同时包含用户语义字段 + `load_index_paths` 所需嵌套结构（对 `data.py` **最小补丁**，仅解析层）：

```json
{
  "run_root": "/abs/path/to/run",
  "samples_csv": "/abs/.../tmp/step13_kl_nl_samples.csv",
  "m_l3_pairs": "/abs/.../tmp/m_l3_pairs.csv",
  "l3_gac": "/abs/.../tmp/l3_gac.json",
  "m_l3_pairs_path": "/abs/.../tmp/m_l3_pairs.csv",
  "kl_nl_ul_index": {
    "path": "/abs/.../tmp/step13_kl_nl_samples.csv",
    "format": "csv"
  },
  "ligand_embedding": {
    "backend": "ecfp",
    "mode": null,
    "path": "/abs/.../data/L3_embedding/L3_embedding_ecfp.npz",
    "dir": "/abs/.../data/L3_embedding",
    "format": "npz",
    "variants": {
      "ecfp": {
        "path": "/abs/.../data/L3_embedding/L3_embedding_ecfp.npz",
        "file": "L3_embedding_ecfp.npz"
      }
    },
    "npz_keys": ["dids", "smiles", "embeddings"]
  },
  "metal_embedding": {
    "path": "/abs/.../data/metal_embedding/element_features_zscore.csv",
    "format": "csv",
    "lookup_column": "element"
  }
}
```

- **方案 B**：重写 `load_index_paths` 仅认扁平结构（改动面更大，不推荐）。

`l3_gac`：纳入 **required_inputs 存在性检查**；若当前 `training.data` 未消费，仅作索引元数据 / 未来扩展，不在 Phase 1C 改采样逻辑。

#### （4）`prepare_training_config`

| 类型 | 路径 |
|------|------|
| **required_inputs** | `training/index.json` |
| **expected_outputs** | `training/config.yaml` |

**路径字段（均绝对路径，禁止 `null`）**：

| 字段 | 值 |
|------|-----|
| `data.training_dir` | `{run_root}/training` |
| `data.index_path` | `{run_root}/training/index.json` |
| `train.ckpt_dir` | `{run_root}/training/ckpts` |

**非路径字段**：`d_m`/`d_l` 可为 `null`；`embedding: ecfp`；集成测试覆盖 `epochs: 1`、`device: cpu`、小 `batch_size`。

模板来源：`app/templates/training_config_template.yaml`。

#### （5）`training_data`

| 类型 | 路径 |
|------|------|
| **required_inputs** | `training/index.json`；index 引用的 samples / m_l3 / ecfp npz / metal zscore 均存在 |
| **expected_outputs** | `training/split_index.json`、`training/train_records.pkl`、`training/val_records.pkl`、`training/test_records.pkl` |

**命令草案**：

```text
python -m training.data \
  --training-dir {training_dir} \
  --index {training_dir}/index.json \
  --seed 42
```

**要求**：

- `ensure_split_index` 仅在 `{training_dir}/split_index.json` 不存在时重建；**不得**读取 `ChemDB/training/` 下已有 pkl。
- `cwd`：`ChemDB`（包根，含 `training/` 包）；`PYTHONPATH` 含 `ChemDB`（与 restructured 一致）。

#### （6）`training_train`

| 类型 | 路径 |
|------|------|
| **required_inputs** | `training/config.yaml`、`training/train_records.pkl`、`training/val_records.pkl`、`training/test_records.pkl` |
| **expected_outputs** | `training/ckpts/best.pt`、`training/ckpts/history.json` |
| **optional_outputs** | `training/ckpts/last.pt` |

**命令草案**：

```text
python -m training.train --config {training_dir}/config.yaml
```

**要求**：`resolve_config` 不得因 `null` 回退到 `ChemDB/training`；生成的 config 已全部填绝对路径。

### 2.3 依赖 DAG（Phase 1C 段）

```text
[核心链] step1 → … → step12 → step13
                    │
    ┌───────────────┼───────────────┐
    v               v               │
zscore_metal    build_l3_ecfp       │（可并行）
    │               │               │
    └───────┬───────┘               │
            v                       │
    prepare_training_index          │
            v                       │
    prepare_training_config         │
            v                       │
      training_data                 │
            v                       │
      training_train                │
```

---

## 3. 需要新增的文件

| 路径 | 作用 |
|------|------|
| `app/services/training_manager.py` | `prepare_training_index(ctx)`、`prepare_training_config(ctx)`；路径校验；模板渲染 |
| `app/templates/training_index_template.json` | index 模板（占位符 `{run_root}` 等） |
| `app/templates/training_config_template.yaml` | config 模板（路径字段非 null） |
| `app/services/training_registry.py`（可选） | Phase 1C 六步 `StepSpec` + `TRAINING_PIPELINE_ORDER`；也可并入 `step_registry.py` |
| `ChemDB/training/__init__.py` | 从 `ChemDB_restructured/ChemDB/training/` 迁入最小集 |
| `ChemDB/training/data.py` | 迁入 + **最小** `load_index_paths` 兼容（若采用方案 A 可不改逻辑） |
| `ChemDB/training/train.py` | 迁入；确认 `--config` 已满足 |
| `ChemDB/training/model.py` | 迁入（`train.py` 依赖） |
| `tests/run_isolation/test_phase1c_training_run_based.py` | 默认 10 项轻量测试 |
| `tests/run_isolation/test_phase1c_integration.py`（或同文件 marker 分组） | 3 类 integration marker |
| `tests/run_isolation/fixtures/training_minimal/`（可选） | 极小 csv / 假 npz，用于无 step13 全量的单元测试 |
| `_auditing/phase_1c_implementation.md` | **实施后**填写（本文件为计划） |

**不迁入**（Phase 1C）：`training/evaluation/`、`training/results_*`、大量 legacy 脚本。

---

## 4. 需要修改的文件

| 路径 | 变更要点 |
|------|----------|
| `app/services/step_registry.py` | 注册 Phase 1C 六步；或 `from training_registry import TRAINING_*` |
| `app/services/orchestrator.py` | 支持 `training` pipeline；`prepare_*` 走 `training_manager`；增强 manifest（§6）；`run_step` 对 prepare 类无 subprocess |
| `app/services/orchestrator.init_run` | 复制 `ChemDB/data/metal_embedding/element_features.csv`（建议整目录只读种子）到 `run/data/metal_embedding/` |
| `app/services/run_context.py` | `create_dirs()` 增加 `data/L3_embedding`、`data/metal_embedding`、`training/ckpts` |
| `ChemDB/src/tools/build_L3_embedding_index.py` | MolCLR **延迟导入**；保证 `--backends ecfp` 不触发 MolCLR |
| `ChemDB/data/metal_embedding/zscore_element_features.py` | 增加 `--input`、`--output`；默认行为兼容 |
| `ChemDB/training/data.py` | 迁入后：可选增强 `load_index_paths` 读扁平字段；保持 `ensure_split_index` 行为 |
| `ChemDB/training/train.py` | 迁入后：确认 `ckpt_dir` 仅用 config 内路径（无仓库 fallback） |
| `tests/conftest.py` | 默认 markexpr 排除 `integration_embedding`、`integration_training_data`、`integration_training_train` |
| `pytest.ini` | 注册上述 3 个 integration marker |
| `tests/run_isolation/conftest.py` | 可选：repo 快照 `data/L3_embedding`、`data/metal_embedding`、`ChemDB/training`（与 `repo_tmp` 并列） |
| `_auditing/README.md` | 索引本文档；Phase 1C 状态改为「计划待确认」→ 实施后改「已完成」 |
| `ChemDB/doc/PIPELINE_IO.md` | 补充 Phase 1C 路径表（可选，实施时同步） |

### 4.1 明确不修改

| 类别 | 路径/模块 |
|------|-----------|
| UI | Streamlit 相关 |
| 存储 | SQLite 相关 |
| 评估 | `training/evaluation/*` |
| 服务 | `server.py`、`affinity_api.py`、`proxy.py` |
| MolCLR | `molclr_api.py`（不改为 Phase 1C 依赖；build 脚本仅延迟 import） |
| 算法 | RankModel、loss、step13 采样、划分比例 |
| 入口 | `main.py`（除非仅文档注释） |

---

## 5. `training_manager.py` 职责草案

```text
training_manager.py
├── prepare_training_index(ctx: RunContext) -> dict
│   ├── 校验 required_inputs 存在
│   ├── 读取 app/templates/training_index_template.json
│   ├── 替换 {run_root}、展开绝对路径
│   ├── 断言所有路径 is_relative_to(run_root)
│   └── 写入 ctx.training_dir / "index.json"
├── prepare_training_config(ctx: RunContext, *, integration: bool = False) -> dict
│   ├── 读取 training/index.json
│   ├── 读取 training_config_template.yaml
│   ├── 填充 training_dir / index_path / ckpt_dir（绝对路径）
│   ├── integration 模式：epochs=1, device=cpu, batch_size=小值
│   └── 写入 training/config.yaml
└── validate_paths_under_run(ctx, paths: list[Path]) -> None
```

编排层 `run_step("prepare_training_index")` 调用上述函数并写 manifest（无 `exit_code`，`status` 由校验结果决定）。

---

## 6. Manifest 规范（Phase 1C 增强）

在现有 `manifests/{step_id}.json` 基础上，Phase 1C **统一**为下列结构（subprocess 与 prepare 步骤均适用）：

```json
{
  "step_id": "build_l3_embedding_ecfp",
  "status": "success",
  "command": ["python", "/abs/.../build_L3_embedding_index.py", "..."],
  "cwd": "/abs/.../ChemDBWebVersion/ChemDB",
  "run_root": "/abs/.../run",
  "started_at": "2026-05-20T12:00:00+00:00",
  "finished_at": "2026-05-20T12:01:00+00:00",
  "duration_seconds": 60.5,
  "exit_code": 0,
  "inputs": [
    {
      "path": "/abs/.../run/tmp/repaired_ligand_data.csv",
      "exists": true,
      "size_bytes": 12345
    }
  ],
  "outputs": [
    {
      "path": "/abs/.../run/data/L3_embedding/L3_embedding_ecfp.npz",
      "exists": true,
      "size_bytes": 67890
    }
  ],
  "error": null
}
```

**规则**：

- `inputs` / `outputs` 中每个 `path` 必须为绝对路径。
- 凡**运行产物**必须在 `run_root` 内；`cwd` 可为 `ChemDB` 源码目录（须在 manifest 中如实记录）。
- 失败时 `status=failed`，`error` 为字符串；`exit_code` 非 0（subprocess）或省略（prepare 步骤）。
- 实施时在 `orchestrator.run_step` 中复用/抽取 `_file_stat(path)` 与 `_manifest_payload(...)`，避免六步重复逻辑。

---

## 7. 测试计划

### 7.1 默认 pytest（轻量，必须全部通过）

文件：`tests/run_isolation/test_phase1c_training_run_based.py`

| # | 测试名 | 验证点 |
|---|--------|--------|
| 1 | `test_init_run_copies_metal_embedding_seed` | init-run 后存在 `run/data/metal_embedding/element_features.csv` |
| 2 | `test_zscore_cli_writes_run_output` | 小 fixture 跑 zscore → `run/.../element_features_zscore.csv`；仓库 `ChemDB/data/metal_embedding` 无新写入 |
| 3 | `test_build_l3_embedding_ecfp_outputs_under_run` | mock 或 minimal real（少量 L3 行）→ `run/data/L3_embedding/L3_embedding_ecfp.npz`；仓库 `L3_embedding` 无新写入 |
| 4 | `test_no_molclr_required_for_phase1c` | 不设置 MolCLR 路径；import/执行 ECFP 路径不触发 `molclr_api`（配合 build 脚本延迟 import） |
| 5 | `test_prepare_training_index_uses_ecfp` | `index.json`：`backend=ecfp`、`mode=null`、path 指向 run 内 ecfp npz；所有路径在 run 内 |
| 6 | `test_prepare_training_config_paths_under_run` | `config.yaml`：`index_path`、`training_dir`、`ckpt_dir` 均在 run 内；无 `ChemDB/training` |
| 7 | `test_training_registry_steps` | Registry 含六步 id |
| 8 | `test_training_manifest_paths_under_run_mock` | mock subprocess/prepare；每步有 `logs/`、`manifests/`；manifest I/O 路径均在 run 内 |
| 9 | `test_no_write_to_repo_embedding_training` | 快照 `ChemDB/data/L3_embedding`、`ChemDB/data/metal_embedding`、`ChemDB/training` 前后一致 |
| 10 | `test_two_runs_ecfp_training_isolated` | `run_001` / `run_002` 各自 embedding、index、config、mock records/ckpt 不串 |

**默认收集**：`tests/conftest.py` 排除 `integration_embedding`、`integration_training_data`、`integration_training_train`（与 Phase 1B 做法一致）。

### 7.2 Integration pytest（不进入默认）

| Marker | 内容 | 前置 |
|--------|------|------|
| `integration_embedding` | 真实 ECFP build + 真实 zscore；产物仅在 `run/data` | run 含 `repaired_ligand_data.csv`、`metal_l3_index.csv`；或 1% 子集 run |
| `integration_training_data` | 真实 `training.data`；pkl/split 仅在 `run/training` | 同上 + ecfp npz + zscore csv |
| `integration_training_train` | 真实 `training.train`；`cpu`、`epochs=1`、小 batch；`best.pt` + `history.json` 在 `run/training/ckpts` | 完成 training_data |

**建议命令**：

```bash
conda activate scidata  # 需 RDKit；不需 MolCLR
export CHEMDB_INTEGRATION_PUBCHEM_SOURCE=...  # 若从 1% 构建 run

pytest -m integration_embedding -s ...
pytest -m integration_training_data -s ...
pytest -m integration_training_train -s ...
```

### 7.3 回归

- Phase 1A 默认 8 项 + Phase 1B 1% Tier A/B **不得回退**。
- 新增 repo 快照测试后，注意与 `assert_repo_tmp_unchanged` 共存。

---

## 8. 验收标准（与需求 §七 对齐）

| # | 标准 |
|---|------|
| 1 | 默认 `pytest` 全部通过 |
| 2 | Phase 1A / 1B 测试不回退 |
| 3 | ECFP npz 仅写入 `run_root/data/L3_embedding/` |
| 4 | metal zscore 仅写入 `run_root/data/metal_embedding/` |
| 5 | `training/index.json` 全部路径在当前 `run_root` 内 |
| 6 | `training/config.yaml` 路径字段均在 `run_root` 内且无 `null` 路径 |
| 7 | `split_index.json` 与 `*_records.pkl` 仅在 `run_root/training/` |
| 8 | `best.pt`、`history.json` 仅在 `run_root/training/ckpts/` |
| 9 | 仓库 `ChemDB/data/L3_embedding`、`ChemDB/data/metal_embedding`（运行生成）、`ChemDB/training` 无新运行产物 |
| 10 | `run_001` 与 `run_002` embedding/records/ckpt 互不覆盖 |
| 11 | Phase 1C 不依赖 MolCLR |
| 12 | 不实现 ChemBERT / MolCLR / GIN / GCN |
| 13 | 每个 logical step 有 `logs/{step_id}.log` 与 `manifests/{step_id}.json` |
| 14 | manifest 中运行产物路径均在 `run_root` 内 |

---

## 9. 建议实施顺序

| 序号 | 任务 | 产出 |
|------|------|------|
| 1 | 迁入 `ChemDB/training/{data,train,model}.py` + `pytest` 可 import | 包可加载 |
| 2 | `zscore_element_features.py` CLI + `init-run` 复制 metal 种子 + `run_context.create_dirs` | 2、9 号测试可写 |
| 3 | `build_L3_embedding_index.py` MolCLR 延迟导入 + registry `build_l3_embedding_ecfp` | 3、4 号测试 |
| 4 | `training_manager` + 两模板 + `prepare_*` registry | 5、6 号测试 |
| 5 | orchestrator manifest 增强 + `training_data` / `training_train` 接入 | 7、8 号测试 |
| 6 | `test_no_write_to_repo_*` + `test_two_runs_*` | 9、10 号测试 |
| 7 | integration markers + 3 类集成测试 | 可选验收 |
| 8 | 文档：`phase_1c_implementation.md`、更新 `_auditing/README.md`、`PIPELINE_IO.md` | 闭环 |

---

## 10. 已知风险与限制（计划阶段）

| 项 | 说明 |
|----|------|
| MolCLR 顶层 import | 必须通过延迟 import 满足 `test_no_molclr`；否则需放宽测试定义（不推荐） |
| `build_l3` 依赖 `metal_l3_index.csv` | 必须由 step13 产出；仅 repaired 不够 |
| `integration_training_*` 耗时 | 依赖完整 step13 表；可复用 Phase 1B 1% run，仍可能较长 |
| `l3_gac.json` | 纳入输入校验；训练 data 模块 Phase 1C 可不消费 |
| GIN/GCN | **刻意不做**；`phase_1c_implementation.md` 记为 future backend |
| SQLite Phase 2 | 实施前需确认：run 元数据是否入库、是否仍保留 per-run `index.json` 文件 |

---

## 11. 运行命令（实施后验证用）

```bash
cd ChemDBWebVersion
export CHEMDB_REPO_ROOT="$(pwd)/ChemDB"

# 1) 初始化 run（含 pubchem + metal 种子）
python -m app.services.orchestrator init-run \
  --run-root workspace/projects/demo/runs/run_001 \
  --pubchem-source integration_data/pubchem_1pct

# 2) 核心链至 step13（耗时；可换已有 run）
python -m app.services.orchestrator run-pipeline \
  --run-root workspace/projects/demo/runs/run_001 \
  --through step13 --stream

# 3) Phase 1C training pipeline（CLI 形态实施时确定）
python -m app.services.orchestrator run-pipeline \
  --run-root workspace/projects/demo/runs/run_001 \
  --through training_train \
  --pipeline training --stream

# 4) 默认测试
pytest
pytest tests/run_isolation/test_phase1c_training_run_based.py -v
```

---

## 12. 修订历史

| 版本 | 日期 | 说明 |
|------|------|------|
| 0.1 | 2026-05-20 | Phase 1C 实施计划初稿，待用户确认 |

---

**请确认本计划后，再按 §9 顺序实施代码；实施完成后将结果写入 `phase_1c_implementation.md`。**
