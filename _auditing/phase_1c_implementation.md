# Phase 1C 实施记录

| 字段 | 值 |
|------|-----|
| 状态 | 已完成（pytest + **run_demo 1% 端到端实测**） |
| 计划 | [phase_1c_plan.md](./phase_1c_plan.md) |
| 日期 | 2026-05-21（实施）；2026-05-21（run_demo 补记） |
| 测试环境 | conda **`scidata`**（RDKit + torch + CUDA）；入口 `./scripts/run_pytest.sh` |

## 目标

- ECFP-only L3 embedding、metal zscore、per-run `training/index.json` + `config.yaml`
- `training.data` / `training.train` 全部读写限定在 `{run_root}` 内
- 编排支持 `--pipeline training --through training_train`
- MolCLR 延迟导入；ECFP 路径不依赖 `molclr_api`

## 不在范围内

Streamlit、SQLite、API、evaluation、recommendation、Model Registry、MolCLR/ChemBERT/GIN/GCN 服务、RankModel/loss/sampling 算法改动。

## 路径确认（实施前）

| 资源 | 路径 |
|------|------|
| `build_L3_embedding_index.py` | `ChemDB/src/tools/build_L3_embedding_index.py` |
| `zscore_element_features.py` | `ChemDB/data/metal_embedding/zscore_element_features.py` |
| training 包 | 从 `ChemDB_restructured/ChemDB/training/` 迁入 `ChemDB/training/{__init__.py,data.py,train.py,model.py}` |

## 新增文件

| 路径 | 说明 |
|------|------|
| `app/services/training_manager.py` | `prepare_training_index` / `prepare_training_config` |
| `app/templates/training_index_template.json` | index 模板（方案 A） |
| `app/templates/training_config_template.yaml` | config 模板 |
| `ChemDB/training/__init__.py` | 训练包 |
| `ChemDB/training/data.py` | 数据划分与 Dataset |
| `ChemDB/training/train.py` | 训练入口 |
| `ChemDB/training/model.py` | RankModel |
| `tests/run_isolation/test_phase1c_training_run_based.py` | Phase 1C 默认 + integration 测试 |
| `scripts/run_pytest.sh` | 统一在 `scidata` 下执行 pytest |
| `TESTING.md` | 测试分层与 run_demo 操作说明 |

## 修改文件

| 路径 | 说明 |
|------|------|
| `app/services/run_context.py` | `create_dirs` 补齐 L3/metal/training/ckpts |
| `app/services/orchestrator.py` | init-run metal 种子；training pipeline；增强 manifest；`PYTHONPATH`；in-process steps |
| `app/services/step_registry.py` | 六步 training registry + `TRAINING_PIPELINE_ORDER` |
| `ChemDB/data/metal_embedding/zscore_element_features.py` | `--input` / `--output` CLI |
| `ChemDB/src/tools/build_L3_embedding_index.py` | MolCLR 延迟导入 |
| `ChemDB/training/train.py` | `resolve_training_device()`；禁止静默回退 CPU；配置 `cpu` 报错 |
| `pytest.ini` | integration_embedding / training_* markers |
| `tests/conftest.py` | 默认排除 integration markers；**强制 scidata 环境** |
| `tests/conftest.py` | 默认排除新 integration markers |
| `tests/run_isolation/conftest.py` | `mock_run_subprocess` helper |
| `tests/run_isolation/test_orchestrator.py` | init-run metal 断言；mock 修复 |
| `tests/run_isolation/test_no_write_to_repo_tmp.py` | mock 修复 |
| `tests/run_isolation/test_step_registry.py` | 允许 `training/` 路径前缀 |

## 行为变更

### init-run

- 只读复制 `ChemDB/data/metal_embedding/element_features.csv` → `{run_root}/data/metal_embedding/`
- 源缺失时 `run.json` 记录 `metal_embedding_seed_error` 并 `raise FileNotFoundError`

### run-pipeline

```bash
python -m app.services.orchestrator run-pipeline \
  --run-root RUN --pipeline training --through training_train
```

- `prepare_training_index` / `prepare_training_config`：in-process，写 manifest + log
- 其余步骤：subprocess，`cwd=ChemDB`，`PYTHONPATH=ChemDB:ChemDB/src`
- 集成/冒烟训练：`CHEMDB_TRAINING_INTEGRATION=1` → 单组超参、**cuda**、epochs=3、early_stop=0、batch_size=4（不做超参扫描）
- 训练设备：`CHEMDB_TRAINING_DEVICE` 默认 `cuda`（可 `cuda:0`）；`train.py` 拒绝 `device: cpu`，无 GPU 时显式失败

### Manifest（Phase 1C）

每步 `manifests/{step_id}.json` 含：`inputs[]` / `outputs[]`（绝对路径、`exists`、`size_bytes`）、`duration_seconds`、`cwd`、`run_root`、`error`。

## 验收

### 默认 pytest

**须在 `scidata` 环境中运行**（见仓库根目录 [TESTING.md](../TESTING.md)）：

```bash
cd ChemDBWebVersion
./scripts/run_pytest.sh tests/ -q
```

结果（scidata，2026-05-21）：**19 passed, 9 deselected**（约 1s；含真实 ECFP 单测）

- Phase 1C：`test_phase1c_training_run_based.py` 默认项通过
- Phase 1A/1B 相关测试未回退（orchestrator mock 改为 patch `_run_subprocess`，兼容 `--stream`）

### Integration pytest（可选）

```bash
./scripts/run_pytest.sh tests/run_isolation/test_phase1c_training_run_based.py -m integration_training_train -v
```

CI/本机未单独跑全量 integration marker 时，以 **run_demo 手动 E2E** 作为 training 真实路径验收（见下节）。

---

## run_demo 1% 端到端实测（2026-05-21）

**run_root**：

```text
/home/chenghua/ZhuoLi/ChemDBWebVersion/workspace/projects/p001/runs/run_demo
```

**数据**：`init-run --pubchem-source integration_data/pubchem_1pct`（约 80 个金属 CSV，非 `minimal_pubchem.csv`）。

**说明**：pytest 1% 集成写在 `/tmp/pytest-of-.../`，**不会**写入 `run_demo`；下列为 orchestrator 手动跑结果。

### 核心链 `--pipeline core --through step13`

| 项 | 值 |
|----|-----|
| 结果 | `pipeline_core.json` → `"ok": true`，7 步均 `success` |
| 墙钟 | step1 `06:15:33` UTC → 结束 `06:31:43` UTC，约 **16 分 10 秒** |
| step1 占比 | `duration_seconds` ≈ **863.7 s**（约 14.4 min） |
| `tmp/` | 38 文件，约 **134 MB** |

training 所需 tmp（均已生成）：

| 文件 | 大小（约） |
|------|------------|
| `step13_kl_nl_samples.csv` | 1.9M |
| `m_l3_pairs.csv` | 636K |
| `l3_gac.json` | 800K |
| `repaired_ligand_data.csv` | 4.8M |
| `metal_l3_index.csv` | 33M |

### Training 链 `--pipeline training --through training_train`

环境：`CHEMDB_TRAINING_INTEGRATION=1`（冒烟：单组超参、3 epoch、GPU、关闭早停）。

| 项 | 值 |
|----|-----|
| 结果 | `pipeline_training`（manifest 汇总）→ `"ok": true`，6 步均 `success` |
| 结束时间 | `2026-05-21T06:43:08Z` |
| 墙钟（training 六步合计） | manifest `duration_seconds` 合计约 **15.6 s** |

| step_id | duration（约） |
|---------|----------------|
| build_l3_embedding_ecfp | 9.0 s |
| zscore_metal_embedding | 0.4 s |
| prepare_training_index / config | &lt;0.01 s |
| training_data | 1.3 s |
| training_train | 4.8 s |

**config.yaml（冒烟）**：`device: cuda`，`epochs: 3`，`batch_size: 4`，`early_stop_patience: 0`，`scheduler: null`。

**主要产物**（均在 `run_root` 内，未污染 `ChemDB/training` 或 `ChemDB/data/L3_embedding` 运行目录）：

| 路径 | 说明 |
|------|------|
| `data/L3_embedding/L3_embedding_ecfp.npz` | ~211M |
| `data/metal_embedding/element_features_zscore.csv` | zscore 输出 |
| `training/index.json`, `config.yaml` | 绝对路径索引 |
| `training/split_index.json`, `*_records.pkl` | training.data |
| **`training/ckpts/best.pt`** | **~226M，验收用模型** |
| `training/ckpts/last.pt` | 最后一轮 |
| `training/ckpts/history.json` | **3 条 epoch** 记录 |

**history 摘要**（冒烟不考核指标，仅验证能训完）：

- epoch 1–3：`train_loss` 29.4 → 6.9 → 1.9；`val_loss` 波动较大属正常。

**复现命令**：

```bash
conda activate scidata
cd ChemDBWebVersion
RUN=workspace/projects/p001/runs/run_demo
PCT="$(pwd)/integration_data/pubchem_1pct"

python -m app.services.orchestrator init-run --run-root "$RUN" --pubchem-source "$PCT"
python -m app.services.orchestrator run-pipeline --run-root "$RUN" --pipeline core --through step13 --stream

export CHEMDB_TRAINING_INTEGRATION=1
python -m app.services.orchestrator run-pipeline --run-root "$RUN" --pipeline training --through training_train --stream
```

## 已知限制

1. **training 依赖 torch + CUDA**：无 GPU 时 `training_train` 失败（设计如此）；编排 mock 测可不装 GPU。
2. **ECFP / 核心链仍主要在 CPU**（RDKit、多进程 step）；仅 **RankModel 训练** 强制 GPU。
3. **冒烟 3 epoch 不表示模型可用**：1% 数据 + 少轮次仅验证路径与 `best.pt` 产出。
4. **核心链与 training 链分离**：需先 `core` 至 step13 再 `--pipeline training`。
5. **pytest 与 run_demo 分离**：Tier A+B 只到 step4_5；全链 + 模型需按上节手动命令在固定 `run_root` 执行。

## Future backends

| Backend | Mode | 状态 |
|---------|------|------|
| `ecfp` | — | Phase 1C 已实现 |
| `molclr` | `gin` / `gcn` | 未来接入；`build_L3_embedding_index.py` 已延迟 import |
| `chembert` | — | 未来接入 |

## Phase 2（SQLite）前待确认

1. run 元数据是否入库，是否仍保留 per-run `training/index.json` 文件。
2. embedding / records 大文件是否只存路径引用。
3. 多 run 并发写库与文件锁策略。
