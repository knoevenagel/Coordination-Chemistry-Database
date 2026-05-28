# ChemDBWebVersion 测试指南

**所有测试应在 `scidata` conda 环境中运行**（含 RDKit、torch 等依赖）。推荐入口：

```bash
cd ChemDBWebVersion
./scripts/run_pytest.sh          # 默认轻量测试
./scripts/run_pytest.sh -v       # 详细输出
```

等价于：

```bash
conda activate scidata
cd ChemDBWebVersion
python -m pytest
```

若未在 `scidata` 中直接调用 `pytest`，收集阶段会报错退出（可用 `CHEMDB_ALLOW_NON_SCIDATA_TEST=1` 临时跳过，不推荐）。

---

## 1. 默认测试（日常 CI / 提交前）

不含 integration，约 1 秒内完成（视机器而定）。

```bash
./scripts/run_pytest.sh
# 或
./scripts/run_pytest.sh tests/ -q
```

覆盖：Run 隔离、编排 manifest、Phase 1C 轻量 mock、registry 等。

---

## 2. Phase 1C 测试

```bash
./scripts/run_pytest.sh tests/run_isolation/test_phase1c_training_run_based.py -v
```

| 类型 | 命令 |
|------|------|
| 默认 10 项 | `./scripts/run_pytest.sh tests/run_isolation/test_phase1c_training_run_based.py` |
| 真实 ECFP + zscore | `./scripts/run_pytest.sh tests/run_isolation/test_phase1c_training_run_based.py -m integration_embedding -v` |
| 真实 training.data | `... -m integration_training_data -v` |
| 真实 training.train（GPU, 3 epoch 冒烟） | `... -m integration_training_train -v` |

---

## 3. Phase 1B 集成测试（1% PubChem）

先生成或确认子集，并设置环境变量：

```bash
conda activate scidata
cd ChemDBWebVersion

# 若尚未生成 1% 数据：
python scripts/generate_pubchem_1pct_subset.py \
  --source ChemDB/data/pubchem \
  --output integration_data/pubchem_1pct

export CHEMDB_INTEGRATION_PUBCHEM_SOURCE="$(pwd)/integration_data/pubchem_1pct"
```

**Phase 1B 最低验收（Tier A+B，约 20–25 分钟）：**

Tier B 在 Phase 1E 后改为：**step1+step2 → bind GA fixture → apply-ga**（不再 `through step4_5` 自动生成 GA）。

```bash
./scripts/run_pytest.sh -s tests/run_isolation/test_integration_1pct.py -v \
  -m "integration_1pct_tier_a or integration_1pct_tier_b"
```

测试会把 `tests/run_isolation/fixtures/ga_sets/fixture_run_demo` 安装到 `workspace/ga_sets/` 再 `bind-ga-version-to-run`。

其他 marker 见 [integration_data/README.md](integration_data/README.md)。

---

## 4. Phase 1E（GA workspace + binding）

```bash
./scripts/run_pytest.sh tests/run_isolation/test_phase1e_ga_workspace.py -v
```

覆盖：GA version 不可变、bind 物化、checksum、无绑定 pipeline 失败、不自动生 GA、stale report、不写 `ChemDB/tmp`。

手动：

```bash
python -m app.services.ga_registry_manager bind-ga-version-to-run \
  --run-root <run> --ga-set-id fixture_run_demo --ga-version-id v001
python -m app.services.ga_registry_manager apply-ga-to-run --run-root <run>
```

---

## 5. Phase 1A 集成（最小 fixture）

```bash
./scripts/run_pytest.sh -m integration -v
```

---

## 6. 手动跑 training pipeline（非 pytest）

需先完成核心链至 step13，或使用已有 `run_root`：

```bash
conda activate scidata
cd ChemDBWebVersion

python -m app.services.orchestrator init-run --run-root /path/to/run
# 核心链：step1→step2 → bind GA version → apply-ga-to-run → orchestrator through step13

python -m app.services.orchestrator run-pipeline \
  --run-root /path/to/run \
  --pipeline training \
  --through training_train \
  --stream
```

集成/冒烟训练（**单组超参**，不做扫描；默认 **3 个 epoch** / **GPU (cuda)** / 关闭早停）：

```bash
export CHEMDB_TRAINING_INTEGRATION=1
# 可选：export CHEMDB_TRAINING_INTEGRATION_EPOCHS=3   # 默认已是 3
python -m app.services.orchestrator run-pipeline \
  --run-root /path/to/run --pipeline training --through training_train
```

`prepare_training_config(integration=True)` 写入的 `config.yaml` 要点：`epochs: 3`、`early_stop_patience: 0`、`scheduler: null`、`batch_size: 4`、`device: cuda`。跑完后可看 `training/ckpts/history.json` 是否至少有 3 条 epoch 记录。

---

## 7. 环境变量速查

| 变量 | 用途 |
|------|------|
| `CHEMDB_INTEGRATION_PUBCHEM_SOURCE` | 1% PubChem 目录（integration_1pct 必填） |
| `CHEMDB_STREAM_STEP_OUTPUT=1` | 子进程日志打到终端（integration marker 自动开启） |
| `CHEMDB_TRAINING_INTEGRATION=1` | training config：单组超参 / **cuda** / epochs=3 / 关闭早停 / 小 batch |
| `CHEMDB_TRAINING_INTEGRATION_EPOCHS` | 冒烟训练 epoch 数（默认 `3`） |
| `CHEMDB_TRAINING_DEVICE` | 训练设备（默认 `cuda`，可设 `cuda:0`） |
| `CHEMDB_ALLOW_NON_SCIDATA_TEST=1` | 跳过 scidata 环境检查（仅调试） |
| `CHEMDB_TEST_CONDA_ENV` | `run_pytest.sh` 使用的 env 名（默认 `scidata`） |

---

## 8. 实施记录

- Phase 1A：[\_auditing/phase_1a_implementation.md](_auditing/phase_1a_implementation.md)
- Phase 1B：[\_auditing/phase_1b_implementation.md](_auditing/phase_1b_implementation.md)
- Phase 1C：[\_auditing/phase_1c_implementation.md](_auditing/phase_1c_implementation.md)
