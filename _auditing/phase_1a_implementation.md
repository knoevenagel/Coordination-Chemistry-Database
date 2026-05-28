# Phase 1A 实施记录（Run 路径隔离）

| 字段 | 内容 |
|------|------|
| 阶段 | Phase 1A |
| 完成日期 | 2026-05-20 |
| 基线评估 | [initial_assessment.md](./initial_assessment.md) |
| 范围 | Run 目录隔离、编排层、step1/2/6_7/8 路径补丁、step13 去 fallback、轻量 pytest |
| 不在范围 | Streamlit、SQLite、`main.py`/服务 API、`training/` 迁入、step9/11 入 registry、端到端 integration 真跑 |

---

## 1. 目标与确认的调整项

用户确认按 Phase 1A 推进，并附加以下约束（均已落实）：

1. **目录名**：`ChemDB_run_based` → **`ChemDB`**；`CHEMDB_REPO_ROOT` 默认 `ChemDBWebVersion/ChemDB`；`pytest.ini`、orchestrator、registry、测试 fixture 统一使用 `ChemDB`。
2. **init-run**：默认只建 run 目录树 + 复制 `metal_list.txt`、`p_elements_list.txt`；**不**默认复制全量 `data/pubchem`；测试用 `tests/run_isolation/fixtures/minimal_pubchem.csv`；可选 `--pubchem-source`（默认关闭）。
3. **step13**：registry 默认传 `--kl-nl-only`；不新增 `--record-limit`；**删除**对源码树外 `result_analysis/...` 的索引 fallback，索引仅在 `output_dir` 内查找或运行时生成。
4. **测试**：默认 `pytest` 只跑轻量 + mock subprocess；真实 step 标 `@pytest.mark.integration`，需 `pytest -m integration`。
5. **step 补丁**：step1/2/6_7/8 仅路径参数，不改科研算法；CLI 保持 `./data`、`./tmp` 默认，orchestrator **显式**传入 `run_root/data`、`run_root/tmp`。
6. **StepRegistry**：step8 `expected_outputs` 仅核心 neo4j 子集；stats 等放入 `optional_outputs`。

---

## 2. 实施过程（时间线摘要）

| 步骤 | 动作 |
|------|------|
| 1 | 将 `ChemDBWebVersion/ChemDB_run_based` 重命名为 `ChemDB` |
| 2 | 新增 `app/services/`：`run_context.py`、`step_registry.py`、`orchestrator.py`、`__main__.py` |
| 3 | 补丁 `ChemDB/src/step1.py`、`step2.py`、`step6_7.py`、`step8.py`（目录 CLI） |
| 4 | 补丁 `ChemDB/src/step_13_complete.py`：移除 `result_analysis` 索引回退 |
| 5 | 新增 `pytest.ini`、`tests/run_isolation/*`、`workspace/.gitkeep` |
| 6 | 验收：`pytest tests/run_isolation -m "not integration"` → 8 passed |

---

## 3. 目录与命名变更

| 变更前 | 变更后 |
|--------|--------|
| `ChemDBWebVersion/ChemDB_run_based/` | `ChemDBWebVersion/ChemDB/` |

编排与测试中 `WEB_ROOT / "ChemDB"`、`DEFAULT_CHEMDB_REPO`、`CHEMDB_ROOT` fixture 均指向新路径。

**说明**：仓库内 `ChemDB/data/pubchem/` 仍保留全量元素 CSV（约 1.8G），供本地真实运行；`init-run` **不会**自动复制该目录，仅在使用 `--pubchem-source` 时写入 run。

---

## 4. 新增文件清单

| 路径 | 作用 |
|------|------|
| `app/__init__.py` | 包占位 |
| `app/services/__init__.py` | 包占位 |
| `app/services/run_context.py` | `RunContext`：`run_root`、`data/`、`tmp/`、`training/`、`logs/`、`reports/`、`manifests/`；`create_dirs()`、`resolve_under_run()` |
| `app/services/step_registry.py` | `STEP_REGISTRY`：step1→step2→step4_5→step6_7→step8→step12→step13 |
| `app/services/orchestrator.py` | `init-run`、`run-step`、`run-pipeline`；输入/输出检查；log + manifest |
| `app/services/__main__.py` | `python -m app.services.orchestrator` |
| `pytest.ini` | `pythonpath = . ChemDB/src`；`integration` marker |
| `tests/run_isolation/conftest.py` | 会话级守护：`ChemDB/tmp` 不变、禁止 `ChemDB/src/tmp` |
| `tests/run_isolation/test_run_context.py` | Run 路径解析 |
| `tests/run_isolation/test_step_registry.py` | registry 路径均在 `run_root` 下 |
| `tests/run_isolation/test_two_runs.py` | 双 run 互不覆盖 |
| `tests/run_isolation/test_no_write_to_repo_tmp.py` | 编排 mock 不写 repo tmp |
| `tests/run_isolation/test_orchestrator.py` | init-run、manifest、mock subprocess |
| `tests/run_isolation/test_integration.py` | `@pytest.mark.integration`，默认 skip |
| `tests/run_isolation/fixtures/minimal_pubchem.csv` | 测试用最小 pubchem |
| `tests/run_isolation/fixtures/metal_list.txt` | 测试用金属列表 |
| `tests/run_isolation/fixtures/p_elements_list.txt` | 测试用 P 元素列表 |
| `workspace/.gitkeep` | run 工作区占位 |
| `_auditing/README.md` | 本目录索引与记档约定 |
| `_auditing/phase_1a_implementation.md` | 本文 |

---

## 5. 修改文件清单（ChemDB 源码）

| 文件 | 修改性质 |
|------|----------|
| `ChemDB/src/step1.py` | 新增 `--pubchem-dir`、`--metal-list`、`--p-elements-list`、`--tmp-dir` 等；保留无 `-` 时的 legacy positional |
| `ChemDB/src/step2.py` | 新增 `--input-dir`、`--output-dir` |
| `ChemDB/src/step6_7.py` | 新增 `--input-dir`、`--output-dir`；`ga_path`/`irl_path` 由 `input_dir` 推导 |
| `ChemDB/src/step8.py` | 新增 `--input-dir`、`--output-dir`、`--metal-list`；`IRL_filtered.csv` 改为 `input_dir / "IRL_filtered.csv"`（修复裸 `Path("tmp/...")`） |
| `ChemDB/src/step_13_complete.py` | 删除对 `parent.parent/result_analysis/.../metal_l3_index.csv` 的 fallback；默认仅 `output_dir/metal_l3_index.csv` |

**未改**：`step4_5.py`、`step12.py`（已有目录 CLI）；`main.py`、`server.py`、step9/10/11。

---

## 6. 编排层行为说明

### 6.1 Run 目录布局

```
{run_root}/
├── data/
│   ├── metal_list.txt          # init-run 从 ChemDB/data 复制
│   ├── p_elements_list.txt
│   └── pubchem/                # 仅 --pubchem-source 或测试 fixture 填充
├── tmp/                        # 各 step 产物
├── training/                   # 目录已创建，Phase 1A 未迁入内容
├── logs/{step_id}.log
├── manifests/
│   ├── run.json                # init-run
│   ├── {step_id}.json          # 每步
│   └── pipeline.json           # run-pipeline 汇总
└── reports/
```

### 6.2 子进程环境

- `cwd` = `CHEMDB_REPO_ROOT`（默认 `ChemDB/`）
- `PYTHONPATH` = `{CHEMDB_REPO_ROOT}/src`
- `CHEMDB_RUN_ROOT` = 当前 `--run-root`
- 命令行通过 `{data_dir}`、`{tmp_dir}` 占位符展开为 run 内绝对路径

### 6.3 CLI

```bash
cd ChemDBWebVersion

# 初始化（可选 pubchem）
python -m app.services.orchestrator init-run \
  --run-root workspace/projects/p001/runs/run_001 \
  --pubchem-source tests/run_isolation/fixtures/minimal_pubchem.csv

python -m app.services.orchestrator run-step --run-root ... --step step2
python -m app.services.orchestrator run-pipeline --run-root ... --through step13
```

`--run-root` 挂在各子命令上（`init-run` / `run-step` / `run-pipeline`），无需写在全局最前。

### 6.4 失败策略

- **输入缺失**：不写 log 跑 step，直接 `status=failed`，写 `manifests/{step}.json`，`error=missing required inputs`。
- **进程非零退出** 或 **mandatory 输出缺失**：`status=failed`，manifest 记录 `exit_code`、`output_check`。
- **step13**：不再读取 repo 外 `result_analysis`；缺 `tmp/` 内 step12 产物时由输入检查或 step 自身失败。

### 6.5 StepRegistry 要点

| Step | mandatory outputs（节选） |
|------|---------------------------|
| step8 | 6 个核心 `neo4j_*` CSV + 关系表；其余 neo4j 表、`step8_stats.json` 为 optional |
| step13 | `tmp/step13_kl_nl_samples.csv`；stats、`metal_l3_index.csv` 为 optional |
| step13 argv | 固定含 `--kl-nl-only` |

---

## 7. 测试策略

| 类型 | 命令 | 内容 |
|------|------|------|
| 默认（轻量） | `pytest tests/run_isolation -m "not integration"` | RunContext、registry 路径、双 run 隔离、repo tmp 守护、orchestrator mock |
| 集成 | `pytest tests/run_isolation -m integration` | 真实 subprocess 跑 step（需 RDKit 等；当前 `test_integration.py` 为 skip 占位） |

**默认测试覆盖点**（用户要求）：

- RunContext 路径
- repo `ChemDB/tmp` 不被污染
- 不创建 `ChemDB/src/tmp`
- 两个 run 互不覆盖
- StepRegistry required/expected 路径落在 `run_root` 内
- orchestrator 能写 log 和 manifest

**验收命令（Phase 1A）**：

```bash
cd /home/chenghua/ZhuoLi/ChemDBWebVersion
python -m pytest tests/run_isolation -q -m "not integration"
# 8 passed, 1 deselected (integration)
```

---

## 8. 已知限制与后续阶段建议

| 项 | 状态 |
|----|------|
| `ChemDB/training/` | 未迁入；`index.json` / per-run `config.yaml` → Phase 1 后续或 Phase 2 |
| step9 / step11 | 未纳入 `STEP_REGISTRY` |
| 端到端 `run-pipeline --through step13` | 未在 CI/默认测试中真跑；需 integration 或手动 |
| `PIPELINE_IO.md` | 仍位于 `ChemDB/doc/`，可补 Run 物理路径一节 |
| Streamlit + SQLite | 明确不在 Phase 1A |

**建议 Phase 1B 记档时更新**：`phase_1b_implementation.md`，并在此 README 索引表中登记。

---

## 9. 修订历史

| 版本 | 日期 | 说明 |
|------|------|------|
| 1.0 | 2026-05-20 | Phase 1A 落地记录初稿 |
