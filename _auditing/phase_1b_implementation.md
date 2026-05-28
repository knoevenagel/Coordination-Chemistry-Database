# Phase 1B 实施记录（Integration 数据集与真实 step 测试）

| 字段 | 内容 |
|------|------|
| 阶段 | Phase 1B |
| 完成日期 | 2026-05-20 |
| 前置 | [phase_1a_implementation.md](./phase_1a_implementation.md) |
| 范围 | 1% PubChem 子集生成、integration pytest、Tier A/B 最低验收 |
| 不在范围 | Streamlit、SQLite、API 服务、科研算法、training（→ Phase 1C）、step8/12/13 必过 |

---

## 1. 目标

1. **保留** `tests/run_isolation/fixtures/minimal_pubchem.csv`，默认轻量测试不变。
2. **新增** `integration_data/pubchem_1pct/` 作为 1% 子集输出目录；大 CSV **不提交 git**。
3. `init-run --pubchem-source` **复制** CSV 到 `run_root/data/pubchem/`；真实 step1 只读 run 内路径（orchestrator 传 `--pubchem-dir {data_dir}/pubchem`）。
4. **最低验收**：Tier A（init-run + 真 step1 + step2）、Tier B（`run-pipeline --through step4_5`）。
5. 可选：Tier C（through step6_7）；enhanced（step8 等）**非** 1B 必过。
6. 子集：**按元素 CSV 分层抽样**，`seed=42`，`MANIFEST.json`；`output_sha256` 必填；`source_sha256` 可选（`--skip-source-checksum`）。
7. `integration_1pct_enhanced` **不** 带父 marker `integration_1pct`，避免 `pytest -m integration_1pct` 误跑 enhanced。

---

## 2. 新增文件

| 路径 | 作用 |
|------|------|
| `scripts/generate_pubchem_1pct_subset.py` | 生成 1% 子集 + MANIFEST + README |
| `integration_data/README.md` | 生成与 pytest 用法 |
| `integration_data/pubchem_1pct/` | 生成产物（gitignore `*.csv`） |
| `.gitignore` | 忽略 `integration_data/pubchem_1pct/*.csv` |
| `tests/run_isolation/integration_helpers.py` | env 解析、init-run、manifest 断言 |
| `tests/run_isolation/test_integration.py` | `@pytest.mark.integration`，minimal 真跑 step1+2 |
| `tests/run_isolation/test_integration_1pct.py` | Tier A/B/C + enhanced |
| `_auditing/phase_1b_implementation.md` | 本文 |

---

## 3. 修改文件

| 路径 | 变更 |
|------|------|
| `pytest.ini` | 新增 markers |
| `tests/conftest.py` | 未指定 `-m` 时默认排除 integration 系列（与 `-m tier_a` 不冲突） |
| `_auditing/README.md` | 索引 Phase 1B；Phase 1C 说明（training） |

**Phase 1B 验收后补充（step1 并行收尾，2026-05-20）**：

| 路径 | 变更 |
|------|------|
| `ChemDB/src/step1.py` | worker 异常必回传结果；主循环 stall 日志 / worker 全退出时的收尾；`CHEMDB_STEP1_SORT_BY_LENGTH`、`CHEMDB_STEP1_ROW_TIMEOUT_SEC` |
| `ChemDB/src/step2.py` | 可选 `CHEMDB_STEP2_WORKERS`（与 step1 一致；**测试不设置**，用默认 worker 数） |
| `tests/run_isolation/conftest.py` | 集成测试仅设 stream / sort / row timeout，**不设** worker 数 |
| `integration_data/README.md` | 实测耗时与 tqdm 尾部现象说明 |

---

## 4. 环境变量

| 变量 | 说明 |
|------|------|
| `CHEMDB_INTEGRATION_PUBCHEM_SOURCE` | 指向 1% 目录（含 `MANIFEST.json` 与 `*.csv`）或单 CSV；未设置则 skip 全部 `integration_1pct_*` |
| `CHEMDB_STREAM_STEP_OUTPUT` | 集成测试 conftest 自动设为 `1`；`-s` 时可在终端看到 step 子进程输出 |
| `CHEMDB_STEP1_SORT_BY_LENGTH` | 集成测试自动 `1`：按 SMILES 长度降序，减轻 tqdm 在 99% 处“假卡住” |
| `CHEMDB_STEP1_ROW_TIMEOUT_SEC` | 集成测试自动 `600`：单行超过 10 分钟则跳过（`0` 表示关闭） |
| `CHEMDB_STEP1_WORKERS` / `CHEMDB_STEP2_WORKERS` | 可选覆盖默认 `min(cpu, 200)`；**集成测试不设置**（见 §10） |

---

## 5. 生成 1% 子集

```bash
cd ChemDBWebVersion
python scripts/generate_pubchem_1pct_subset.py --skip-source-checksum
export CHEMDB_INTEGRATION_PUBCHEM_SOURCE="$(pwd)/integration_data/pubchem_1pct"
```

---

## 6. Pytest 用法

| 命令 | 说明 |
|------|------|
| `pytest` | 默认：仅 Phase 1A 轻量（8 tests） |
| `pytest -m integration` | minimal 真 step1+2（需 RDKit） |
| `pytest -m "integration_1pct_tier_a or integration_1pct_tier_b"` | **Phase 1B 最低验收** |
| `pytest -m integration_1pct_tier_c` | 可选 through step6_7 |
| `pytest -m integration_1pct_enhanced` | 可选 step8（需先跑到 step6_7） |

**排除 enhanced**（若将来某测试同时打了两个 marker）：

```bash
pytest -m "integration_1pct and not integration_1pct_enhanced"
```

**不要**依赖裸 `pytest -m integration_1pct` 作为唯一 CI 入口（Tier C 无父 marker；enhanced 无父 marker）。

---

## 7. 测试矩阵

| 测试 | Marker | 数据 | 行为 |
|------|--------|------|------|
| `test_integration_minimal_step1_step2` | `integration` | minimal fixture | 真 step1+2 |
| `test_1pct_tier_a_*` | `integration_1pct` + `tier_a` | 1% via env | init-run 复制 + step1+2 |
| `test_1pct_tier_b_*` | `integration_1pct` + `tier_b` | 1% | pipeline → step4_5 |
| `test_1pct_tier_c_*` | `tier_c` only | 1% | pipeline → step6_7 |
| `test_1pct_enhanced_step8` | `enhanced` only | 1% | pipeline → step6_7 + step8 |

共用 `conftest`：`ChemDB/tmp` 不变、`ChemDB/src/tmp` 不存在；集成 marker 下自动 `CHEMDB_STREAM_STEP_OUTPUT=1`、`CHEMDB_STEP1_SORT_BY_LENGTH=1`、`CHEMDB_STEP1_ROW_TIMEOUT_SEC=600`（**不**覆盖 worker 数）。

---

## 8. Phase 1B 验收结果与 step1 尾部耗时（2026-05-20）

### 8.1 Pytest 最低验收（通过）

环境：`scidata`（RDKit）、`CHEMDB_INTEGRATION_PUBCHEM_SOURCE=integration_data/pubchem_1pct`、`-s` 查看子进程输出。

```bash
pytest -s tests/run_isolation/test_integration_1pct.py -v \
  -m "integration_1pct_tier_a or integration_1pct_tier_b"
```

实测（默认 worker，`min(cpu, 200)`）：

```
3 passed, 2 deselected in 1269.72s (0:21:09)
```

即 Tier A（init-run + step1 + step2 ×2）+ Tier B（pipeline through step4_5）合计约 **21 分钟**；与“只等 6～7 分钟”的预期不符，但属正常完成而非 pytest/orchestrator 死锁。

### 8.2 现象：tqdm 在 19436/19437 长时间不动

1% 去重后约 **19437** 条。step1 并行 tqdm 常在 **~50s 内到 99%**，随后在 **19435→19437** 再耗 **数分钟至十余分钟**（取决于机器与当批“慢 SMILES”落在哪条 worker）。

| 验证方式 | 结论 |
|----------|------|
| 直接跑 `ChemDB/src/step1.py` + 1% 数据（非 pytest） | 同样会在末尾停住；约 **6～7+ 分钟** 可结束（视慢行而定） |
| `pytest` + `orchestrator.run_step` | 行为一致；**不是** pytest 独有 |
| 与 `ChemDB_restructured` 对比 | 并行/队列逻辑相同；差异在路径 CLI、CSV 列名、`CHEMDB_STEP1_WORKERS` 等 Phase 1A 补丁 |

**团队结论**：部分 worker 上的个别结构在 RDKit / `ChemicalComplex` 上**确实需要很长时间**，属数据与算法特性，不是集成层引入的回归。

### 8.3 曾采取的排查与代码加固

1. **怀疑** `CHEMDB_STEP1_WORKERS=8` 可缩短尾部等待 → 集成测试曾强制 8 worker；验收通过后按团队意见**改回默认**，不再在 conftest 里设 worker 数。
2. **`step1.py` 加固**（保留）：
   - worker 外层 `except` 后仍 `result_queue.put((None, True))`，避免主进程永久 `queue.get`；
   - 全部 worker 退出但 `completed_tasks < total` 时打 error 并 break，避免真死锁；
   - 无进展每 120s 打 `stall` 日志；
   - 集成测试启用 `CHEMDB_STEP1_SORT_BY_LENGTH=1`（长 SMILES 优先，减少“最后一条 straggler”拖住 tqdm）；
   - 集成测试启用 `CHEMDB_STEP1_ROW_TIMEOUT_SEC=600`（极端结构最多等 10 分钟）。
3. **孤立测试**：曾在 `/tmp/chemdb_step1_1pct_isolated_*` 下 nohup 跑 step1；轮询约 6.5 min 时仍在 19436/19437，后进程正常退出；临时目录已清理。

### 8.4 决策记录

| 项 | 决定 |
|----|------|
| Phase 1B 最低验收 | **通过**（Tier A + Tier B，约 21 min） |
| 集成测试 worker 数 | **不覆盖**，与生产一致 `min(cpu, 200)` |
| 集成测试保留的环境项 | `CHEMDB_STREAM_STEP_OUTPUT`、`CHEMDB_STEP1_SORT_BY_LENGTH`、`CHEMDB_STEP1_ROW_TIMEOUT_SEC` |
| CI / 本地等待时间 | 建议按 **≥25 分钟** 规划 1% Tier A+B；不要以 step1 tqdm 到 99% 作为结束信号 |

---

## 9. Phase 1C 预告（training，非 Phase 2）

- per-run `training/index.json`、`config.yaml`、`records.pkl`、`ckpts` 隔离
- 完成后再进入 SQLite / Streamlit

---

## 10. 修订历史

| 版本 | 日期 | 说明 |
|------|------|------|
| 1.0 | 2026-05-20 | Phase 1B 初版 |
| 1.1 | 2026-05-20 | §8 验收结果、step1 尾部耗时与团队结论；step1/step2/conftest 补充；集成测试 worker 改回默认 |
