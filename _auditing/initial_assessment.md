# ChemDB Run-based 改造初始评估（Initial Assessment）

> **状态**：Phase 1A 已落地（2026-05-20）。实施细节见 [phase_1a_implementation.md](./phase_1a_implementation.md)；后续阶段变更见 [_auditing/README.md](./README.md)。

| 字段 | 内容 |
|------|------|
| 审计对象（基线扫描时） | `ChemDBWebVersion/ChemDB_run_based` → **已重命名为** `ChemDBWebVersion/ChemDB` |
| 目标架构 | `ChemDBWebVersion/` + `ChemDB/` + `app/` + `workspace/projects/.../runs/run_001/` |
| 当前阶段范围 | Run 路径隔离；**不含** Streamlit、SQLite |
| 关联文档 | `ChemDB/doc/PIPELINE_IO.md`；训练链参考 `ChemDB_restructured/ChemDB/training/` |
| Phase 1A 实施记录 | [phase_1a_implementation.md](./phase_1a_implementation.md) |

---

## 0. 执行摘要

| 结论 | 详情 |
|------|------|
| 是否已是可 import 包 | **否**（无 `pyproject.toml`，`src/` 无 `__init__.py`，平铺 `from utils import`） |
| Run 路径隔离成熟度 | **低**：多数 step 写死 `./tmp`、`./data`；仅 step4_5/12/13、`build_L3_embedding`、部分 training 脚本有目录参数 |
| 仅靠 `subprocess(cwd=run_root)` | **仅部分可行**；step1/2/6_7/9/11、zscore、`main.py`、服务 API **必须改参数或包一层** |
| 第一阶段最小改造 | 引入 **RunContext + 编排层**；为 step1/2/6_7/8/9/11 补目录参数；**每 Run 生成 `index.json`/`config.yaml`**；pytest 断言产物不写回 repo `tmp/` |

**仓库缺口**：`ChemDB`（原 `ChemDB_run_based`）**不含 `training/` 目录**；训练/评估审计以 `ChemDB_restructured/ChemDB/training/` 为准，合并进目标 `ChemDB/` 时一并纳入改造。Phase 1A 已引入 `app/services` 编排层，见实施记录。

---

## 1. 路径读写扫描

### 1.1 `src/step*.py` — 模块级常量（相对进程 cwd）

| 文件 | 行号区间 | 硬编码路径 | 读写 |
|------|----------|------------|------|
| `src/step1.py` | L9–13, L155–180, L220–224 | `INPUT_DATA_DIR=./data/pubchem`, `METAL_LIST_PATH=./data/metal_list.txt`, `OUTPUT_DIR=./tmp` | 读 pubchem/metal；写 tmp CSV/JSON |
| `src/step2.py` | L9–14, L123–146, L154 | `INPUT_DIR=./tmp`, `OUTPUT_DIR=./tmp` | 读 `ligand_data.csv`；写 `repaired_ligand_data.csv` |
| `src/step4_5.py` | L23–24, L507–526 | 构造器 `input_dir`/`output_dir` 默认同上 | 读 `repaired_ligand_data`；写 GA/IRL/GAC 等 |
| `src/step6_7.py` | L7–14, L229–233, L277–306, L338 | `GA_PATH=./tmp/GA_with_id.csv`, `IRL_PATH=./tmp/IRL_filtered.csv`, `INPUT_DIR=./tmp` | 读 GA/IRL/repaired；写 `fragments.csv` |
| `src/step8.py` | L9–12, L50–85, L207–211 | `INPUT_DIR=./tmp`, `METAL_LIST=./data/metal_list.txt` | 读 tmp 多表；**L211** `irl_file = Path("tmp/IRL_filtered.csv")`（非 `INPUT_DIR`） |
| `src/step9.py` | L10–12, L104–145 | `INPUT_DIR=./tmp`, `DATA_DIR=./data` | 读 tmp 四表 + **`data/IRL_filtered.csv`** |
| `src/step10.py` | L41–43, L118–155 | `INPUT_DIR=./tmp`, `DATA_DIR=./data` | 读 tmp + `data/IRL_filtered`（未进 `main.py`） |
| `src/step11.py` | L12–13, L84–107 | `INPUT_DIR`/`OUTPUT_DIR=./tmp` | 读 repaired + neo4j 关系；写 embeddings |
| `src/step12.py` | L18–19, L45–73 | 默认同上；**CLI** L530–557 可覆盖 | 读 neo4j 关系 + gac；写 JSON/CSV |
| `src/step_13_complete.py` | L13–14, L815–879, L1418–1419 | 默认同上；**CLI** L1364–1457；**fallback** `parent.parent/result_analysis/...` | 读 Step12 产物；写 step13 CSV |

### 1.2 `src/` 其他

| 文件 | 位置 | 路径行为 |
|------|------|----------|
| `src/main.py` | L22–23, L33–96 | `sys.path.insert(0, _script_dir/"src")`；`STEP_DEFINITIONS` 中 `./tmp/...` 状态检测 |
| `src/server.py` | L20–40, L188–192 | `DATASETS[*].csv_path = "./tmp/..."`；stats 列表 `tmp/*.json` |
| `src/affinity_api.py` | L58–59, L373–374 | `BASE_DIR = Path(__file__).parent.parent`；`INPUT_DIR = BASE_DIR/"tmp"`（相对源码树，不随 cwd） |
| `src/molclr_api.py` | L42, L107–123, L160–166 | `base_dir = dirname(dirname(__file__))` → ChemDB 根；MolCLR 克隆到 `{root}/tmp/MolCLR` |
| `src/comparison.py` | L29–32, L284–325 | `PROJECT_ROOT = abspath(join(dirname(__file__), pardir))`；写 `PROJECT_ROOT/tmp/comparison` |
| `src/proxy.py` | — | 无本地数据路径（反代） |

### 1.3 `src/tools/*.py`

| 文件 | 路径逻辑 |
|------|----------|
| `build_L3_embedding_index.py` | L26–28 `REPO_ROOT=parent(src)`；L50–51 默认 `{REPO_ROOT}/tmp`、`{REPO_ROOT}/data/L3_embedding`；**CLI** `--tmp-dir`, `--out-dir` |
| `DID_embedding.py` | L22–28 `PROJECT_ROOT=dirname(dirname(file))` → **`src/`（错误层级）**；`DEFAULT_DID_CSV` → **`src/tmp/repaired_ligand_data.csv`** |
| `L4_create_independent.py` | L11 硬编码 `/home/chenghua/ZhuoLi/ChemDB_restructured/ChemDB/src` |
| `CC_split.py` | L722–730 `__main__` 内 `/home/weilin/...` 示例路径 |
| `distance.py` | `load_distance_calculator(data_dir)` 参数化 |

### 1.4 `training/`（`ChemDB_restructured/ChemDB`，待迁入 `ChemDB/`）

| 文件 | 路径逻辑 |
|------|----------|
| `training/index.json` | **绝对路径** 指向 `ChemDB_restructured/.../tmp`、`data/L3_embedding` |
| `training/data.py` | `ensure_split_index(training_dir, ...)`；`load_index_paths(index_path)`；**CLI** `--training-dir`, `--index` |
| `training/train.py` | `resolve_config`：`training_dir`、`index_path`、`ckpt_dir` 默认 `training/`、`training/ckpts` |
| `training/config.yaml` | `training_dir: null` → 解析为 `training/` 包目录 |
| `training/evaluation/run_all_quests.py` | 默认 `training/models`、`evaluation/quests`、`evaluation/results`；**CLI** 均可覆盖 |
| `training/evaluation/core.py` | L17–24 `sys.path` 插入 `CHEMDB_ROOT`、`SRC_DIR` |
| `training/evaluation/legacy/*` | 大量 `TRAINING_ROOT`、`tmp_aromatic_small_rings` 等实验路径 |
| `data/metal_embedding/zscore_element_features.py` | `ROOT=Path(__file__).parent`；固定同目录 in/out |

### 1.5 `tools/` 与根目录

| 文件 | 路径 |
|------|------|
| `tools/neo4j/import.sh` | L7 `-v "../../tmp":/import` |
| `tools/chemtool.py` | `sys.path.insert(SRC_PATH)` |
| `test_mingxin_Ni_polymer/*.py` | 默认 `root/tmp`；输出 JSON 含旧绝对路径 `ChemDB_restructured` |

---

## 2. 各脚本目录参数支持矩阵

图例：**CLI** = 命令行；**Ctor** = Processor 构造器；**cwd** = 仅当 `chdir(run_root)` 且 run 内有 `data/`、`tmp/` 子目录时有效。

| 步骤 ID | 脚本 | CLI 目录参数 | Ctor 目录参数 | cwd=run_root alone |
|---------|------|--------------|---------------|-------------------|
| step1 | `step1.py` | ❌（仅 workers/limit） | ❌ | ❌ 需 `data/pubchem` |
| step2 | `step2.py` | ❌ | ❌ | ⚠️ 需已有 `tmp/ligand_data` |
| step4_5 | `step4_5.py` | ✅ `-i/-o`, workers, limit | ✅ | ✅ 推荐显式传参 |
| step6_7 | `step6_7.py` | ❌ | ❌ | ⚠️ 依赖 `GA_PATH`/`IRL_PATH`=`./tmp/...` |
| step8 | `step8.py` | ❌ | ❌ | ⚠️ + `data/metal_list.txt` |
| step9 | `step9.py` | ❌ | ❌ | ❌ 还需 `data/IRL_filtered.csv` |
| step11 | `step11.py` | ❌ | ❌ | ⚠️ 全在 tmp 时可工作 |
| step12 | `step12.py` | ✅ `--input-dir`, `--output-dir`, `-K` | CLI 覆盖实例 | ✅ |
| step13 | `step_13_complete.py` | ✅ input/output, `-K`, `-D`, metal index | CLI 覆盖 | ✅（注意 fallback） |
| build_L3 | `build_L3_embedding_index.py` | ✅ `--tmp-dir`, `--out-dir` | — | ✅ out 应指向 run/data |
| zscore_metal | `zscore_element_features.py` | ❌ | ❌ | ❌ |
| training.data | `training/data.py` | ✅ `--training-dir`, `--index` | ✅ | ❌ 依赖 index 路径 |
| training.train | `training/train.py` | ✅ `--config` | ✅ | ❌ |
| evaluation | `run_all_quests.py` | ✅ `--results-dir`, `--index`, … | — | ❌ |

---

## 3. cwd / sys.path / 硬编码依赖

### 3.1 相对路径（隐式 cwd）

所有 `OUTPUT_DIR = "./tmp"` 在未 `chdir` 时解析为**启动进程时的 cwd**。  
`python src/main.py`（仓库根）→ `tmp/` 在根下；`cd src && python step1.py` → `src/tmp/`（已知陷阱）。

### 3.2 `sys.path.insert`

| 文件 | 片段 |
|------|------|
| `src/main.py:23` | `sys.path.insert(0, str(_script_dir / "src"))` — `_script_dir` 为 `src/`，实际可能插入 **`src/src`** |
| `src/molclr_api.py:42,168` | 插入 ChemDB 根、`tmp/MolCLR` |
| `src/tools/build_L3_embedding_index.py:29-30` | 插入 `SRC_DIR` |
| `src/tools/DID_embedding.py:23-24` | 插入 **`src/`**（非仓库根） |
| `src/tools/L4_create_independent.py:11` | 绝对路径旧仓库 |
| `tools/chemtool.py:27` | 插入 `SRC_PATH` |
| `test_mingxin_Ni_polymer/*` | 插入 `SRC` / `CHEMDB_ROOT` |
| `training/evaluation/core.py:21-24` | 插入 `CHEMDB_ROOT`、`SRC` |

**未发现** `os.getcwd()` 调用。

### 3.3 硬编码绝对路径

| 文件 | 内容 |
|------|------|
| `training/index.json` | 全部 path 指向 `/home/chenghua/ZhuoLi/ChemDB_restructured/...` |
| `src/tools/L4_create_independent.py:11` | `/home/chenghua/ZhuoLi/ChemDB_restructured/ChemDB/src` |
| `step_13_complete.py:1419` | `parent.parent/result_analysis/analysis_candT/metal_l3_index.csv` |
| `test_mingxin_Ni_polymer/outputs/*.json` | 旧 ChemDB_restructured ckpt/tmp 路径 |
| `CC_split.py:730` | 仅 `__main__` 演示 |

### 3.4 与 Run 布局的特殊冲突

| 问题 | 说明 |
|------|------|
| Run `data/` vs 源码 `data/` | step1/9 读 `./data`；编排需 **cwd=run_root** 或显式传参 |
| 全局 `data/L3_embedding`（约 12G） | `build_L3` 默认写仓库；Run 须 `--out-dir run/data/L3_embedding` |
| `affinity_api` / `server` | 绑定源码树根 `tmp/`，不随 Run 切换 |
| `molclr_api` | 在仓库 `tmp/MolCLR` 克隆；多 Run 争用 |

---

## 4. 最小 RunContext 设计（第一阶段）

```text
RunContext
├── run_root          # workspace/projects/{pid}/runs/{rid}/
├── data_dir          # run_root/data/     (pubchem, metal_list, L3_embedding, metal_embedding)
├── tmp_dir           # run_root/tmp/
├── training_dir      # run_root/training/ (index.json, pkl, ckpts, evaluation/results)
├── log_dir           # run_root/logs/
├── report_dir        # run_root/reports/
└── manifest_dir      # run_root/manifests/
```

**派生规则：**

```text
data_dir     = run_root / "data"
tmp_dir      = run_root / "tmp"
training_dir = run_root / "training"
log_dir      = run_root / "logs"
report_dir   = run_root / "reports"
manifest_dir = run_root / "manifests"
```

**编排层职责（不接 Streamlit/SQLite）：**

1. 创建目录树；复制或 symlink 共享资源（`metal_list.txt`、`p_elements_list.txt`、`element_features.csv`）到 `data_dir`。
2. 用户 PubChem CSV → `data_dir/pubchem/`。
3. 生成 `training_dir/index.json`（路径均在 run 内）。
4. 生成 `training_dir/config.yaml`（`training_dir`、`ckpt_dir` 指向 run）。
5. 调用 step：`--input-dir` / `--output-dir` 或统一环境变量 `CHEMDB_RUN_ROOT` 等。

**流水线是否从 Run 读参？**  
第一阶段：**是，经编排层间接读**（RunContext → CLI/env → 现有脚本）；**不要求** step 内查数据库。

---

## 5. Step Registry 草案

```text
REGISTRY: step_id | script | cwd | command | requires[] | produces[] | checks[]
```

| step_id | script | 建议 cwd | 命令要点 | requires | produces | 检查规则 |
|---------|--------|----------|----------|----------|----------|----------|
| step1 | `src/step1.py` | ChemDB 或 run_root* | 需 `--data-dir --tmp-dir` | `data/pubchem/*.csv`, `data/metal_list.txt` | `tmp/complex_data.csv`, `tmp/ligand_data.csv` | 行数>0 |
| step2 | `src/step2.py` | 同上 | 需 `--input-dir --output-dir` | `tmp/ligand_data.csv` | `tmp/repaired_ligand_data.csv` | 含 `ligand_new_did` |
| step4_5 | `src/step4_5.py` | ChemDB | `step4_5.py -i {tmp} -o {tmp}` | `tmp/repaired_ligand_data.csv` | `GA_with_id.csv`, `ligand_with_gac.csv`, `IRL_filtered.csv` | GA 行数>0 |
| step6_7 | `src/step6_7.py` | **run_root**‡ | 需 Ctor/env 覆盖 GA/IRL 路径 | repaired, GA, IRL | `tmp/fragments.csv` | 非空 |
| step8 | `src/step8.py` | **run_root**‡ | 修 L211 为 `{tmp}/IRL_filtered` | step1–7 产出 + metal_list | `neo4j_*.csv` | l3_l4 关系存在 |
| step9 | `src/step9.py` | **run_root**‡ | 需 `data_dir`；IRL 建议用 tmp | 多 tmp CSV | `tmp/did_index.csv` | 可选 |
| step11 | `src/step11.py` | **run_root**‡ | 需目录 CLI | step8 关系 | `l3_embeddings.csv` | 可选 |
| step12 | `src/step12.py` | ChemDB/src | `--input-dir --output-dir --K` | neo4j 关系 + gac | `l3_l5.json`, `m_l3_pairs.csv`, … | JSON 可解析 |
| step13 | `step_13_complete.py` | ChemDB/src | `--kl-nl-only` + in/out dir | Step12 全套 | `step13_kl_nl_samples.csv` | 列 T,M,label |
| build_L3_embedding | `tools/build_L3_embedding_index.py` | ChemDB | `--tmp-dir --out-dir {data}/L3_embedding` | metal_l3_index, repaired | `L3_embedding_*.npz` | npz keys 正确 |
| zscore_metal_embedding | `zscore_element_features.py` | metal_embedding | 需 `--in --out` | `element_features.csv` | `element_features_zscore.csv` | 行数一致 |
| training.data | `python -m training.data` | ChemDB | `--training-dir --index --seed` | index→step13 + m_l3_pairs | `split_index.json`, `*_records.pkl` | 按 run 隔离 |
| training.train | `python -m training.train` | ChemDB | `--config {training}/config.yaml` | pkl + embeddings | `best.pt`, `history.json` | history 有指标 |
| evaluation | `run_all_quests.py` | ChemDB | `--results-dir --index --models-dir` | models, quests | `summary_all_tasks.csv` | manifest 存在 |

\* 理想：`python -m chemdb.steps.step1 --run-root {run_root}`。  
‡ `cwd=run_root` 且存在 `run_root/data`、`run_root/tmp` 时，`./tmp`、`./data` 才一致。

**依赖 DAG：**

```text
step1 → step2 → step4_5 → step6_7 → step8 → step12 → step13
                              ↓
                    build_L3_embedding (常需 step13 的 metal_l3_index)
step13 + m_l3_pairs → training.data → training.train → evaluation
zscore_metal_embedding → training (metal path in index)
```

---

## 6. subprocess 兼容 vs 必须改参

### 6.1 可暂时 subprocess + 现有 CLI

| 脚本 | 条件 |
|------|------|
| `step4_5.py` | `-i/-o` = `{tmp_dir}`；`PYTHONPATH=ChemDB/src` |
| `step12.py` / `step_13_complete.py` | 同上；禁用/重写 `result_analysis` fallback |
| `build_L3_embedding_index.py` | `--tmp-dir`、`--out-dir` 指向 run |
| `training.data` / `train` / `run_all_quests` | `--training-dir`、`--index`、`--config`、`--results-dir`；**须预生成 index** |

### 6.2 必须改代码或薄包装

| 脚本 | 原因 | 最小改法 |
|------|------|----------|
| `step1.py` | 无 tmp/data 参数 | CLI：`--pubchem-dir`, `--tmp-dir`, `--metal-list` |
| `step2.py` | 无 in/out | `--input-dir`, `--output-dir` |
| `step6_7.py` | `GA_PATH`/`IRL_PATH` 模块常量 | Ctor 或 `CHEMDB_TMP_DIR` |
| `step8.py` | 无 CLI；L211 裸 `tmp/` | 统一 `INPUT_DIR`；CLI |
| `step9.py` | `DATA_DIR` 分离 | `--data-dir`；IRL 改读 tmp |
| `step11.py` | 无目录 CLI | 同 step8 |
| `zscore_element_features.py` | `ROOT=__file__.parent` | `--input`/`--output` |
| `main.py` | 写死 `./tmp` | Phase 1 不用；用 Registry |
| `server` / `affinity_api` / `molclr_api` | 绑定仓库 tmp | Phase 2：`--run-root` |
| `DID_embedding.py` | PROJECT_ROOT→`src/` | 改为 REPO_ROOT |
| `L4_create_independent.py` | 绝对 path | 包内 import |

### 6.3 建议子进程环境模板

```text
env:
  CHEMDB_REPO_ROOT = /path/to/ChemDBWebVersion/ChemDB
  CHEMDB_RUN_ROOT   = /path/to/workspace/.../run_001
  PYTHONPATH        = {CHEMDB_REPO_ROOT}/src
cwd:
  优先 run_root（step6_7–11、step8 修完后）
  或 ChemDB/src（step12/13 现状）
```

---

## 7. pytest 方案（不写回源码目录）

**目录**：`ChemDB/tests/run_isolation/`（目标布局下新建）

| 测试 ID | 目的 |
|---------|------|
| `test_run_workspace_fixture` | `tmp_path/run_001/{data,tmp,training,logs,manifests}` |
| `test_run_context_paths` | 断言均在 `tmp_path` 下 |
| `test_step4_5_writes_only_under_tmp_dir` | 小 fixture + step4_5；repo `tmp/` 无新文件 |
| `test_training_data_respects_training_dir` | pkl 仅在 run `training/` |
| `test_index_json_paths_inside_run` | index 路径前缀 ⊂ run_root |
| `test_no_write_to_src_tmp` | 回归 `DID_embedding` 误写 `src/tmp/` |
| `test_subprocess_cwd_run_root` | 编排层 cwd 断言 |

**要点：**

- `conftest.py`：`repo_root`、`run_root`、`repo_root/tmp/.pytest_guard`。
- 前后对比 repo `tmp/` 文件集合相等（白名单 guard）。
- 大 step 用 `@pytest.mark.integration` 可选。
- CI 强制 `CHEMDB_RUN_ROOT=tmp_path`。

---

## 8. 分阶段改造计划

### Phase 0 — 仓库整理

- 将 `ChemDB_restructured/ChemDB/training` 迁入 `ChemDBWebVersion/ChemDB/training`。
- 建立 `workspace/`、`app/services/`、`app/templates/` 骨架（无 Streamlit）。
- 按需清理 §10 冗余数据。

### Phase 1 — Run 路径隔离（当前目标）

| 任务 | 产出 |
|------|------|
| P1.1 | `app/services/run_context.py` |
| P1.2 | `app/services/step_registry.py` |
| P1.3 | `app/services/orchestrator.py` |
| P1.4 | `app/templates/run_index.json.j2`, `run_config.yaml.j2` |
| P1.5 | ChemDB 补丁：step1/2/6_7/8/9/11、step8 L211、zscore、index 生成器 |
| P1.6 | `tests/run_isolation/` |
| P1.7 | 更新 `doc/PIPELINE_IO.md` Run 物理路径一节 |

**不做**：Streamlit、SQLite、pip 包（Phase 2+）。

### Phase 2 — 元数据与 UI

- `workspace/chemdb.sqlite` + Project/Run/Artifact。
- Streamlit → orchestrator。
- API 按 `run_id` 选数据。

### Phase 3 — 可安装包与多 Run 服务

- `pyproject.toml`、`import chemdb.*`。
- MolCLR 缓存 per-run 或共享只读。

---

## 9. 目标目录映射

```text
ChemDBWebVersion/
├── ChemDB/                    # 原 ChemDB_run_based（已重命名）+ training 待合并
│   ├── src/
│   ├── training/
│   ├── data/                  # 模板/默认资源，非 Run 产物
│   └── docs/
├── app/
│   ├── services/              # RunContext, Registry, Orchestrator
│   ├── templates/
│   └── streamlit_ui/          # Phase 2
└── workspace/
    ├── chemdb.sqlite          # Phase 2
    └── projects/.../runs/run_001/
        ├── data/
        ├── tmp/
        ├── training/
        ├── logs/
        ├── reports/
        └── manifests/
```

**原则**：源码树 `ChemDB/data/L3_embedding` 视为历史全局缓存；新 Run 的 npz **只写入** `run_001/data/L3_embedding/`。

---

## 10. 建议清除的文件（ChemDB，扫描时目录名为 ChemDB_run_based）

删除前请确认无未备份实验。

### 10.1 大型实验副本

| 路径 | 约大小 | 说明 |
|------|--------|------|
| `tmp/tmp_aromatic_small_rings/` | **35G** | 完整流水线重复 + pydentate |
| `tmp/tmp_adamantyl/` | **12G** | 另一套实验 tmp |
| `data/L3_embedding/*.npz` | **12G** | 应迁至 Run 或外部存储 |

### 10.2 根 `tmp/` 与工具缓存

| 路径 | 说明 |
|------|------|
| `tmp/fragments.csv`（约 3.5G）等 | 可与子目录去重 |
| `tmp/MolCLR/`、`tmp/ckpt/` | 迁缓存目录或按需下载 |
| `src/tmp/` | 误跑污染（若存在） |

### 10.3 日志与无关文档

| 路径 | 说明 |
|------|------|
| `step13_full_20260404_103649.log` | 一次性日志 |
| `t-SMILES.pdf` | 非代码资产 |
| `tmp/**/pydentate/` | 独立分析，非核心流水线 |
| `tmp/**/*.log` | 运行日志 |

### 10.4 Neo4j 运行时

| 路径 | 说明 |
|------|------|
| `tools/neo4j/neo4j/data/`、`logs/` | 应在 `.gitignore`，不进 git |

### 10.5 测试输出

| 路径 | 说明 |
|------|------|
| `test_mingxin_Ni_polymer/outputs/*.json` | 含过期绝对路径 |
| `test_mingxin_Ni_polymer/__pycache__/` | 可忽略 |

### 10.6 缺失项

| 项 | 说明 |
|----|------|
| `training/` | run_based **缺失**；从 restructured 复制 |
| `app/`、`workspace/` | 目标新目录，尚未存在 |

### 10.7 建议保留

| 路径 | 说明 |
|------|------|
| `data/pubchem/` | 默认样本或迁至 `workspace/shared_inputs/` |
| `data/metal_list.txt`、`p_elements_list.txt` | 模板输入 |
| `doc/PIPELINE_IO.md` | 逻辑 IO |
| `src/`、`tools/chemtool.py`、`requirements.txt` | 核心代码 |

---

## 11. 第一阶段验收标准

1. 在 `workspace/.../run_001/` 完成 step1→8→12→13（小 `record_limit`）后，**repo 根 `tmp/` mtime 不变**（pytest 守护）。
2. 两个 run 目录互不覆盖 `step13_kl_nl_samples.csv`。
3. `training_dir/index.json` 路径均落在对应 run 下。
4. `manifests/` 记录每步命令、退出码、产物列表。
5. 无需 Streamlit/SQLite：`python -m app.services.orchestrator run --run-root ... --through step13`。

---

## 修订记录

| 版本 | 日期 | 说明 |
|------|------|------|
| 1.0 | 2026-05-20 | 初始评估，基于 ChemDB_run_based 源码扫描 |
| 1.1 | 2026-05-20 | 补充 Phase 1A 状态与路径更名说明；实施见 phase_1a_implementation.md |
