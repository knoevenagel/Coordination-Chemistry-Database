# Phase 1E 审计与实施计划：Workspace-level GA Sets + Run GA Binding

**日期**：2026-05-21  
**状态**：计划（本文件仅审计与设计，**不实现代码**）  
**前置**：Phase 1A（run orchestrator + core pipeline）、1C（training）、1D（model-after）已完成。

---

## 目标摘要

| 原则 | 说明 |
|------|------|
| GA 资源 workspace 级 | 可复用、可编辑、可版本化；目录 `workspace/ga_sets/` |
| Run 绑定 GA version | 每个 run 显式绑定某一 `ga_set_id` + `ga_version_id` |
| 物化快照 | 绑定/应用时 **复制** GA CSV 到 `run/tmp/GA_with_id.csv`，不引用可变 workspace 文件 |
| 下游不变 | step6_7 及之后仍读 `run/tmp/GA_with_id.csv` |
| Stale 语义 | GA version 变更 → GAC/IRL 及全部下游产物视为过期，需重跑；生成报告不自动删文件 |
| 本阶段边界 | 仅 file-based registry；**无 SQLite、无 Streamlit**；不改 training/model_after 化学/训练逻辑 |

---

## 一、审计 `ChemDB/src/step4_5.py`

### 1.1 是否耦合在同一次运行中？

**是。** `Step4_5Processor.process_all()` 在单次调用中顺序执行：

| 内部步骤 | 方法 | 产物 |
|----------|------|------|
| Step 0 | `step0_deduplicate()` | `ligand_data_deduplicated.csv`（可选中间） |
| Step 1 | `step1_extract_GA()` | `ga_patterns_temp.csv`（临时） |
| Step 2 | `step2_generate_GA_ID()` | **`GA_with_id.csv`** |
| Step 3 | `step3_calculate_GAC()` | **`ligand_with_gac.csv`** |
| Step 4 | `step4_extract_IRL()` | **`IRL_filtered.csv`** |
| Step 5 | `step5_deduplicate_IRL()` | `IRL_filtered_cleaned.csv`（可选） |

`main()` 仅暴露 `--input-dir` / `--output-dir` / `--workers` / `--limit`，**无 mode 开关**，默认 `process_all()` = 全流程。

GA extraction、GA_ID generation、GAC calculation、IRL filtering **逻辑上可分、实现上绑死在一次 `process_all()`**。

### 1.2 函数职责映射

| 能力 | 函数 / 类 | 输入 | 输出 |
|------|-----------|------|------|
| 从 `repaired_ligand_data.csv` 提取 GA 候选 | `Step4_5Processor._get_smiles_DID_from_csv()` | CSV 行 | `(DID, SMILES)` 列表 |
| 单分子 GA 提取 | `process_single_molecule_for_GA()` + `FindGA` | 配体 SMILES | `GA_dict: ring_smiles → {DID}` |
| 合并全库 GA | `step1_extract_GA()` | 并行 `process_single_molecule_for_GA` | `ga_patterns_temp.csv`（`GA_SMILES`, `GA_source`） |
| 生成 GA_ID | `generate_ga_id()`；`step2_generate_GA_ID()` | `GA_SMILES` | `GA_ID` 字符串 |
| 写 `GA_with_id.csv` | `step2_generate_GA_ID()` | temp CSV | **`GA_with_id.csv`**（仅 `GA_SMILES`, `GA_ID`） |
| 计算 GAC | `calculate_gac()`；`_process_gac_worker()`；`step3_calculate_GAC()` | 配体 + `ga_list`（来自 `GA_with_id` 的 SMILES 列） | **`ligand_with_gac.csv`**（`DID`, `SMILES`, `GAC`） |
| 生成 IRL | `ExtractIRLObject`；`step4_extract_IRL()` | `ligand_with_gac.csv` | **`IRL_filtered.csv`**（`DID`, `SMILES`, `GAC`） |
| IRL 去重 | `step5_deduplicate_IRL()` | `IRL_filtered.csv` | `IRL_filtered_cleaned.csv` |

**注意**：`step3` 只使用 `GA_with_id.csv` 的 **`GA_SMILES` 列表** 做子结构匹配计数，**不读 `GA_ID`**。`GA_ID` 主要为 step6_7 原子标注服务。

### 1.3 最小改动：`--mode` 拆分方案

**可行**，且无需改化学算法，仅改 **调度入口**（`main()` + 薄包装方法）。

建议新增 CLI：

```text
--mode {full-auto,generate-ga,apply-ga}   # 默认 full-auto
```

| mode | 调用链 | 输入要求 | 输出 |
|------|--------|----------|------|
| **generate-ga** | `step0`（可选，建议保留）→ `step1` → `step2` | `repaired_ligand_data.csv` | `GA_with_id.csv`（+ 可选 `ga_patterns_temp` 删除逻辑不变） |
| **apply-ga** | `step0`（可选）→ `step3` → `step4` → `step5` | `repaired_ligand_data.csv` + **已存在** `{output_dir}/GA_with_id.csv` | `ligand_with_gac.csv`, `IRL_filtered.csv`, … |
| **full-auto** | 现有 `process_all()` | 同现网 | 同现网 |

实现要点（最小 diff）：

1. 在 `Step4_5Processor` 增加 `process_generate_ga()`、`process_apply_ga()`，分别从 `process_all()` 抽出对应步骤子集。
2. `main()` 根据 `--mode` 分支调用；`full-auto` 仍调用 `process_all()`，保证 orchestrator 现有 `step4_5` 子进程 **零行为变化**（默认参数下）。
3. `apply-ga` 在 `step3` 前 **断言** `self.ga_output.exists()`，若缺失则失败并明确报错（不隐式 generate）。
4. orchestrator `step_registry` 中 `step4_5` 的 `command_argv` 可增加 `--mode full-auto`（显式文档化，可选）。

**不建议** 在 Phase 1E 拆成两个独立脚本文件；单文件 + mode 更易测试与复用 `generate_ga_id` / `calculate_gac`。

### 1.4 mode 语义核对（与设计要求对齐）

| 要求 | 方案 |
|------|------|
| generate-ga 只生成候选 GA，不算 GAC/IRL | ✅ 仅 step0–2 |
| apply-ga 读已有 `GA_with_id.csv` 再算 GAC/IRL | ✅ 仅 step3–5，GA 文件只读 |
| full-auto 保持 pipeline 兼容 | ✅ 默认 `process_all()` |

**Phase 1E 与 mode 的关系**：workspace 操作 `create-ga-set-from-run` 调用 `generate-ga`（输出写到 ga version 目录）；`apply-ga-to-run` 调用 `apply-ga`（读写均在 run `tmp/`）。

---

## 二、审计所有 GA consumers

### 2.1 搜索汇总表

| file_path | function/class | read/write | current path assumption | supports run_root? | Phase 1E impact |
|-----------|----------------|------------|-------------------------|-------------------|-----------------|
| `ChemDB/src/step4_5.py` | `Step4_5Processor` | **W** GA_with_id；R repaired | `{output_dir}/GA_with_id.csv`，默认 `./tmp` | ✅ `--input-dir`/`--output-dir` | 增加 `--mode`；generate/apply 拆分；orchestrator 传 `{tmp_dir}` |
| `ChemDB/src/step6_7.py` | `Step6_7Processor._load_ga_data` | **R** GA_with_id | `{input_dir}/GA_with_id.csv`（模块常量 `GA_PATH` 仅文档遗留） | ✅ ctor `input_dir` | **无改**；继续读 run 物化文件 |
| `ChemDB/src/step6_7.py` | `MoleculeMarker` | R（内存 dict） | 来自 GA CSV 的 `GA_SMILES→GA_ID` | 间接 | 无改 |
| `ChemDB/src/tools/L4_create.py` | `L4_Create` | R `ga_list` SMILES | 调用方传入 list，非直接读文件 | N/A | 无改（由 step6_7 喂数据） |
| `ChemDB/src/tools/L4_create_independent.py` | 同上 | R | 同上 | N/A | 无改 |
| `ChemDB/src/step8.py` | Step8 processor | R `IRL_filtered`, `fragments` | `{input_dir}/...` | ✅ | 间接依赖 GA（经 IRL/fragments） |
| `ChemDB/src/step12.py` | `extract_l3_gac` | R `ligand_with_gac.csv` | `{input_dir}/ligand_with_gac.csv` | ✅ | 间接 stale |
| `ChemDB/src/step_13_complete.py` | KL/NL 采样 | R `l3_gac.json` 等 | `{input_dir}/` | ✅ orchestrator | 间接 stale |
| `ChemDB/src/step9.py` / `step10.py` | legacy 图构建 | R `IRL_filtered` | `./tmp` 或 `data_dir` | ⚠️ 非 Web orchestrator 主路径 | 1E 不纳入 |
| `ChemDB/src/main.py` | pipeline 描述 | 文档 | `./tmp` | ❌ | 无改 |
| `ChemDB/src/server.py` | 服务配置 | R IRL path | `./tmp/IRL_filtered.csv` | ❌ | 1E 不做 Streamlit/server |
| `ChemDB/doc/PIPELINE_IO.md` | 文档 | — | `./tmp` | — | 实施后更新文档 |
| `app/services/step_registry.py` | `step4_5`, `step6_7` | 契约 | `tmp/GA_with_id.csv` 等相对 run_root | ✅ | 可选新增 logical steps；`step4_5` 仍可 full-auto |
| `app/services/orchestrator.py` | `run_step` | 调度 | 通过 `{tmp_dir}` 占位符 | ✅ | 不塞 GA 逻辑；可选调用 `ga_registry_manager` |
| `tests/run_isolation/test_integration_1pct.py` | 集成断言 | R | run `tmp/` | ✅ | 后续可加 1E 测试 |
| `ChemDB/training/**` | — | — | — | — | **无 GA 直接依赖** |
| `ChemDB/training/model_after/**` | — | — | — | — | **无 GA 直接依赖** |
| `workspace/.../manifests/step4_5.json` | 运行记录 | — | 绝对路径 | ✅ | 绑定后 step4_5 manifest 语义变（apply vs full-auto） |

### 2.2 重点结论

**直接读 `run/tmp/GA_with_id.csv`**

- `step6_7.py`（必需；缺失时 warning 并 `ga_data={}`，pipeline 会退化）

**写 `GA_with_id.csv`**

- 仅 `step4_5.py` `step2_generate_GA_ID()`（Phase 1E 增加：**ga_registry bind** 复制写入，不算“生成”）

**依赖 GA 列表的 downstream artifacts（间接链）**

```text
GA_with_id.csv
  → ligand_with_gac.csv          (step4_5 step3)
  → IRL_filtered.csv             (step4_5 step4; step6_7 也读)
  → fragments.csv                (step6_7; 依赖 GA+IRL 标注)
  → neo4j_*.csv                  (step8)
  → l3_l5.json, l5_l3.json, l3_gac.json, m_l3_pairs.csv, … (step12)
  → step13_kl_nl_samples.csv     (step13)
  → data/L3_embedding/L3_embedding_ecfp.npz (training)
  → training/index.json, config.yaml, *.pkl, ckpts/*.pt
  → workspace/model_after_results/** (评估结果；绑定/重训后应视为过期)
```

**training / model_after**

- 代码库内 **无** `GA_with_id` / `GA_SMILES` 字符串引用。
- 依赖链：**GA → GAC → IRL → fragments → 图/样本 → embedding → 训练 → checkpoint → model_after**。
- Phase 1E **不修改** `training/`、`model_after/` 逻辑；通过 stale 报告提示重跑 pipeline。

---

## 三、审计 GA schema 与 ID 生成

### 3.1 当前真实字段（`run_demo` 实测）

文件：`workspace/projects/p001/runs/run_demo/tmp/GA_with_id.csv`

- 列：**`GA_SMILES`, `GA_ID`**（仅此两列写入最终 CSV）
- 行数：126（无表头重复）
- 中间文件 `ga_patterns_temp.csv` 含 `GA_SMILES`, `GA_source`（**不**进入最终 downstream 文件）

### 3.2 下游最低列需求

| 消费者 | 必需列 | 备注 |
|--------|--------|------|
| step4_5 `step3` | `GA_SMILES`（从 CSV 读入为 list） | 不读 `GA_ID` |
| step6_7 | `GA_SMILES` + `GA_ID`（或兼容 `SMILES`/`ID`） | 建 dict 时 **同 SMILES 后者覆盖前者** |
| Phase 1E materialize | 至少上述两列 | 额外列可在 workspace 版本保留，**物化到 run 时建议只写两列**（或保留额外列但 step6_7 忽略） |

**结论**：物化到 `run/tmp/GA_with_id.csv` 时 **必须保留 `GA_SMILES`, `GA_ID`**；workspace 版本 CSV 可含 `source`, `active`, `note` 等，但 **bind 物化** 时应提供兼容视图（默认 strip 为两列，或在文档中规定 step6_7 已忽略额外列）。

### 3.3 `generate_ga_id` 是否 deterministic？

**是（对给定 SMILES 字符串）。** 实现：`MD5(smiles.encode())` → 取 hex 前 12 位 → 转 int → `G{value % 1000000:06d}`。

- 同一 `GA_SMILES` 字节串 → 同一 `GA_ID`。
- **未** 对输入再做 canonicalize；依赖上游 `MolToSmiles(..., canonical=True)` 已规范化 ring SMILES。

### 3.4 GA_ID collision 风险

- 空间：约 \(10^6\) 个 ID（`G000000`–`G999999`）。
- 不同 SMILES 可能映射到同一 ID（MD5 截断 + mod）；**存在理论碰撞**。
- 当前 `run_demo`：**0** 个重复 `GA_ID`、**0** 个重复 `GA_SMILES`。
- Phase 1E `create-ga-version-from-csv` 应：**导入时检测 `GA_ID` 唯一**；若新生成 ID 与已有冲突则报错或递增重试（实施时二选一，建议 **检测 + fail** 并提示人工处理）。

### 3.5 GA_SMILES 是否 canonicalized？

**在提取阶段是。**

- 芳香环 / 三元环 / 四元环：`Chem.MolToSmiles(submol, canonical=True)`（`FindGA` / `process_single_molecule_for_GA`）。
- workspace 手工编辑 CSV 时：**必须在 `ga_registry_manager` 校验阶段** 用 RDKit `MolFromSmiles` + `MolToSmiles(canonical=True)` 规范化；非法或无法解析 → 拒绝新版本。

### 3.6 非法 SMILES 处理（现状）

| 阶段 | 行为 |
|------|------|
| 父配体 SMILES 无效 | `FindGA.__init__` 抛 `ValueError`；`process_single_molecule_for_GA` **捕获后跳过该分子**（静默） |
| 环 SMILES 无效 | 记入 `non_kekulizable_smiles`，不进入 GA_dict |
| Kekulize 失败（芳香环） | 不加入 GA_dict |
| `generate_ga_id` 异常 | 返回 `None`，该行 **不写入** `GA_with_id.csv` |
| `calculate_gac` 无效分子 | 返回 GAC=0 或 worker 返回 `None` 跳过 |

Phase 1E 导入 CSV 应 **比现网更严格**（拒绝非法行），避免把脏数据写入不可变 version。

### 3.7 重复 GA_SMILES 处理（现状）

- **step1**：`merged_GA_dict` 以 `ring_smiles` 为 key，天然去重。
- **step2**：逐行 `generate_ga_id`，**无** 显式 drop_duplicates；若 temp 有重复 SMILES，会生成 **多行相同 SMILES + 相同 GA_ID**。
- **step6_7**：读入时 `dict[GA_SMILES]=GA_ID`，**后者覆盖前者**。

Phase 1E `create-ga-version-from-csv` 应：**按 canonical SMILES 去重**（保留首行或合并策略需在实施 spec 中固定为“首行 wins + 警告日志”）。

---

## 四、Workspace-level GA registry 设计

### 4.1 目录结构

```text
workspace/
└── ga_sets/
    └── {ga_set_id}/                    # 例如 p001_default, ni_library_v1
        ├── ga_set.json                 # set 级元数据（可变）
        └── versions/
            └── {version_id}/           # 例如 v001, v002；不可变
                ├── GA_with_id.csv      # 不可变快照
                └── ga_version.json     # 版本元数据
```

**约定**

- `ga_set_id`：`[a-z0-9][a-z0-9_-]{2,63}`（建议从 project/run 派生，可人工指定）。
- `version_id`：`v` + 三位递增数字（`v001`…）；**禁止覆盖**已有 version 目录。
- 无 `workspace/ga_sets/current.csv` 之类可变指针；run 只记录 `ga_set_id` + `ga_version_id` + checksum。

### 4.2 `ga_set.json`（建议 schema）

```json
{
  "ga_set_id": "p001_run_demo_derived",
  "name": "GA set from run_demo",
  "description": "",
  "created_at": "2026-05-21T12:00:00+00:00",
  "updated_at": "2026-05-21T12:00:00+00:00",
  "tags": ["auto", "p001"],
  "source": {
    "type": "from_run",
    "run_root": "workspace/projects/p001/runs/run_demo"
  }
}
```

### 4.3 `ga_version.json`（建议 schema）

```json
{
  "ga_set_id": "p001_run_demo_derived",
  "ga_version_id": "v001",
  "parent_version_id": null,
  "created_at": "2026-05-21T12:00:00+00:00",
  "created_by": "ga_registry_manager",
  "source": {
    "type": "from_run_generate_ga",
    "run_root": "workspace/projects/p001/runs/run_demo",
    "repaired_ligand_sha256": "..."
  },
  "num_ga": 126,
  "checksum": "sha256:abcdef...",
  "columns": ["GA_SMILES", "GA_ID"],
  "notes": ""
}
```

**checksum**：对 **物化用** CSV 计算（建议：按 `GA_SMILES` 排序后 UTF-8 序列化再 SHA256），与 `run/manifests/ga_binding.json` 中 checksum 算法一致。

### 4.4 `GA_with_id.csv`（workspace 版本）

**必需列**：`GA_SMILES`, `GA_ID`

**可选列**（仅 workspace 保留，不写入 run 物化文件）：`source`, `active`, `note`, `created_at`

**不可变规则**：`versions/{version_id}/` 一旦创建，**不原地修改**；编辑 → `create-ga-version-from-csv` → 新 `v00N`。

---

## 五、GA operations 设计

所有 CLI 建议：`python -m app.services.ga_registry_manager <subcommand> ...`  
实现位置：**`app/services/ga_registry_manager.py`**（不放入 `orchestrator.py`）。

### 5.1 `create-ga-set-from-run`

| 项 | 内容 |
|----|------|
| **输入** | `--run-root`（含 `tmp/repaired_ligand_data.csv`）；可选 `--ga-set-id`、`--version-id`（默认自动分配） |
| **行为** | 1) 在 run `tmp/` 或临时目录调用 `step4_5 --mode generate-ga -i/-o`；2) 将生成的 `GA_with_id.csv` **复制**到 `workspace/ga_sets/{ga_set_id}/versions/v001/`；3) 写 `ga_set.json`、`ga_version.json`；4) 可选：同时 `bind-ga-version-to-run`（见兼容策略） |
| **输出** | 新 ga set + v001 |
| **不写** | ChemDB 仓库 `tmp/`、`data/`、`result_analysis/` |

可与“从 run 已有 `GA_with_id.csv` 晋升”区分：若 run 已有 GA 且用户指定 `--from-existing-run-ga`，则跳过 generate，直接复制（便于迁移旧 run）。

### 5.2 `create-ga-version-from-csv`

| 项 | 内容 |
|----|------|
| **输入** | `--ga-set-id`；`--input-csv`（编辑后的 GA 表）；可选 `--parent-version-id` |
| **行为** | 校验 SMILES；canonicalize；去重 `GA_SMILES`；对缺 `GA_ID` 的行调用 `generate_ga_id`；校验 `GA_ID` 唯一；分配下一 `v00N`；写入新版本目录 |
| **输出** | `versions/vXXX/GA_with_id.csv` + `ga_version.json`；更新 `ga_set.json` 的 `updated_at` |
| **约束** | **不覆盖**旧 version |

### 5.3 `bind-ga-version-to-run`

| 项 | 内容 |
|----|------|
| **输入** | `--run-root`；`--ga-set-id`；`--ga-version-id` |
| **行为** | 1) 读 workspace 版本 CSV；2) **复制**到 `run_root/tmp/GA_with_id.csv`（两列物化）；3) 写 `manifests/ga_binding.json`；4) 若检测到 downstream stale 文件 → 写 `manifests/ga_stale_report.json` |
| **禁止** | symlink 到 workspace；不持有 workspace 文件句柄作为“当前 GA” |

### 5.4 `apply-ga-to-run`

| 项 | 内容 |
|----|------|
| **输入** | `--run-root`（需 `tmp/repaired_ligand_data.csv` + `tmp/GA_with_id.csv`） |
| **行为** | 调用 `step4_5 --mode apply-ga -i/-o {run}/tmp`；**不**写 workspace GA；**不**覆盖 workspace version |
| **输出** | `ligand_with_gac.csv`, `IRL_filtered.csv`, `IRL_filtered_cleaned.csv`, `step4_5_stats.json` |
| **前置** | 建议已 `bind`；若 GA 文件存在但无 `ga_binding.json`，允许 apply 但 warning |

---

## 六、Run GA binding metadata

### 6.1 `run_root/manifests/ga_binding.json`

```json
{
  "run_root": "/abs/path/to/run",
  "ga_set_id": "p001_run_demo_derived",
  "ga_version_id": "v001",
  "source_ga_csv": "workspace/ga_sets/p001_run_demo_derived/versions/v001/GA_with_id.csv",
  "materialized_path": "tmp/GA_with_id.csv",
  "checksum": "sha256:...",
  "num_ga": 126,
  "bound_at": "2026-05-21T12:00:00+00:00",
  "bound_by": "ga_registry_manager",
  "status": "bound"
}
```

| 字段 | 用途 |
|------|------|
| `checksum` | 与 workspace version 一致；用于检测 run 内文件是否被手工篡改 |
| `source_ga_csv` | 审计追溯；**运行时逻辑不得 open 此路径** |
| `status` | `bound` / `superseded`（重新绑定时可将旧 binding 标为 superseded 并追加历史，或单文件覆盖——实施时选 **单文件覆盖 + 可选 `binding_history` 数组**） |

### 6.2 是否需要 `run_config.json`？

**Phase 1E 不引入。** 理由：

- 当前 run 配置分散在 `manifests/*.json` + `training/config.yaml`；
- GA 绑定信息低频变更，**manifest 足够**；
- 避免与 orchestrator `run.json` 职责重叠。

若未来 Phase 2 需要统一 run 元数据，可再聚合；1E 仅 `ga_binding.json` + `ga_stale_report.json`。

---

## 七、Stale detection 设计

### 7.1 触发条件

- `bind-ga-version-to-run` 且新版本 checksum ≠ 上次 `ga_binding.json.checksum`（或首次 bind 但 downstream 已存在）
- `create-ga-version-from-csv` **不自动** stale run；仅影响未来 bind

### 7.2 应视为 stale 的 artifacts（相对 `run_root`）

**直接 GA/GAC/IRL（step4_5 apply 产出）**

- `tmp/ligand_with_gac.csv`
- `tmp/IRL_filtered.csv`
- `tmp/IRL_filtered_cleaned.csv`
- `tmp/ligand_data_deduplicated.csv`（若存在且与 step0 相关，可选标 stale）
- `tmp/step4_5_stats.json`

**Fragments / 图导出（step6_7–8）**

- `tmp/fragments.csv`
- `tmp/step6_7_stats.json`
- `tmp/neo4j_*.csv`（全部 neo4j 前缀输出）
- `tmp/step8_stats.json`

**Step12–13**

- `tmp/l3_l5.json`, `tmp/l5_l3.json`, `tmp/l5_freq_weight.json`
- `tmp/l3_gac.json`, `tmp/l5_l3_filtered_K30.json`
- `tmp/m_l3_pairs.csv`, `tmp/metal_l3_index.csv`
- `tmp/step12_stats.json`
- `tmp/step13_kl_nl_samples.csv`, `tmp/step13_kl_nl_samples.stats.csv`, `tmp/step13_stats.json`

**Training（run 内）**

- `data/L3_embedding/L3_embedding_ecfp.npz`
- `training/index.json`, `training/config.yaml`
- `training/split_index.json`, `training/train_records.pkl`, `training/val_records.pkl`, `training/test_records.pkl`
- `training/ckpts/best.pt`, `training/ckpts/last.pt`, `training/ckpts/history.json`

**Model-after（workspace，按 run 关联提示）**

- `workspace/model_after_results/**/{batch}/` 下凡引用该 `model_run_root` 的评估结果（路径模式：`**/p001_run_demo/**` 或 manifest 中记录的 run）

**不标 stale**

- `tmp/repaired_ligand_data.csv`（GA 提取输入，与 GA **集合** 选择无关，除非业务要求重提 GA）
- `tmp/complex_data.csv`, `tmp/ligand_data.csv`（step1–2）
- `data/pubchem/*`, `data/metal_embedding/*`（init-run 种子）

### 7.3 `manifests/ga_stale_report.json`（最小机制）

```json
{
  "generated_at": "...",
  "trigger": "bind-ga-version-to-run",
  "previous_binding": { "ga_version_id": "v001", "checksum": "..." },
  "new_binding": { "ga_version_id": "v002", "checksum": "..." },
  "stale_files": [
    {"path": "tmp/ligand_with_gac.csv", "exists": true, "reason": "downstream_of_ga"},
    {"path": "tmp/IRL_filtered.csv", "exists": true, "reason": "downstream_of_ga"}
  ],
  "recommended_actions": [
    "run: ga_registry_manager apply-ga-to-run",
    "run: orchestrator run-pipeline --through step13 (or training)"
  ]
}
```

- **不自动删除** 任何文件。
- orchestrator `run-step` **可不** 在 1E 强制拦截 stale（仅报告）；Phase 2 可加 `--require-fresh-ga`。

---

## 八、StepRegistry / Orchestrator 影响

### 8.1 建议新增 logical steps（名称供 registry / 文档使用）

| logical_step_id | 实现方式 | 是否进入 `PIPELINE_ORDER` |
|-----------------|----------|---------------------------|
| `create_ga_set_from_run` | `ga_registry_manager` in-process 或 subprocess | **否**（人工/工具触发） |
| `bind_ga_version_to_run` | 同上 | **否** |
| `apply_ga_to_run` | 调用 `step4_5 --mode apply-ga` | **可选**（见下） |
| `create_ga_version_from_csv` | `ga_registry_manager` | **否** |

### 8.2 与现有 `step4_5` 的关系

**方案 A（推荐，最小破坏）**

- 保留 `step_registry` 中 `step4_5` 为 **full-auto**（兼容现网 `run-pipeline --through step4_5`）。
- 文档约定新工作流：`bind` → `apply_ga_to_run`（单独命令）→ `run-pipeline --through step6_7`（或 `--through step13`）。
- `step4_5` full-auto 与“先 bind 再 apply”**互斥**；若 `ga_binding.json` 存在且用户跑 full-auto step4_5，stale 报告 + warning（1E 可选实现）。

**方案 B（更干净，改动略大）**

- 将 registry `step4_5` 拆为 `step4_5_apply_ga`（仅 apply-ga）+ 文档化 generate 仅用于 registry。
- core pipeline 顺序变为：`step2` → **`bind_ga`（前置条件）** → **`apply_ga`** → `step6_7` …
- 需要默认 bind 策略（见第十节）。

Phase 1E 实施推荐 **方案 A**；方案 B 放 Phase 1E.1。

### 8.3 `ga_registry_manager.py` 职责清单

- 路径解析：`WEB_ROOT` / `WORKSPACE_ROOT` / `GA_SETS_ROOT`
- CRUD：create set、create version、list versions
- 校验：RDKit SMILES、唯一 ID、checksum
- `bind`：物化 + `ga_binding.json`
- `stale_scan`：生成 `ga_stale_report.json`
- 子进程调用：`step4_5` generate-ga / apply-ga（`CHEMDB_REPO_ROOT` + run `tmp`）
- **禁止写**：`ChemDB/tmp/GA_with_id.csv`、`ChemDB/src/tmp/`、`ChemDB/data/`、`ChemDB/result_analysis/`

### 8.4 `orchestrator.py` 变更范围

- **不** 在 orchestrator 内实现 GA registry 业务逻辑。
- 可选：在 `init-run` 后文档提示绑定 GA；或提供 `run-pipeline` 前置 hook（**1E 不做**）。
- `run_step(step4_5)` 仍可原样调用；generate/bind 由独立 CLI 完成。

---

## 九、测试计划

**文件**：`tests/run_isolation/test_phase1e_ga_workspace.py`  
**隔离**：使用临时 `workspace` + tmp run_root；mock 或 mini CSV；禁止写 ChemDB 仓库路径。

| # | 测试名 | 验证点 |
|---|--------|--------|
| 1 | `test_create_ga_set_from_run_outputs_workspace_ga_version` | `workspace/ga_sets/{id}/versions/v001/GA_with_id.csv` 存在且 `num_ga>0` |
| 2 | `test_ga_version_is_immutable` | 编辑后 `create-ga-version-from-csv` → `v002`；`v001` checksum 不变 |
| 3 | `test_bind_ga_version_to_run_materializes_snapshot` | `run/tmp/GA_with_id.csv` 字节等于 version CSV（或 canonical 相等） |
| 4 | `test_run_binding_does_not_reference_mutable_file` | `ga_binding.json` 含 `ga_version_id`+`checksum`；删除 workspace CSV 后 run 内仍可用 |
| 5 | `test_create_ga_version_from_csv_validates_smiles` | 非法 SMILES → 非零退出且无新 version 目录 |
| 6 | `test_create_ga_version_from_csv_generates_missing_ga_id` | 仅 `GA_SMILES` 行 → 自动 `GA_ID` |
| 7 | `test_apply_ga_to_run_uses_materialized_ga` | patch/monitor：apply 读 `run/tmp/GA_with_id.csv`；改 workspace 不影响 GAC 输出 |
| 8 | `test_binding_new_ga_version_marks_downstream_stale` | 预建 `ligand_with_gac.csv`、`IRL_filtered.csv`、`step13_kl_nl_samples.csv` → rebind → `ga_stale_report.json` 列出 |
| 9 | `test_two_runs_can_bind_different_ga_versions` | run_a / run_b 不同 checksum，互不影响 |
| 10 | `test_no_write_to_repo_ga_registry` | 断言 ChemDB 仓库内无新增/修改 GA 相关文件 |

**pytest 标记**：`@pytest.mark.integration_ga_workspace`（可选 deselect 默认 CI）。

**依赖**：需 RDKit；generate-ga 可 `@pytest.mark.slow` 用小 CSV fixture。

---

## 十、兼容性要求

| # | 要求 | 设计对策 |
|---|------|----------|
| 1 | 旧 full-auto pipeline 尽量可用 | `step4_5` 默认 `--mode full-auto`；现有 orchestrator 命令不变 |
| 2 | 用户不选手动 GA set 时可自动推导 | **推荐流程**：首次 `run-pipeline --through step4_5` 后执行 `create-ga-set-from-run` + `bind`（文档/脚本）；或提供 `--auto-publish-ga` 一次性标志（1E 可选） |
| 3 | 后续 step 仍读 `run/tmp/GA_with_id.csv` | bind 物化保证 |
| 4 | 无 SQLite | 仅 JSON + CSV |
| 5 | 无 Streamlit | 不涉及 |
| 6 | 不改 training/model_after | 仅 stale 提示重跑 |
| 7 | 无 DB 表 | 无 |

### 10.1 推荐用户工作流（Phase 1E 后）

**新 run（显式 GA 管理）**

```text
init-run → core through step2
  → create-ga-set-from-run（或 bind 已有 set）
  → bind-ga-version-to-run
  → apply-ga-to-run
  → orchestrator run-pipeline --through step13
  → training pipeline ...
```

**旧 run 迁移**

```text
自 run_demo/tmp/GA_with_id.csv 创建 ga_set v001 → bind 同版本 → 标记已有 downstream 是否 stale
```

**兼容旧命令**

```text
orchestrator run-pipeline --through step4_5   # 仍一次性生成 GA+GAC+IRL
# 之后可 create-ga-set-from-run --from-existing-run-ga 晋升为 workspace 版本
```

---

## 十一、实施顺序建议（供开发阶段使用）

1. **`step4_5.py`**：实现 `--mode` 三套入口 + 单元测试（mock RDKit 可选）。
2. **`ga_registry_manager.py`**：路径、checksum、create/bind/version-from-csv。
3. **`ga_stale_report`**：扫描逻辑 + bind 集成。
4. **`test_phase1e_ga_workspace.py`**：十项测试。
5. **文档**：更新 `_auditing/README.md`、`ChemDB/doc/PIPELINE_IO.md`（GA 章节）、`step_registry` 注释。
6. **（可选）** orchestrator 文档字符串 + `phase_1e_implementation.md` 实施记录。

---

## 十二、风险与开放问题

| 风险 | 缓解 |
|------|------|
| `GA_ID` 哈希碰撞 | 导入时唯一性检查；冲突 fail |
| 手工编辑 SMILES 未 canonical | registry 强制 canonicalize |
| full-auto 与 bind 双源 GA 不一致 | 文档 + stale；长期用方案 B 拆分 step4_5 |
| `step6_7` 缺失 GA 文件时静默 `{}` | bind 作为必需前置；测试覆盖 |
| 大库 generate-ga 耗时 | workspace 版本复用；run 只 bind |

**开放问题（实施前需产品确认）**

1. `create-ga-set-from-run` 是否默认 **同时 bind 回源 run**？
2. 物化 CSV 是否 **严格两列**，还是保留额外列给未来 UI？
3. rebind 后是否 **强制** orchestrator 拒绝 stale step（1E 仅报告 vs 硬失败）？

---

## 附录：step4_5 mode 与 orchestrator 命令对照

| 用户意图 | 命令 |
|----------|------|
| 从配体库 **发现** GA 模式 | `step4_5 --mode generate-ga` 或 `ga_registry_manager create-ga-set-from-run` |
| 将 workspace 版本 **固定到 run** | `ga_registry_manager bind-ga-version-to-run` |
| 用 run 内 GA **算 GAC/IRL** | `step4_5 --mode apply-ga` 或 `ga_registry_manager apply-ga-to-run` |
| 一条命令全做（旧行为） | `step4_5` 或 `orchestrator run-step --step step4_5` |

---

*本文档仅审计与计划；代码实现见后续 Phase 1E 开发任务。*
