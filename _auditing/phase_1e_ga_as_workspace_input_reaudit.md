# Phase 1E 重审：GA 作为 Workspace 级输入资源（语义变更审计与改造计划）

**日期**：2026-05-21  
**状态**：仅审计与计划（**不改代码、不改测试、不实现 Phase 1E**）  
**依据**：当前仓库真实代码与测试（Phase 1A–1D 已落地），非旧计划臆测。

---

## 语义变更摘要

| 维度 | 旧语义（当前实现） | 新语义（目标） |
|------|-------------------|----------------|
| GA 身份 | `step4_5` 在 run 内**自动生成**的 `tmp/GA_with_id.csv` 中间产物 | **Workspace 级**可版本化、可编辑、可复用的**输入资源** |
| Run 关系 | 无显式绑定；pipeline 跑到 step4_5 即有 GA | 每个 run **必须绑定**某一 `ga_version`；物化为 `run/tmp/GA_with_id.csv` 冻结快照 |
| Core pipeline | `step1 → step2 → step4_5 → …`（step4_5 必经且含 GA 生成） | **不含**自动 GA 生成；缺绑定时**明确失败** |
| step4_5 角色 | Core 默认一步 | **Authoring tool**：`generate-ga` / `apply-ga`；`full-auto` 仅 legacy/debug |

---

## 一、当前 GA 依赖现状（全仓库审计）

搜索范围：`ChemDBWebVersion/`（含 `ChemDB/src`、`app/services`、`tests`、`_auditing`、workspace 样例 manifest）。  
符号：`ga_output` 仅出现在 `step4_5.py` 的 `self.ga_output` 属性（路径 `GA_with_id.csv`）。

### 1.1 总表（按文件）

| file_path | symbol/function | read/write | current path assumption | run-based? | depends on auto-generated GA? | required change |
|-----------|-----------------|------------|-------------------------|------------|------------------------------|-----------------|
| **直接 GA_with_id.csv** |
| `ChemDB/src/step4_5.py` | `Step4_5Processor.step2_generate_GA_ID` | **W** | `{output_dir}/GA_with_id.csv`，默认 `./tmp` | ✅ `-i/-o` | 自身生成 | 拆 `mode`；core 不再默认调用 generate |
| `ChemDB/src/step4_5.py` | `step3_calculate_GAC` | **R** | 读 `self.ga_output` | ✅ | 依赖同次 run 的 step2 输出 | `apply-ga` 只读已有文件 |
| `ChemDB/src/step6_7.py` | `Step6_7Processor._load_ga_data` | **R** | `{input_dir}/GA_with_id.csv` | ✅ orchestrator 传 `{tmp_dir}` | **是**（core 先跑 step4_5） | 无路径改；前置 `require_bound_ga` |
| `ChemDB/src/step6_7.py` | `GA_PATH` 模块常量 | 文档/遗留 | `./tmp/GA_with_id.csv` | ❌ 常量 | — | 无功能改（ctor 已覆盖） |
| `app/services/step_registry.py` | `step4_5` expected_outputs | 契约 **W** | `tmp/GA_with_id.csv` | ✅ | **是** | 从 core 移除或改为 apply 产出 |
| `app/services/step_registry.py` | `step6_7` required_inputs | 契约 **R** | `tmp/GA_with_id.csv` | ✅ | **是** | 增加 `manifests/ga_binding.json` 前置检查 |
| `tests/.../test_integration_1pct.py` | `test_1pct_tier_b_*` | 断言存在 | `ctx.tmp_dir/*.csv` | ✅ | **是** | 重写：bind fixture + apply-ga |
| `workspace/.../manifests/step4_5.json` | 运行记录 | — | 记录绝对路径 | ✅ | 样例来自旧语义 | 新 run 应有 `ga_binding.json` |
| **间接（GAC / IRL / fragments / …）** |
| `ChemDB/src/step4_5.py` | `step3` → `ligand_with_gac.csv` | **W** | `{output_dir}/ligand_with_gac.csv` | ✅ | **是** | `apply-ga` 产出；GA 变 → stale |
| `ChemDB/src/step4_5.py` | `step4/5` → `IRL_filtered*.csv` | **W** | `{output_dir}/` | ✅ | **是** | 同上 |
| `ChemDB/src/step6_7.py` | `_load_irl_data` | **R** | `{input_dir}/IRL_filtered.csv` | ✅ | 依赖 step4_5 apply | bind+apply 后再跑 |
| `ChemDB/src/step6_7.py` | `process_all_ligands` → `fragments.csv` | **W** | `{output_dir}/fragments.csv` | ✅ | 间接 | stale |
| `ChemDB/src/step8.py` | 读 `fragments`, `IRL_filtered` | **R** | `{input_dir}/` | ✅ | 间接 | stale |
| `ChemDB/src/step12.py` | `extract_l3_gac` | **R** | `ligand_with_gac.csv` | ✅ | 间接 | stale |
| `ChemDB/src/step12.py` | 多 JSON/CSV 产出 | **W** | `tmp/l3_l5.json` 等 | ✅ | 间接 | stale |
| `ChemDB/src/step_13_complete.py` | `l3_gac.json`, samples | **R/W** | `{input_dir}/` | ✅ | 间接 | stale |
| `app/services/step_registry.py` | step8/12/13 契约 | R/W 链 | 相对 run_root | ✅ | 间接 | 顺序调整 |
| `app/services/training_manager.py` | `prepare_training_index` | **R** | `tmp/step13_*`, `l3_gac.json` 等 | ✅ | 间接 | 无 GA 字段；可选 metadata |
| `ChemDB/training/**` | 训练脚本 | — | run `training/` | ✅ | 间接 | 1E 不改逻辑；stale 提示 |
| `ChemDB/training/model_after/**` | bundle/scoring | — | `training/index.json`, ckpt | ✅ | 间接 | 1E 不改逻辑；metadata 可选 |
| **GA 算法 / ID（无直接 CSV 路径）** |
| `ChemDB/src/step4_5.py` | `FindGA`, `process_single_molecule_for_GA` | — | — | — | generate 路径 | 供 `generate-ga` / authoring |
| `ChemDB/src/step4_5.py` | `generate_ga_id`, `calculate_gac` | — | — | — | generate/apply | 复用，不改算法 |
| `ChemDB/src/tools/L4_create.py` | `ga_list` 参数 | **R**（内存） | 调用方传入 | — | step6_7 间接 | 无 |
| `ChemDB/src/tools/L4_create_independent.py` | 同上 | **R** | — | — | — | 无 |
| **Legacy / 非 WebVersion 主路径** |
| `ChemDB/src/main.py` | `STEP_ORDER` 含 step4_5 | — | `./tmp` | ❌ | **是** | Web 不用；文档标注 legacy |
| `ChemDB/src/step9.py` | `data/IRL_filtered.csv` | **R** | `./data` + `./tmp` | ❌ | — | 不在 `PIPELINE_ORDER` |
| `ChemDB/src/step10.py` | `IRL_filtered` | **R** | `data_dir` | ❌ | — | 不在 orchestrator |
| `ChemDB/src/server.py` | IRL csv_path | **R** | `./tmp/IRL_filtered.csv` | ❌ | — | 1E 不做 Streamlit/server |
| `ChemDB/doc/PIPELINE_IO.md` | 文档 | — | `./tmp` | — | 描述旧语义 | 实施后更新 |
| `_auditing/*.md` | 历史审计 | — | — | — | — | 本文件 supersede 部分 1E 计划 |
| **仅测试/文档引用** |
| `tests/.../test_integration_1pct.py` | tier_b/c 断言 GA 文件 | — | tmp | ✅ | **是** | **必须重写** |
| `tests/.../test_orchestrator.py` | 仅 mock step2 | — | — | ✅ | 否 | 保留 + 新增 1E 测试 |
| `tests/.../test_phase1c_*.py` | 手工 seed，无 GA | — | — | ✅ | 否 | 大部分保留 |
| `tests/.../test_phase1d_*.py` | 无 GA | — | — | ✅ | 否 | 保留；metadata 可选 |
| **WebVersion 实际会调用** |
| `app/services/orchestrator.py` | `run_pipeline` / `run_step` | 调度 | `{tmp_dir}` 展开 | ✅ | **是**（经 step4_5） | **核心改造点** |
| `app/services/__main__.py` | → `orchestrator.main` | — | — | ✅ | — | 无 |
| `app/services/model_after_manager.py` | evaluate/recommend | — | workspace results | ✅ | 否（间接依赖训练链） | metadata 可选 |

### 1.2 分类结论

1. **直接读写 `GA_with_id.csv`**：仅 `step4_5`（写+读）、`step6_7`（读）、`step_registry`（契约）、1% 集成测试（断言）。
2. **间接依赖**：`ligand_with_gac` → `IRL_filtered` → `fragments` → `neo4j_*` → step12/13 → embedding → training → model_after。
3. **测试/文档**：`test_integration_1pct` tier B/C；`PIPELINE_IO.md`、`_auditing` 旧文。
4. **Legacy/unused（WebVersion orchestrator）**：`main.py` step9/10/11、`server.py`；`GA_PATH` 常量已被 `input_dir` 覆盖。
5. **WebVersion 实际路径**：`python -m app.services.orchestrator` + `python -m app.services.model_after_manager`；**不**通过 `ChemDB/src/main.py`。

---

## 二、审计 `ChemDB/src/step4_5.py`

### 2.1 是否同时完成 GA extraction、ID、GAC、IRL 及后半段产物？

**是。** `process_all()` 顺序执行：

```806:836:ChemDBWebVersion/ChemDB/src/step4_5.py
    def process_all(self) -> Dict:
        ...
        if success:
            success = self.step0_deduplicate()
        if success:
            success = self.step1_extract_GA()
        if success:
            success = self.step2_generate_GA_ID()
        if success:
            success = self.step3_calculate_GAC()
        if success:
            success = self.step4_extract_IRL()
        if success:
            success = self.step5_deduplicate_IRL()
```

单次运行产出：`GA_with_id.csv`、`ligand_with_gac.csv`、`IRL_filtered.csv`、`IRL_filtered_cleaned.csv`（及可选 `ligand_data_deduplicated.csv`、`step4_5_stats.json`）。

### 2.2–2.7 函数职责（实测）

| # | 问题 | 答案 |
|---|------|------|
| 2 | 从 `repaired_ligand_data.csv` 提取 GA candidates | `_get_smiles_DID_from_csv()` → `step1_extract_GA()` → 并行 `process_single_molecule_for_GA()` + `FindGA` |
| 3 | 生成 GA_ID | `generate_ga_id()`；`step2_generate_GA_ID()` |
| 4 | 写 `GA_with_id.csv` | `step2_generate_GA_ID()` → `ga_id_df.to_csv(self.ga_output)`，`ga_output = output_dir / "GA_with_id.csv"` |
| 5 | 读 `GA_with_id.csv` | `step3_calculate_GAC()`：`pd.read_csv(self.ga_output)` 取 `GA_SMILES` 列表 |
| 6 | 写 `ligand_with_gac.csv` | `step3_calculate_GAC()` → `gac_df.to_csv(self.gac_output)` |
| 7 | IRL 等后半段 | `step4_extract_IRL()` → `IRL_filtered.csv`；`step5_deduplicate_IRL()` → `IRL_filtered_cleaned.csv` |

### 2.8 当前 CLI 是否支持 mode？

**不支持。** `main()` 仅有 `--workers`、`--limit`、`--input-dir`、`--output-dir`，固定调用 `processor.process_all()`：

```873:902:ChemDBWebVersion/ChemDB/src/step4_5.py
def main():
    ...
    parser.add_argument('--input-dir', '-i', ...)
    parser.add_argument('--output-dir', '-o', ...)
    ...
    stats = processor.process_all()
```

### 2.9 最小改造方案（不改化学算法）

| mode | 方法 | 步骤 | 约束 |
|------|------|------|------|
| `generate-ga` | 新增 `process_generate_ga()` | step0（建议保留）→ step1 → step2 | **不**调用 step3–5；**不**要求事先存在 `GA_with_id.csv` |
| `apply-ga` | 新增 `process_apply_ga()` | step0（可选）→ step3 → step4 → step5 | step3 前 `assert self.ga_output.is_file()`；**禁止** step1–2 覆盖已有 `GA_with_id.csv` |
| `full-auto` | 保留 `process_all()` | 全流程 | legacy/debug；**orchestrator core 不得默认调用** |

实现量：约 30–50 行调度 + `argparse --mode` + `main()` 分支；**零**改动 `FindGA` / `calculate_gac` / `ExtractIRLObject` 内部逻辑。

---

## 三、Phase 1A–1D 受影响范围

### 3.1 Phase 1A

**审计文件**：`run_context.py`、`step_registry.py`、`orchestrator.py`、`__main__.py`、`tests/run_isolation/*`

| 问题 | 现状（代码） |
|------|----------------|
| Core pipeline 默认顺序 | `PIPELINE_ORDER = ["step1","step2","step4_5","step6_7","step8","step12","step13"]`（`step_registry.py` L29–37） |
| step4_5 是否必经？ | **是**。`pipeline_through("step6_7")` 含 step4_5；`run_pipeline(..., "step4_5")` 为 1B Tier B 验收 |
| required/expected 假设 | step4_5：`required_inputs=["tmp/repaired_ligand_data.csv"]`，`expected_outputs` 含 **`tmp/GA_with_id.csv`** + gac + irl（L80–83） |
| 若 GA 必须先绑定，registry 如何变？ | 见第六节；建议：`step4_5` **移出** `PIPELINE_ORDER`；新增 `require_bound_ga`（in-process）、`apply_ga_to_run`（subprocess `step4_5 --mode apply-ga`） |
| 是否新增 logical steps？ | **是**：`require_bound_ga`、`apply_ga_to_run`；**否**（不把 `materialize_bound_ga` 单独列为 pipeline step，合并进 bind CLI） |
| Phase 1A 测试需改？ | `test_orchestrator.py`、`test_step_registry.py`、`test_run_context.py`、`test_no_write_*`：**可保留**；**必须改/增**：`test_integration_1pct.py` tier B/C；新增 `test_phase1e_ga_workspace.py` |

### 3.2 Phase 1B（1% integration）

| 问题 | 现状 |
|------|------|
| 是否依赖 step4_5 自动生成 GA？ | **是**。`test_1pct_tier_b_pipeline_through_step4_5`：`run_pipeline(ctx, "step4_5")` 并断言 `GA_with_id.csv` 等（L65–71） |
| Tier C | `run_pipeline(ctx, "step6_7")` **仍包含** step1→step2→**step4_5**→step6_7 |
| 改造方向 | 拆为：preprocessing（step1+2）→ **bind workspace GA fixture** → **apply-ga** → downstream；或 `through step6_7` 前注入 bind+apply |
| GA fixture | **需要** `tests/run_isolation/fixtures/ga_sets/.../v001/GA_with_id.csv`（小表即可）或 tier B 后 `generate-ga-from-run` 晋升 |
| 不建议 | 继续默认 `through step4_5` 作为 core 语义 |

### 3.3 Phase 1C

| 问题 | 现状 |
|------|------|
| 直接读 `GA_with_id.csv`？ | **否**。`ChemDB/training/**` 无 GA 字符串匹配 |
| index/config 记录 GA？ | **否**。`training_index_template.json` 含 `l3_gac`、`samples_csv` 等，**无** `ga_set_id` / `ga_version_id` |
| 模型 metadata 应记录 GA？ | **建议**（1E-standard）：`manifests/ga_binding.json` 副本或 `training/manifest` 字段；MVP 可仅 run 级 binding |
| GA binding 变更后 stale | `data/L3_embedding/L3_embedding_ecfp.npz`；`training/index.json`；`config.yaml`；`split_index.json`；`*_records.pkl`；`ckpts/best.pt`；`last.pt`；`history.json` |
| 测试 | `test_phase1c_*` 用 `_seed_training_inputs()` **绕过** 真实 GA pipeline → **可保留**；集成训练测试不依赖 step4_5 |

### 3.4 Phase 1D

| 问题 | 现状 |
|------|------|
| 直接读 GA？ | **否** |
| results 记录 GA version？ | **否**。`evaluate_model_manifest.json` 仅有 `model_run_root`、`checkpoint_path` 等 |
| summary/best_model 加 GA？ | **建议** 1E-standard；MVP 不改 `model_after` 代码，仅在文档说明间接依赖 |
| evaluate-models 展示 GA？ | 可选：读各 run 的 `ga_binding.json` 写入 summary 列（registry 层，非 scoring 逻辑） |
| model_runs.csv | 当前列：`model_id, model_run_root, checkpoint_path, ...`（batch 实测）；**建议** 可选列 `ga_set_id, ga_version_id, ga_binding_checksum` |

### 3.5 总结表

| phase | affected? | code changes needed | metadata/schema changes needed | tests to update | priority |
|-------|-----------|---------------------|-------------------------------|-----------------|----------|
| **1A** | **是** | `step_registry` 顺序；`orchestrator` 前置检查；新增 `ga_registry_manager.py`；`step4_5 --mode` | `ga_binding.json`、`ga_stale_report.json` | tier B/C；新增 1E 套件 | **P0** |
| **1B** | **是** | 集成流程 bind+apply；GA fixture | fixture `ga_version.json` | `test_integration_1pct.py` | **P0** |
| **1C** | 间接 | 默认无；可选 training manifest 引用 binding | 可选 index/manifest 字段 | 多数保留 | P1 |
| **1D** | 间接/可选 | MVP 无；standard 为 results 增 GA 列 | manifest/summary 可选 | 保留；可选增断言 | P2 |

---

## 四、Workspace-level GA registry（file-based）

```text
workspace/
└── ga_sets/
    └── {ga_set_id}/
        ├── ga_set.json
        └── versions/
            └── {ga_version_id}/          # 建议 v001, v002… 只增不改
                ├── GA_with_id.csv
                └── ga_version.json
```

### 4.1 `ga_set.json`（最小字段）

`ga_set_id`, `name`, `description`, `created_at`, `updated_at`, `tags`, `source`（object：`type`, `run_root`, `notes` 等）

### 4.2 `ga_version.json`（最小字段）

`ga_set_id`, `ga_version_id`, `parent_version_id`, `created_at`, `created_by`, `source`, `num_ga`, `checksum`, `notes`, `columns`

### 4.3 `GA_with_id.csv`

- **必需**：`GA_SMILES`, `GA_ID`（与现网 `run_demo` 一致，126 行实测无重复）
- **可选**（仅 workspace）：`source`, `active`, `note`, `created_at`
- **物化规则**：复制到 `run/tmp/GA_with_id.csv` 时默认 **仅两列**（或保留额外列但 step6_7/step3 忽略——推荐 strip 保证兼容）

### 4.4 checksum

建议：对物化用两列 CSV 按 `GA_SMILES` 排序后 UTF-8 序列化再 SHA256；写入 `ga_version.json` 与 `ga_binding.json`。

---

## 五、Run-level GA binding

路径：`{run_root}/manifests/ga_binding.json`

| 字段 | 说明 |
|------|------|
| `run_root` | 绝对路径 |
| `ga_set_id` / `ga_version_id` | 绑定版本 |
| `source_ga_csv` | workspace 内版本文件路径（**审计用**，运行时只读 `materialized_path`） |
| `materialized_path` | 固定 `tmp/GA_with_id.csv` |
| `checksum` | 与 version 一致 |
| `num_ga` | 行数 |
| `bound_at` / `bound_by` / `status` | `bound` / `superseded` |

**bind-ga-version-to-run** 行为（与需求对齐）：

1. 读 workspace 版本 CSV  
2. 校验存在、schema、checksum  
3. **复制**到 `run/tmp/GA_with_id.csv`（非 symlink）  
4. 写 `ga_binding.json`  
5. 若 downstream 存在 → `ga_stale_report.json`  
6. 不写 ChemDB 源码目录  

**不引入 `run_config.json`**（当前无此文件；`run.json` 仅 init-run）。

---

## 六、新 Core Pipeline 语义

### 6.1 现状

```text
step1 → step2 → step4_5 → step6_7 → step8 → step12 → step13
```

`orchestrator run-pipeline --through step4_5` 会执行 step1+2+4_5 并**自动生成** GA/GAC/IRL。

### 6.2 目标

```text
step1 → step2 → require_bound_ga → apply_ga_to_run → step6_7 → step8 → step12 → step13
```

| # | 问题 | 建议 |
|---|------|------|
| 1 | step4_5 是否从 core 移除？ | **是**（从 `PIPELINE_ORDER` 移除 `step4_5`） |
| 2 | apply_ga_to_run 是否替代原 step4_5？ | **是**（功能上替代 **apply** 半段；不替代 generate） |
| 3 | require_bound_ga 是否独立 logical step？ | **是**（in-process：检查 `ga_binding.json` + `tmp/GA_with_id.csv` + checksum 一致） |
| 4 | generate-ga-from-run 是否在 core 外？ | **是**（`ga_registry_manager` authoring tool） |
| 5 | `--through step4_5` 兼容？ | **废弃为 core through**；可保留 alias 映射到 `apply_ga_to_run` 并 **要求已 bind**；或仅 `run-step --step step4_5 --mode full-auto` 作 debug |
| 6 | 测试命令更新 | 见第九节；文档中 `through step4_5` → `through apply_ga_to_run` 或分步 bind+apply |

### 6.3 无绑定时行为

`require_bound_ga` 失败示例：

- 缺少 `manifests/ga_binding.json`  
- 缺少 `tmp/GA_with_id.csv`  
- checksum 与 binding 不一致  

**禁止**：静默调用 `generate-ga` 或 `full-auto`。

---

## 七、GA Authoring Tool（`ga_registry_manager.py`）

建议：`python -m app.services.ga_registry_manager <cmd>`

| 操作 | CLI 名 | 输入 | 输出 | 备注 |
|------|--------|------|------|------|
| 1 | `generate-ga-from-run` | `run/tmp/repaired_ligand_data.csv` | workspace `ga_sets/.../v001/` + json | 调 `step4_5 --mode generate-ga`；**默认不 bind**；`--bind-to-run` 才 bind |
| 2 | `create-ga-version-from-csv` | 编辑 CSV + `ga_set_id` | 新 `v00N/` | 校验/canonicalize/去重/补 GA_ID；不覆盖旧版 |
| 3 | `bind-ga-version-to-run` | run + set + version | 物化 + `ga_binding.json` + 可选 stale | 见第五节 |
| 4 | `apply-ga-to-run` | run + 已物化 GA | gac/irl 等 | `step4_5 --mode apply-ga`；不碰 workspace |

**与 orchestrator 边界**：registry 管 workspace + bind + stale；orchestrator 管 core/training pipeline 步骤调度。

---

## 八、Stale 规则

**触发**：`bind-ga-version-to-run` 且 checksum 相对上次 binding 变化（或首次 bind 但 downstream 已存在）。

**列入 `ga_stale_report.json`（存在则记录，不删除）**：

- `tmp/ligand_with_gac.csv`, `IRL_filtered.csv`, `IRL_filtered_cleaned.csv`, `step4_5_stats.json`
- `tmp/fragments.csv`, `step6_7_stats.json`
- `tmp/neo4j_*.csv`, `step8_stats.json`
- `tmp/l3_l5.json`, `l5_l3.json`, `l5_freq_weight.json`, `l3_gac.json`, `l5_l3_filtered_K*.json`, `m_l3_pairs.csv`, `metal_l3_index.csv`, `step12_stats.json`
- `tmp/step13_kl_nl_samples.csv`, `step13_*.stats.csv`, `step13_stats.json`
- `data/L3_embedding/L3_embedding_ecfp.npz`
- `training/index.json`, `config.yaml`, `split_index.json`, `*_records.pkl`, `ckpts/*`
- `workspace/model_after_results/**` 下可解析到该 `model_run_root` 的目录（若存在）

**不列入**：`repaired_ligand_data.csv`（除非产品要求重提 GA）、`complex_data.csv`、`ligand_data.csv`、pubchem 种子。

**Phase 1E**：只报告，不自动删，不强制阻止 orchestrator（hard fail 可放 1E-standard）。

---

## 九、测试矩阵（新增/修改）

文件：`tests/run_isolation/test_phase1e_ga_workspace.py`

### A. GA registry

| 测试 | 目的 |
|------|------|
| create GA set from run | workspace v001 存在 |
| immutable version | v002 不覆盖 v001 |
| invalid SMILES rejected | 无新 version |
| missing GA_ID generated | 补 ID |
| checksum recorded | ga_version.json |

### B. Run binding

| 测试 | 目的 |
|------|------|
| bind materializes snapshot | `tmp/GA_with_id.csv` |
| ga_binding.json written | 字段完整 |
| checksum mismatch detected | 篡改 run 文件失败 |
| same version two runs | 独立副本 |
| different versions | 互不影响 |

### C. Core pipeline

| 测试 | 目的 |
|------|------|
| pipeline without binding **fails clearly** | `require_bound_ga` |
| **does not auto-generate GA** | 无 bind 时无新 GA 文件 |
| bind + apply + step6_7 | 可继续 |
| step6_7 requires materialized GA | 缺文件失败 |
| no write ChemDB/tmp | 与现 guard 一致 |

### D. Stale

| 测试 | 目的 |
|------|------|
| rebind → ga_stale_report.json | 列出 gac/irl/step13 等 |
| no auto deletion | 文件仍在 |

### E. Backward compatibility

| 测试 | 目的 |
|------|------|
| step4_5 full-auto legacy | `run-step` + `--mode full-auto` 仍产出三文件 |
| identify old tests | **仅** `test_1pct_tier_b`、tier C 依赖 auto GA |
| non-GA 1A–1D tests | orchestrator mock、phase1c/d 保留 |

### F. Integration

| 测试 | 目的 |
|------|------|
| minimal bound run past apply-ga | 小 fixture |
| 1% integration | tier B：step1–2 → bind fixture → apply → assert；tier C：同前 + step6_7 |

### 需修改的现有测试（明确列表）

| 文件 | 测试 | 原因 |
|------|------|------|
| `test_integration_1pct.py` | `test_1pct_tier_b_pipeline_through_step4_5` | 假设 auto GA |
| `test_integration_1pct.py` | `test_1pct_tier_c_pipeline_through_step6_7` | pipeline 含 step4_5 |
| `test_integration_1pct.py` | `test_1pct_enhanced_step8` | pre pipeline 含 step4_5 |

**可保留**：`test_orchestrator.py`、`test_step_registry.py`（registry 变更后需更新路径断言若 step id 变化）、`test_phase1c_*`、`test_phase1d_*`、`test_integration.py`（仅 step1–2）、`test_two_runs.py`、`test_no_write_*`。

---

## 十、最终输出建议

### 10.1 影响的 Phase

| Phase | 影响程度 |
|-------|----------|
| 1A | **结构性**（pipeline 顺序、失败语义） |
| 1B | **集成路径**（tier B/C） |
| 1C | **间接 + metadata 可选** |
| 1D | **间接 + metadata 可选** |

### 10.2 必须做的修改（P0）

1. `step4_5.py`：`--mode generate-ga | apply-ga | full-auto`  
2. `ga_registry_manager.py`：generate-from-run、create-version-from-csv、bind、apply  
3. `step_registry.py`：移除 core 中 `step4_5`；加 `require_bound_ga`、`apply_ga_to_run`  
4. `orchestrator.py`：`run_in_process_step` 处理 `require_bound_ga`；`apply_ga` 子进程  
5. `test_phase1e_ga_workspace.py` + **重写** `test_integration_1pct` tier B/C  

### 10.3 仅 metadata/schema（可 MVP 后做）

- `ga_binding.json` / `ga_version.json`（MVP **必须**有 binding；version metadata 完整字段 MVP 可有）  
- training / model_after manifest 中的 `ga_set_id`（**standard**）  
- `model_runs.csv` 可选 GA 列（**standard**）  

### 10.4 必须重写的测试

- `test_1pct_tier_b_pipeline_through_step4_5`（改为 bind+apply 或拆分 pipeline）  
- `test_1pct_tier_c_*`、`test_1pct_enhanced_step8`（pipeline 起点变更）  

### 10.5 可保留的测试

- `test_orchestrator`（mock step2）  
- `test_run_context`、`test_two_runs`  
- `test_phase1c_*`（人工 seed，不走过 GA）  
- `test_phase1d_*`  
- `test_integration_minimal_step1_step2`  
- repo tmp guard（`conftest` autouse）  

### 10.6 推荐最小实施路线

```text
Phase 1E-MVP（按序）
  1. step4_5 --mode（generate / apply / full-auto）
  2. ga_registry_manager：bind + generate-ga-from-run + apply-ga-to-run
  3. step_registry：require_bound_ga + apply_ga_to_run；移除 PIPELINE_ORDER 中的 step4_5
  4. orchestrator in-process 检查
  5. test_phase1e（A+B+C 核心）+ 修 tier B integration + GA fixture
  6. 文档：TESTING.md、phase_1e_implementation.md

Phase 1E-standard（后续）
  7. create-ga-version-from-csv
  8. ga_stale_report 完整扫描 + model_after_results 关联
  9. training/model_after manifest GA 字段
  10. 废弃 --through step4_5 的迁移说明与 CLI alias
```

### 10.7 Phase 1E-MVP 边界

| 在 MVP 内 | 不在 MVP 内 |
|----------|------------|
| workspace `ga_sets/` 单 set 多 version | UI / Streamlit |
| bind + materialize + checksum | SQLite |
| core 无 bind 即失败 | 自动从 run 晋升 GA（无 `--bind-to-run` 则不绑） |
| apply-ga + require_bound_ga | model_after 代码改动 |
| step4_5 mode 三件套 | 强制 stale 阻断 pipeline |
| 1E 测试 + tier B 修复 | CSV 在线编辑工具 |
| 小 GA fixture | 全量 1% run 上自动生成 ga_set |

### 10.8 Phase 1E-standard（以后）

- `create-ga-version-from-csv` 完整校验 UX  
- `model_runs.csv` / evaluation manifest GA 列  
- rebind 硬失败选项  
- 从既有 `run_demo` 批量迁移 ga_set  
- `PIPELINE_IO.md` / operator runbook  

### 10.9 进入 SQLite（Phase 2）前必须完成

1. File-based `ga_sets` + `ga_binding` **语义稳定**且集成测试通过  
2. Core pipeline **不再**隐式生成 GA（避免 DB 与文件双源冲突）  
3. Stale 规则文档化（即便 Phase 2 用 DB 存 binding，失效链不变）  
4. 至少一个 **固定 GA fixture version** 供 CI 1% 集成使用  
5. `step4_5` mode 与 authoring CLI 稳定（DB 层只替换 JSON/CSV 存储，不替换化学逻辑）

### 10.10 是否先 1E-MVP 再 SQLite Phase 2？

**建议：是。**

理由（基于现状）：

- 当前 GA 与 `step4_5` **硬编码耦合**在 `PIPELINE_ORDER`（实测 L32）；若不先完成文件级 registry，SQLite 将把“中间产物”误建模为可变的 run 状态。  
- Phase 1C/1D **已不读 GA 文件**，但**依赖** GA 下游产物；先固化 bind+stale 文件契约，DB 仅持久化相同字段。  
- 1% 测试 **硬依赖** auto step4_5（实测）；先修集成再开 DB 可避免迁移期双套 pipeline。

---

## 附录 A：新语义 vs 当前实现差距清单

| 需求 | 当前实现 | 差距 |
|------|----------|------|
| GA 为 workspace 输入 | run 内 step4_5 生成 | 无 `workspace/ga_sets` |
| run 必须绑定 version | 无 binding | 无 `ga_binding.json` |
| 物化快照 | 即生成文件 | bind 复制未实现 |
| core 不自动 GA | `PIPELINE_ORDER` 含 step4_5 | **需改 registry** |
| 无绑定报错 | 仅 step6_7 缺 GA 时 warning→空 dict | **需 fail-fast** |
| step4_5 拆分 | 仅 `process_all` | **无 --mode** |

---

## 附录 B：`generate_ga_id` 与 schema（实测补充）

- **Deterministic**：同一 `GA_SMILES` 字符串 → 同一 `G{md5%1e6:06d}`。  
- **Collision**：不同 SMILES 可能同 ID（mod 1e6）；导入/version 创建需唯一性检查。  
- **Canonicalization**：提取时用 RDKit canonical SMILES；手工 CSV 需在 registry 再 canonicalize。  
- **非法 SMILES**：父分子无效则整分子跳过；`generate_ga_id` 失败则该行不写入。  

---

*本报告替代/补充 `_auditing/phase_1e_ga_workspace_plan.md` 中与“core 仍可默认 step4_5”相关的兼容假设；新语义以 **GA 必须先绑定、core 不自动生成** 为准。*
