# Phase 2 Readiness Assessment（GA Workspace 之后）

| 字段 | 值 |
|------|-----|
| 类型 | 审计 / Phase 2 就绪评估（**不修改代码、不实现 SQLite**） |
| 日期 | 2026-05-28 |
| 前置阶段 | Phase 1A–1E 已完成 |
| 依据 | 当前 `ChemDBWebVersion` 真实源码、`workspace/` 目录、`run_demo` 产物、pytest（49 passed） |
| 关联 | [phase_1e_implementation.md](./phase_1e_implementation.md)、[phase_1d_implementation.md](./phase_1d_implementation.md)、[phase_1abc_overall_assessment.md](./phase_1abc_overall_assessment.md) |

---

## 执行摘要

Phase 1A–1E 的 **file-based 系统已足够稳定，可以进入 Phase 2A（只读 ingest SQLite registry）**，但存在三类已知缺口需在 ingest 层处理，而非阻塞 Phase 2 启动：

1. **历史 run 与 Phase 1E 不同步**：`workspace/projects/p001/runs/run_demo` 为 Phase 1E 之前跑通，`pipeline_core.json` 仍含 `step4_5`，**无** `manifests/ga_binding.json`。
2. **无 project-level manifest**：`project_id` / `run_id` 仅能从目录名推断；`run.json` 不含显式 ID 字段。
3. **model / model-after 结果未关联 GA**：metrics、ranking 不含 `ga_set_id` / `ga_version_id`；需 ingest 时 join `run_ga_bindings`。

**建议**：Phase 2A 先做 **只读 ingest CLI + 临时 SQLite**，不改 orchestrator；legacy run 用路径解析 + 可选 `ga_binding` 缺失标记。

---

## 一、当前 workspace 总结构

### 1.1 关键目录树（审计时真实存在）

```text
workspace/
├── projects/
│   └── p001/
│       └── runs/
│           └── run_demo/                    # 唯一完整 E2E run（core@step13 + training@train）
│               ├── data/
│               │   ├── pubchem/*.csv        # 79 个 1% 文件
│               │   ├── metal_list.txt
│               │   ├── p_elements_list.txt
│               │   ├── metal_embedding/
│               │   │   ├── element_features.csv
│               │   │   └── element_features_zscore.csv
│               │   └── L3_embedding/
│               │       ├── L3_embedding_ecfp.npz
│               │       └── build.log
│               ├── tmp/                       # core 中间产物 + GA_with_id.csv（legacy 自动生成）
│               ├── training/
│               │   ├── index.json, config.yaml, split_index.json
│               │   ├── *_records.pkl
│               │   └── ckpts/{best.pt,last.pt,history.json}
│               ├── logs/*.log
│               ├── manifests/*.json           # 16 个 JSON（见第四节）
│               └── reports/                 # 空
├── ga_sets/
│   └── fixture_run_demo/
│       ├── ga_set.json
│       └── versions/v001/
│           ├── ga_version.json
│           └── GA_with_id.csv
├── model_after_tasks/                       # git 跟踪
│   ├── ni_lkb_p_ni/
│   │   ├── task.json, candidates.csv, positives.csv
│   └── pd_lkb_p_cluster/
│       ├── task.json, candidates.csv, positives.csv, negatives.csv
└── model_after_results/                     # .gitignore；本地存在
    ├── ni_lkb_p_ni/
    │   ├── p001_run_demo/                   # 单模型 evaluate-model
    │   └── batch_dummy_ckpts/               # 多模型 evaluate-models
    └── pd_lkb_p_cluster/
        └── p001_run_demo/
```

**说明**：

- `tests/run_isolation/fixtures/ga_sets/fixture_restructured/` 存在于测试 fixture，**未**写入 `workspace/ga_sets/`（仅 integration 测试临时安装）。
- `.gitignore`：`workspace/projects/**/runs/*`、`workspace/model_after_results/**` 被忽略；`ga_sets/`、`model_after_tasks/` 可进 git。

### 1.2 各类目录职责

| 目录 | 存什么 | 输入 / 产物 / 复用 / 一次性 |
|------|--------|------------------------------|
| `projects/{pid}/runs/{rid}/` | 一次 pipeline 运行的全部 run 级数据 | **混合**：`data/pubchem` 等为输入快照；`tmp/`、`training/`、`logs/`、`manifests/` 为运行产物 |
| `ga_sets/{ga_set_id}/` | GA 集合元数据 + 不可变 version | **可复用输入资源**（跨 run 共享）；version 目录创建后不可变 |
| `model_after_tasks/{task_id}/` | 评估任务定义（JSON + CSV） | **可复用输入资源**（workspace 级，跨 project） |
| `model_after_results/{task_id}/{batch_id}/` | 排序、指标、选模结果 | **一次性结果**（按 task + batch 归档） |

### 1.3 分类汇总

| 类别 | 路径示例 |
|------|----------|
| **输入资源** | `run/data/pubchem/`、`run/data/metal_list.txt`、`ga_sets/.../GA_with_id.csv`、`model_after_tasks/.../*.csv` |
| **运行产物** | `run/tmp/`、`run/training/`、`run/logs/`、`run/manifests/`、`model_after_results/` |
| **可复用资源** | `ga_sets/`、`model_after_tasks/`、已 bind 的 GA version（workspace 级） |
| **一次性结果** | 某次 run 的 `tmp/step13_*.csv`、某 batch 的 `ranking.csv` / `metrics.json` |

---

## 二、当前 Project / Run 抽象是否稳定

### 2.1 审计结论

| 问题 | 现状（源码 + 文件） |
|------|---------------------|
| `run_root` 如何定义 | `RunContext.from_run_root(path)` → `Path(run_root).resolve()`；子目录固定：`data/`, `tmp/`, `training/`, `logs/`, `reports/`, `manifests/`（`app/services/run_context.py`） |
| `project_id` / `run_id` 如何得到 | **仅从路径**：`workspace/projects/{project_id}/runs/{run_id}/`；代码中无独立解析函数；`model_after_paths.batch_id_from_run()` 从路径取 `{project_id}_{run_id}` |
| 是否有 `run.json` | **有**：`init-run` 写入 `manifests/run.json`（`orchestrator.init_run` L286） |
| 是否有 `project.json` | **无**（全仓库 0 个 `project.json`） |
| `init-run` 写了什么 | `run.json`：`action`, `run_root`, `chemdb_repo_root`, `created_at`, `seed_files_copied`, `pubchem_files`, `status`（可选 `metal_embedding_seed_error`） |
| run 状态如何判断 | 无统一 `run_status.json`；需推断：`run.json.status` + `pipeline_core.json.ok` + `pipeline_training.json.ok` + 各 step manifest 的 `status` |
| 是否记录输入数据 | **部分**：`run.json.pubchem_files`（文件名列表）；无 pubchem 源路径、无 MANIFEST 哈希 |
| 是否记录 pipeline 步骤 | **是**：每步 `manifests/{step_id}.json` + `pipeline_{core\|training}.json` |
| 是否记录 GA binding | **Phase 1E 设计有** `manifests/ga_binding.json`；**run_demo 无此文件**（legacy run） |
| 是否记录 training 产物 | **是**：`training_train.json` outputs 含 `best.pt`、`history.json`；`training/index.json` 含绝对路径 |
| 是否记录 model-after 产物 | **否**（不在 run manifest 内）；结果在 `workspace/model_after_results/`，通过 `model_run_root` 字符串关联 |

### 2.2 问答

1. **SQLite `projects` / `runs` 可从哪些文件 ingest？**
   - `projects`：目录 `workspace/projects/{project_id}/`（无 manifest → 仅 ID + 路径）
   - `runs`：`manifests/run.json` + `pipeline_core.json` + `pipeline_training.json` + 路径 `workspace/projects/{pid}/runs/{rid}/`

2. **是否缺 project-level manifest？**  
   **是**。当前无 `project.json`、无 project 创建时间/描述/owner。

3. **Phase 2A 是否需要补 `project.json` 或扩展 `run.json`？**  
   **建议补，但不阻塞 ingest**：
   - 最小：`run.json` 增加 `project_id`、`run_id`（init-run 写入，Phase 2A 可选补丁）
   - 可选：`projects/{pid}/project.json`（name、created_at、notes）
   - Phase 2A ingest 可 **先从路径解析 ID**，manifest 字段作为增强

4. **`run_id` 应来自目录名还是 manifest？**  
   **主键 = 目录名**（filesystem 契约）；manifest 可冗余存储便于搬迁校验。当前 manifest 无 `run_id` 字段。

5. **`project_id` 应来自目录名还是 manifest？**  
   **同上：目录名为主**；无 manifest 字段。

### 2.3 run_demo 特殊说明（legacy）

`run_demo/manifests/pipeline_core.json` 步骤列表为：`step1` → `step2` → **`step4_5`** → `step6_7` → … → `step13`（非 Phase 1E 的 `require_bound_ga` / `apply_ga_to_run`）。  
`tmp/GA_with_id.csv` 由 legacy `step4_5 --mode full-auto` 生成，**无** `ga_binding.json`、**无** workspace GA 溯源。

---

## 三、当前 GASet / GAVersion / RunGABinding 是否稳定

### 3.1 文件布局（Phase 1E 契约）

```text
workspace/ga_sets/{ga_set_id}/ga_set.json
workspace/ga_sets/{ga_set_id}/versions/{ga_version_id}/ga_version.json
workspace/ga_sets/{ga_set_id}/versions/{ga_version_id}/GA_with_id.csv
run_root/manifests/ga_binding.json
run_root/tmp/GA_with_id.csv                    # materialized snapshot
run_root/manifests/ga_stale_report.json        # 可选，rebind 时生成
```

### 3.2 Schema（来自 `ga_registry_manager.py` + 真实 fixture）

**`ga_set.json`**（`fixture_run_demo` 实测）：

| 字段 | 示例 / 说明 |
|------|-------------|
| `ga_set_id` | `"fixture_run_demo"` |
| `name` | `"run_demo GA (small)"` |
| `description` | string |
| `created_at`, `updated_at` | ISO8601 UTC |
| `tags` | `["fixture"]` |
| `source` | `{"type": "fixture"}` 或 `from_run` / `from_csv` |

**`ga_version.json`**（`v001` 实测）：

| 字段 | 说明 |
|------|------|
| `ga_set_id`, `ga_version_id` | 如 `fixture_run_demo`, `v001` |
| `parent_version_id` | null 或 `v001` |
| `created_at`, `created_by` | ISO8601；`ga_registry_manager` / `fixture` |
| `source` | `{"type": "fixture"}` 或 `from_run_generate_ga` / `from_csv` + 路径 |
| `num_ga` | `126` |
| `checksum` | `sha256:...`（canonical） |
| `notes` | string |
| `columns` | `["GA_SMILES", "GA_ID"]` |

**`GA_with_id.csv`**：

- 列：`GA_SMILES`, `GA_ID`（仅两列；materialize 时保证有 ID）
- checksum 算法：按 `GA_SMILES` 排序后 CSV 的 SHA256（`ga_csv_checksum()`）

**`ga_binding.json`**（`bind_ga_version_to_run()` L536–547，run_demo **无**，测试覆盖）：

| 字段 | 说明 |
|------|------|
| `run_root` | 绝对路径 |
| `ga_set_id`, `ga_version_id` | |
| `source_ga_csv` | workspace version CSV 绝对路径 |
| `materialized_path` | 相对 run：`"tmp/GA_with_id.csv"` |
| `checksum`, `num_ga` | 与 materialized 文件一致 |
| `bound_at`, `bound_by` | ISO8601；默认 `ga_registry_manager` |
| `status` | `"bound"` |

**`ga_stale_report.json`**（rebind 时，L341–357）：

| 字段 | 说明 |
|------|------|
| `generated_at`, `trigger` | |
| `previous_binding`, `new_binding` | 摘要 |
| `stale_files` | `[{path, exists, reason}]`；含 downstream tmp/training/ckpt 及 `model_after_results` 中匹配 `model_run_root` 的 manifest |
| `recommended_actions` | CLI 提示字符串 |

### 3.3 问答

| # | 问题 | 结论 |
|---|------|------|
| 1 | `ga_set.json` 字段 | 见上表 |
| 2 | `ga_version.json` 字段 | 见上表 |
| 3 | CSV schema | `GA_SMILES`, `GA_ID` |
| 4 | `ga_binding.json` 字段 | 见上表 |
| 5 | checksum 是否稳定 | **是**；canonical sort + sha256；bind 与 verify 一致（pytest 覆盖 mismatch） |
| 6 | `num_ga` 是否记录 | **是**（version + binding） |
| 7 | `source_ga_csv` 是否记录 | **是**（binding 绝对路径） |
| 8 | `materialized_path` 是否记录 | **是**（相对路径 `tmp/GA_with_id.csv`） |
| 9 | run 能否追踪 GA | **Phase 1E run 可以**（binding）；**run_demo 不能** |
| 10 | stale report 可否 ingest | **可以**；记录 rebind 事件、stale 路径列表、previous/new binding |
| 11 | SQLite ingest 来源 | `ga_sets/*/ga_set.json`；`versions/*/ga_version.json`；`run/manifests/ga_binding.json`；可选 `ga_stale_report.json` |
| 12 | 是否缺字段 | version 已有 `parent_version_id`, `source`, `notes`, `created_at`；可选补：`ga_csv_path`（可推导）、`content_type` |

### 3.4 GA version 不可变性

- `create_ga_version_from_csv` / `_write_ga_version_files`：version 目录 `exist_ok=False`（已存在则失败）
- pytest：`test_create_ga_version_from_csv_immutable` 验证 v001 不被覆盖
- 修改 workspace GA version **不会**改变已 bind run 的 snapshot（run 内 `tmp/GA_with_id.csv` 独立；checksum verify 可检测篡改）

---

## 四、当前 StepExecution / Artifact 是否稳定

### 4.1 run_demo 现有 manifest 列表

| step_id | 文件 | 备注 |
|---------|------|------|
| (init) | `run.json` | 非 step |
| step1 | `step1.json` | |
| step2 | `step2.json` | |
| step4_5 | `step4_5.json` | **legacy**；新 pipeline 应为 `apply_ga_to_run.json` |
| — | `require_bound_ga.json` | **run_demo 无** |
| — | `apply_ga_to_run.json` | **run_demo 无** |
| step6_7 | `step6_7.json` | |
| step8 | `step8.json` | |
| step12 | `step12.json` | |
| step13 | `step13.json` | |
| build_l3_embedding_ecfp | `build_l3_embedding_ecfp.json` | |
| zscore_metal_embedding | `zscore_metal_embedding.json` | |
| prepare_training_index | `prepare_training_index.json` | in-process |
| prepare_training_config | `prepare_training_config.json` | in-process |
| training_data | `training_data.json` | |
| training_train | `training_train.json` | |
| (summary) | `pipeline_core.json`, `pipeline_training.json` | |

### 4.2 Manifest schema（`orchestrator._finalize_manifest`）

**统一字段（subprocess 与 in-process 共有）**：

| 字段 | 存在 | 说明 |
|------|------|------|
| `step_id` | ✅ | |
| `status` | ✅ | `success` / `failed` / 初始 `running` |
| `started_at`, `finished_at` | ✅ | ISO8601 UTC |
| `duration_seconds` | ✅ | |
| `run_root` | ✅ | 绝对路径 |
| `cwd` | ✅ | 固定 `ChemDB` repo |
| `command` | ✅ | argv 或 `in_process:{step_id}` |
| `input_check` | ✅ | `{ok, missing, checked}` |
| `inputs` | ✅ | `[{path, exists, size_bytes}]` **绝对路径** |
| `outputs` | ✅ | 同上（finalize 时刷新） |
| `output_check` | ✅ | expected outputs 存在性 |
| `exit_code` | ✅ | subprocess；in-process 0/1 |
| `error` | ✅ | string 或 null |
| `log_file` | ✅ | subprocess 有；in-process 有 log 路径 |
| `paths_under_run_root` | ✅ | bool |
| `in_process_result` | 可选 | `require_bound_ga`, `prepare_training_*` |

**pipeline 汇总**（`pipeline_*.json`）：`action`, `pipeline`, `through`, `run_root`, `ok`, `steps[{step_id, status, manifest}]`, `finished_at`

### 4.3 问答

| # | 问题 | 结论 |
|---|------|------|
| 1 | schema 是否统一 | **基本统一**；in-process 多 `in_process_result`；`run.json` / `pipeline_*.json` 为不同 action |
| 2 | 关键字段齐全 | **是**（见上表） |
| 3 | inputs/outputs 绝对路径 | **是**（`_path_entry` → `path.resolve()`） |
| 4 | outputs 是否够登记 artifacts | **够** expected + mandatory；optional 默认不在 outputs 列表（除非 expand） |
| 5 | 能否判断 step 成功 | **是**：`status == "success"` 且 `output_check.ok` |
| 6 | 能否判断 artifact 属于 run | **是**：路径前缀 `run_root` + `paths_under_run_root: true` |
| 7 | 未记录的 important 产物 | `logs/*.log`；大量 `optional_outputs`（`*_stats.json`、neo4j 可选 CSV）；`data/pubchem` 输入文件；`model_after_results`；`reports/`；`last.pt` 为 optional 但在 manifest 若未跑 optional check 可能不在 outputs |
| 8 | SQLite ingest 方式 | 每 `manifests/{step_id}.json` → `step_executions`；`inputs`+`outputs` → `artifacts`；`pipeline_*.json` → run 级 summary 字段 |

---

## 五、当前 Model 抽象是否稳定

### 5.1 审计对象：`run_demo`（training 成功）

| 文件 | 作用 |
|------|------|
| `training/index.json` | run 内绝对路径索引 |
| `training/config.yaml` | 超参 + `embedding: ecfp` |
| `training/ckpts/best.pt` | ~227MB checkpoint |
| `training/ckpts/last.pt` | optional |
| `training/ckpts/history.json` | 每 epoch val metrics |
| `training/split_index.json` | seed + train/val/test target 列表 |
| `training/*_records.pkl` | 划分后训练样本 |
| `manifests/training_train.json` | step 执行记录 |

**`training/index.json` 关键字段**（实测）：`run_root`, `samples_csv`, `m_l3_pairs`, `l3_gac`, `kl_nl_ul_index.path`, `ligand_embedding.{backend,path,variants}`, `metal_embedding.path`

**`history.json` 字段**：`epoch`, `train_loss`, `val_loss`, `recall_at_1`, `recall_at_5`, `mrr`

### 5.2 问答

| # | 问题 | 结论 |
|---|------|------|
| 1 | model 如何定义 | **一次 run 内的一个 checkpoint**（`best.pt` 或 `last.pt` 等）；RankModel + 路径上下文 |
| 2 | `model_id` 是否存在 | **run 级不存在**；model-after 侧默认 `Path(run_root).name` 或 CLI `--model-id` |
| 3 | Phase 2A 如何生成 `model_id` | 建议：`{project_id}_{run_id}_{checkpoint_stem}` 如 `p001_run_demo_best`；或 UUID |
| 4 | model 与 run 关系 | **一对多**（同一 run 可有 best + last + 未来多 ckpt） |
| 5 | `checkpoint_path` | `training/ckpts/best.pt`（或 manifest outputs / model_runs.csv） |
| 6 | `config_path` | `training/config.yaml` |
| 7 | `index_path` | `training/index.json` |
| 8 | `embedding_backend` | `config.yaml` → `data.embedding`；`index.json` → `ligand_embedding.backend` |
| 9 | `l3_embedding_path` | `index.json` → `ligand_embedding.path` |
| 10 | `metal_embedding_path` | `index.json` → `metal_embedding.path` |
| 11 | training metrics | `training/ckpts/history.json` 最后 epoch；无单独 `metrics.json` |
| 12 | 是否应记录 GA | **应记录**（复现必需）；**当前 model 文件无 GA 字段** → ingest 时 join `run_ga_bindings` |
| 13 | SQLite ingest 来源 | `training_train.json` + `training/index.json` + `config.yaml` + `history.json` + checkpoint 文件 stat |
| 14 | 是否缺 `model.json` | **缺**；Phase 2A 可由 ingest 合成，或 Phase 2A 可选在 training 完成后写 `training/model.json`（非必须） |

---

## 六、当前 model-after task 是否稳定

### 6.1 `task.json` schema（实测）

```json
{
  "task_id": "ni_lkb_p_ni",
  "metal": "Ni",
  "embedding": "ecfp",
  "candidate_file": "candidates.csv",
  "positive_file": "positives.csv",
  "negative_file": "" | "negatives.csv",
  "notes": "..."
}
```

### 6.2 CSV schema

- `candidates.csv`：`candidate_id, smiles, did, name, source`
- `positives.csv`：`positive_id, ...`（同列集）
- `negatives.csv`：`negative_id, ...`（可选）

`load_task_bundle()`（`ChemDB/training/model_after/io.py`）运行时计算行数，**不写入 task.json**。

### 6.3 问答

| # | 问题 | 结论 |
|---|------|------|
| 1 | task.json schema | 见上 |
| 2 | `task_id` 如何确定 | `task.json.task_id` 或 **目录名** |
| 3 | metal | `task.json.metal` |
| 4 | CSV 路径 | 相对 task 目录的文件名；绝对路径 = `task_dir / candidate_file` |
| 5 | count 是否记录 | **否**（ingest 时可读 CSV 计算） |
| 6 | 是否独立于 project | **是**（workspace 级；`model_after_paths.MODEL_AFTER_TASKS_ROOT`） |
| 7 | tags/notes/source | 有 `notes`；无 `tags`/`source` 字段（可加 optional） |
| 8 | SQLite ingest | `model_after_tasks/{task_id}/task.json` + CSV 路径 + 可选 row counts |
| 9 | 是否缺 task manifest | 无独立 `task_manifest.json`；**task.json 足够** |

---

## 七、当前 model-after result 是否稳定

### 7.1 目录模式（实测两种）

**A. 单模型 `evaluate-model`** → `model_after_results/{task_id}/{project_id}_{run_id}/`（默认加 `_eval` 逻辑在 `default_output_dir`；实测 `p001_run_demo` **无** `_eval` 后缀，因显式 `--output-dir`）

**B. 多模型 `evaluate-models`** → `model_after_results/{task_id}/{batch_id}/`（如 `batch_dummy_ckpts`）

### 7.2 文件清单

| 文件 | 单模型 batch | 多模型 batch |
|------|-------------|-------------|
| `ranking.csv` | ✅ 根目录 | ✅ `models/{model_id}/` |
| `metrics.json` | ✅ | ✅ |
| `report.json` | ✅ | ✅ |
| `heldout_ranking.csv` | ✅ | ✅ |
| `evaluation_manifest.json` | ✅（小） | ✅（在 models/ 下） |
| `evaluate_model_manifest.json` | ✅ | — |
| `evaluate_models_manifest.json` | — | ✅ |
| `model_selection_summary.csv` | — | ✅ |
| `best_model.json` | — | ✅ |
| `model_selection_manifest.json` | — | ✅ |
| `model_runs.csv` | — | ✅（输入副本/引用） |

### 7.3 关键字段

**`metrics.json`**：`model_id`, `task_id`, `metal`, `num_candidates`, `num_positives`, `num_negatives`, `mrr`, `hit_at_5/10/20/50`, `positive_score_mean`, `negative_score_mean`, `pos_neg_margin`, `status`, `warnings`, `context_source`, `model_run_root`, `checkpoint_path`, `ranking_path`

**`ranking.csv`**：`rank, score, candidate_id, did, smiles, name, source, is_known_positive, is_known_negative, valid, error, model_id, model_run_root, checkpoint_path`

**`best_model.json`**：`model_id`, `model_run_root`, `checkpoint_path`, `selection_score`, `selection_reason`, `metrics_path`, `ranking_path`

**`model_selection_summary.csv`**：`model_id, model_name, model_run_root, checkpoint_path, embedding_backend, status, mrr, hit_at_5, ..., selection_score, rank_among_models, error`

**`batch_id` 确定**：CLI `--batch-id`；默认 `batch_default`（多模型）或 `batch_id_from_run()` → `{project_id}_{run_id}`

### 7.4 问答

| # | 问题 | 结论 |
|---|------|------|
| 1 | batch_id | **目录名** 或 CLI 参数 |
| 2 | batch 是否记录 task_id | **间接**：路径 `{task_id}/{batch_id}`；manifest 含 `task_dir` 可解析 |
| 3 | batch 是否记录 model_runs_file | **evaluate_models_manifest.json** 有 |
| 4 | batch 是否记录 output_dir | **是**（各 manifest） |
| 5–8 | CSV/JSON 字段 | 见上 |
| 9 | 是否记录 GA | **否**；需 join run → `ga_binding` |
| 10 | SQLite ingest | `evaluate_*_manifest.json` + `metrics.json` + `best_model.json` + `model_selection_summary.csv` |
| 11 | Streamlit 缺字段 | 缺 `project_id`, `run_id`（可从 `model_run_root` 解析）；缺 `ga_*`；有 `model_name`, `task`/`metal`, `selection_score` |

---

## 八、当前 no-write 和 reproducibility 状态

| # | 检查项 | 结论 |
|---|--------|------|
| 1 | 不写 `ChemDB/tmp` | **pytest autouse** `assert_repo_tmp_unchanged`（`tests/run_isolation/conftest.py`） |
| 2 | 不写 `ChemDB/src/tmp` | **同上** |
| 3 | 不写 `ChemDB/training/evaluation/results` | **REPO_ARTIFACT_PATTERNS** 守护（Phase 1C test）；model-after 不写该路径 |
| 4 | GA registry 不写 ChemDB 源码 | **是**；仅写 `workspace/ga_sets` + run `manifests/` + run `tmp/` |
| 5 | model-after 不写 task_dir | **是**（pytest + `RESULT_FILENAMES` 设计） |
| 6 | model-after results 只写 output_dir | **是** |
| 7 | GA version immutable | **是**（见第三节） |
| 8 | `run/tmp/GA_with_id.csv` 是 snapshot | **是**（materialize 复制 + 独立 checksum） |
| 9 | 修改 workspace GA 不影响已绑定 run | **是**（除非 rebind） |
| 10 | absolute path 搬迁风险 | **有**：`run.json`、`index.json`、`config.yaml`、manifest、`metrics.json` 均含绝对路径；**搬迁 workspace 需 re-ingest 或路径重写** |

---

## 九、Phase 2 SQLite schema 初步建议

> 原则：DB 存路径 + 摘要 + 外键；不存大文件内容。

### 9.1 `projects`

| 项 | 内容 |
|----|------|
| PK | `project_id` TEXT |
| 字段 | `project_root`, `name`, `created_at`, `updated_at`, `notes` |
| 来源 | 目录 `workspace/projects/{id}/`；`project.json`（Phase 2A 可选） |
| 2A 补充 | `name`, `created_at` 可先空或默认 = project_id |
| 可空 | `notes`, `description` |

### 9.2 `runs`

| 项 | 内容 |
|----|------|
| PK | `run_id` TEXT |
| FK | `project_id` → projects |
| 字段 | `run_root`, `created_at`, `chemdb_repo_root`, `status`, `pubchem_file_count`, `pipeline_core_ok`, `pipeline_training_ok`, `last_finished_at` |
| 来源 | `manifests/run.json`, `pipeline_core.json`, `pipeline_training.json` |
| 2A 补充 | `project_id`/`run_id` 从路径解析写入 |
| 可空 | `pipeline_training_ok`（未跑 training） |

### 9.3 `ga_sets`

| 项 | 内容 |
|----|------|
| PK | `ga_set_id` TEXT |
| 字段 | `name`, `description`, `created_at`, `updated_at`, `tags_json`, `source_json` |
| 来源 | `ga_sets/{id}/ga_set.json` |

### 9.4 `ga_versions`

| 项 | 内容 |
|----|------|
| PK | `(ga_set_id, ga_version_id)` |
| FK | `ga_set_id` → ga_sets |
| 字段 | `parent_version_id`, `created_at`, `created_by`, `num_ga`, `checksum`, `ga_csv_path`, `source_json`, `notes` |
| 来源 | `versions/{vid}/ga_version.json` + CSV path |

### 9.5 `run_ga_bindings`

| 项 | 内容 |
|----|------|
| PK | `binding_id` INTEGER 或 `(run_id, bound_at)` |
| FK | `run_id`, `ga_set_id`, `ga_version_id` |
| 字段 | `checksum`, `num_ga`, `source_ga_csv`, `materialized_path`, `bound_at`, `bound_by`, `status`, `is_current` |
| 来源 | `manifests/ga_binding.json`；history 来自 `ga_stale_report.json.previous_binding` |
| 2A 补充 | legacy run 无 binding → `status=unknown` 或 NULL |
| 可空 | 全表对 legacy run |

### 9.6 `step_executions`

| 项 | 内容 |
|----|------|
| PK | `(run_id, step_id, started_at)` 或 surrogate `execution_id` |
| FK | `run_id` → runs |
| 字段 | `step_id`, `status`, `command_json`, `cwd`, `started_at`, `finished_at`, `duration_seconds`, `exit_code`, `error`, `log_path`, `manifest_path`, `in_process_result_json` |
| 来源 | `manifests/{step_id}.json` |

### 9.7 `artifacts`

| 项 | 内容 |
|----|------|
| PK | `artifact_id` |
| FK | `run_id`, optional `execution_id` |
| 字段 | `rel_path`, `abs_path`, `role` (input/output), `step_id`, `size_bytes`, `exists`, `artifact_type` |
| 来源 | manifest `inputs`/`outputs`；可扩展扫描 `training/ckpts/*.pt` |
| 2A 补充 | `artifact_type` 枚举（csv/json/npz/pt/pkl/log） |

### 9.8 `models`

| 项 | 内容 |
|----|------|
| PK | `model_id` TEXT |
| FK | `run_id` |
| 字段 | `checkpoint_path`, `config_path`, `index_path`, `embedding_backend`, `l3_embedding_path`, `metal_embedding_path`, `history_path`, `best_epoch`, `val_mrr`, `val_recall_at_1`, `val_recall_at_5`, `ga_set_id`, `ga_version_id`, `ga_binding_checksum` |
| 来源 | `training/*` + `training_train.json` + join `run_ga_bindings` |
| 2A 补充 | `model_id` 合成；GA 字段 join |
| 可空 | GA 字段（legacy） |

### 9.9 `model_after_tasks`

| 项 | 内容 |
|----|------|
| PK | `task_id` TEXT |
| 字段 | `task_dir`, `metal`, `embedding`, `candidate_file`, `positive_file`, `negative_file`, `notes`, `candidate_count`, `positive_count`, `negative_count` |
| 来源 | `task.json` + CSV 行数 |

### 9.10 `model_after_batches`

| 项 | 内容 |
|----|------|
| PK | `(task_id, batch_id)` |
| FK | `task_id` |
| 字段 | `output_dir`, `command`, `model_runs_file`, `finished_at`, `n_models`, `n_success`, `best_model_id`, `summary_path` |
| 来源 | `evaluate_models_manifest.json` / 目录扫描 |

### 9.11 `model_after_model_results`

| 项 | 内容 |
|----|------|
| PK | `(task_id, batch_id, model_id)` |
| FK | batch, optional `model_id` → models |
| 字段 | `model_run_root`, `checkpoint_path`, `project_id`, `run_id`, `mrr`, `hit_at_*`, `selection_score`, `rank_among_models`, `status`, `ranking_path`, `metrics_path`, `report_path` |
| 来源 | `metrics.json`, `model_selection_summary.csv`, `best_model.json` |
| 2A 补充 | `project_id`/`run_id` 从 `model_run_root` 解析 |

---

## 十、Phase 2A 最小实施路线

### 10.1 原则（用户给定）

只做 SQLite registry；不改 pipeline；orchestrator 暂不写 DB；filesystem ingest；DB 不存大文件；测试用 tmp_path；不污染真实 workspace；无 Streamlit。

### 10.2 推荐模块

```text
app/services/registry_db/
  schema.sql          # DDL
  ingest.py           # 纯函数：path → rows
  cli.py              # init-db, register-*, rebuild-index
  connection.py       # sqlite path 管理
tests/registry_db/    # tmp workspace + tmp chemdb.sqlite
```

### 10.3 最小 CLI 设计

| CLI | 读取 | 写入表 | 说明 |
|-----|------|--------|------|
| `init-db` | `schema.sql` | 全部 DDL | `--db-path` 默认 `{workspace}/chemdb.sqlite` 或测试临时路径 |
| `register-project` | `workspace/projects/{pid}/` [+ 可选 `project.json`] | `projects` | 无 manifest 则仅 ID + root |
| `register-run` | `runs/{rid}/manifests/{run,pipeline_*}.json` | `runs`, `step_executions` | 可选 `--include-steps` |
| `register-ga-set` | `ga_sets/{id}/ga_set.json`, `versions/*/ga_version.json` | `ga_sets`, `ga_versions` | 可 `--ga-set-id` 或 scan 全部 |
| `register-run-ga-binding` | `run/.../ga_binding.json` [+ stale] | `run_ga_bindings` | 缺失则 skip 或 `status=missing` |
| `register-task` | `model_after_tasks/{id}/task.json` + CSV | `model_after_tasks` | 计算 counts |
| `register-model` | `run/training/` + `training_train.json` + binding join | `models`, `artifacts` | 每个 `.pt` 一行 |
| `register-model-after-batch` | `model_after_results/{task}/{batch}/` | `model_after_batches`, `model_after_model_results` | scan manifest + summary |
| `rebuild-index` | 上述全部 | UPSERT 全表 | `--workspace-root`；幂等 |

**幂等键建议**：

- `runs`: `(project_id, run_id)`
- `ga_versions`: `(ga_set_id, ga_version_id)`
- `step_executions`: `(run_id, step_id, started_at)` 或 manifest mtime + step_id
- `models`: `(run_id, checkpoint_path)`
- `model_after_model_results`: `(task_id, batch_id, model_id)`

### 10.4 实施顺序

1. DDL + `init-db`
2. ingest 函数（单表单元测试）
3. `register-run` + `register-ga-set`（最高价值）
4. `register-model` + binding join
5. model-after ingest
6. `rebuild-index` 聚合

---

## 十一、Phase 2A 测试建议

| # | 测试 | 断言 |
|---|------|------|
| 1 | 初始化临时 SQLite | `tmp_path/chemdb.sqlite` 存在；表齐全 |
| 2 | 注册 project/run | `runs` 行数；`run_root` 正确 |
| 3 | 注册 GA set/version | checksum/num_ga 匹配 fixture |
| 4 | 注册 run GA binding | bind 后 ingest；checksum 一致 |
| 5 | 注册 step_executions/artifacts | 与 `step1.json` inputs/outputs 数量一致 |
| 6 | 注册 model | `best.pt` 路径；history val_mrr |
| 7 | 注册 workspace task | counts 与 CSV 一致 |
| 8 | 注册 model-after batch | summary 行数 = n_models |
| 9 | 查询 task 下所有 results | SQL join batches + results |
| 10 | 查询 best model | `best_model.json` → `model_id` 匹配 |
| 11 | 重复 ingest 不重复 | UPSERT；行数不变 |
| 12 | 文件缺失 | `status=missing` / `ingest_error` 列；不抛 uncaught |
| 13 | 不写真实 workspace | 全在 `tmp_path/workspace` |
| 14 | 不写 ChemDB 源码 | 复用 `assert_repo_tmp_unchanged` |

**Fixture 策略**：复制 `tests/run_isolation/fixtures/` + 最小 manifest 集合；GA 用 `fixture_run_demo`；不必依赖 18M pubchem。

---

## 十二、最终判断

### 12.1 是否可以进入 Phase 2A？

**可以。** Phase 1A–1E file-based 契约清晰，manifest 字段足够支撑 **只读 ingest**；pytest 49 passed 表明隔离与 GA 逻辑稳定。

### 12.2 若不能，必须先补哪些 file-based metadata？

不阻塞 2A，但建议在 **Phase 2A 并行或紧随其后** 补充：

| 优先级 | 项 | 原因 |
|--------|-----|------|
| P1 | 新 run 必须有 `ga_binding.json` | Phase 1E 已强制；legacy run 需标记 |
| P2 | `run.json` 增加 `project_id`, `run_id` | 减少纯路径解析 |
| P3 | `training/model.json` 或 `manifests/model_best.json` | 显式 model 登记（ingest 可合成则非必须） |
| P4 | `project.json` | project 元数据 |
| P5 | model-after `metrics.json` 增加 `ga_set_id`/`ga_version_id` | 跨 join 简化（ingest 时可 join 替代） |

### 12.3 最小必须补的 manifest

**无硬阻塞。** ingest 层必须处理：

- 无 `ga_binding.json` 的 legacy run
- 无 `project.json`
- 无 `model_id` 的 run 级 model

### 12.4 哪些表可以第一版先不做？

| 表 | 建议 |
|----|------|
| 可先简化 | `artifacts` 可先只存 key artifacts（ckpt/npz/index） |
| 可延后 | stale binding 历史（单独 `ga_binding_events`） |
| 可延后 | pubchem 输入文件级 artifact（仅 count） |

### 12.5 Phase 2A 最容易出错的地方

1. **绝对路径** vs 搬迁后的 workspace root  
2. **legacy run_demo** 与 Phase 1E pipeline 步骤 ID 不一致（`step4_5` vs `apply_ga_to_run`）  
3. **model_after 两种 batch 布局**（单模型根目录 vs `models/{id}/`）  
4. **幂等 UPSERT** 键选择不当导致重复行  
5. **一对多 checkpoint**（best + last）与 `models` 行  
6. **ingest 误扫** `ChemDB/tmp` 或 repo 内 training 产物  

### 12.6 是否建议先只读 ingest，再做 `--write-db` 集成？

**强烈建议。**

- **Phase 2A**：只读 ingest CLI + `rebuild-index`；orchestrator / ga_registry / model_after **不改**  
- **Phase 2B**（可选）：orchestrator 每步完成后 optional hook 写 DB  
- 好处：验证 schema 与查询 API，不引入 pipeline 回归风险  

### 12.7 Phase 3 Streamlit 主要读取哪些表？

| UI 能力 | 主要表 |
|---------|--------|
| Project / Run 列表 | `projects`, `runs`, `step_executions` |
| Run 详情 / 产物 | `artifacts`, `step_executions` |
| GA 选择与绑定 | `ga_sets`, `ga_versions`, `run_ga_bindings` |
| 训练模型 | `models` |
| 评估任务 | `model_after_tasks` |
| 评估结果 / 选模 | `model_after_batches`, `model_after_model_results` |

Streamlit **不应**直接扫 filesystem；应通过 SQLite 查路径后再读 CSV/JSON/加载 checkpoint。

---

## 附录 A：Phase 1E 新 core pipeline vs run_demo

| 项 | Phase 1E 新 run | run_demo（legacy） |
|----|-----------------|-------------------|
| Core 步骤 | step1→2→require_bound_ga→apply_ga→6_7→8→12→13 | step1→2→**step4_5**→6_7→… |
| GA 来源 | workspace bind | step4_5 自动生成 |
| `ga_binding.json` | 有 | **无** |
| SQLite ingest | 完整 GA 溯源 | GA 字段 NULL + 标记 `legacy_auto_ga` |

## 附录 B：pytest 状态（审计日）

```text
49/61 tests collected (12 deselected)
默认：49 passed（~1s）
```

## 附录 C：修订记录

| 版本 | 日期 | 说明 |
|------|------|------|
| 1.0 | 2026-05-28 | Phase 2 readiness 初版（post Phase 1E） |
