# ChemDB Web Phase 3：Resource-Centric 前端设计草案

> **Phase 3A 就绪审计**（只读 Streamlit dashboard 范围、repository 缺口、实施计划）：见 [phase_3a_readiness_audit.md](./phase_3a_readiness_audit.md)。本文档描述完整 Phase 3 愿景；Phase 3A 为其只读展示子集。

## 1. 总体定位

Phase 3 的目标不是把已有命令行步骤简单搬到网页上，而是把 ChemDB 包装成一个面向非程序用户的操作系统。

用户不应该直接面对 `step1.py`、`step2.py`、`step4_5.py`、`step13.py`、`training.train` 等底层脚本，而应该面对更自然的研究对象：

```text
GA 列表
Task 列表
模型库
模型构建
模型选择
推荐结果
```

因此，Phase 3 的前端设计应采用 **resource-centric** 结构，而不是 run-centric 或 step-centric 结构。

后台仍然保留：

```text
Project
Run
StepExecution
Artifact
Model
Task
GASet
GAVersion
```

但这些主要作为系统内部记录和调试信息存在。用户主流程应围绕资源和操作展开。

---

## 2. 当前已完成基础

目前系统已经完成以下基础设施：

```text
Phase 1A–1E：
ChemDB 已基本改造成 run-based 文件系统。
每个 run 有独立 run_root。
GA 列表已经从 run 中间产物改为 workspace-level ga_sets 资源。
run 必须绑定某个 GA version 后才能继续下游流程。

Phase 2A：
SQLite read-only registry 已完成。
可以从 workspace ingest 项目、run、GA、模型、task、评估结果等元数据。
legacy_auto_ga run 和 bound GA run 均已验证。
model-after result 可以通过 registry_model_id 关联回 models。
```

当前可用的核心对象包括：

```text
workspace/projects/
workspace/ga_sets/
workspace/model_after_tasks/
workspace/model_after_results/
workspace/chemdb.sqlite
```

---

## 3. 前端核心抽象

最终前端应围绕以下 6 个主要对象组织：

```text
1. GA Library
2. Task Library
3. Model Builder
4. Model Library
5. Model Selector
6. Recommendation
```

其中：

```text
GA Library：
管理和注入 GA 列表。

Task Library：
管理和注入模型后评估任务。

Model Builder：
选择数据集和 GA 列表，自动运行模型前流程，产出模型。

Model Library：
展示已经训练好的模型。

Model Selector：
选择一组 task 和一组模型，自动评估并排序。

Recommendation：
使用选定模型对候选配体进行真实推荐。
```

超参数库暂时推迟。原因是当前部分训练参数可能仍然写死在 `src` 或 `training` 中，后续需要在整体设计稳定后，再统一抽象为 workspace-level hyperparameter resource。

---

## 4. 最终用户工作流

用户视角下，完整流程应是：

```text
1. 注入或选择 GA 列表
2. 注入或选择 Task
3. 选择数据集
4. 使用选定 GA 列表构建模型
5. 模型自动加入模型库
6. 选择任务范围和模型范围
7. 系统自动评估模型并排序
8. 选择最合适模型
9. 使用该模型进行真实推荐
```

对应关系为：

```text
GA List + Dataset
        ↓
Model Builder
        ↓
Model Library

Task List + Model Scope
        ↓
Model Selector
        ↓
Best Model / Ranking

Best Model + Candidate List
        ↓
Recommendation
```

---

## 5. 页面结构建议

左侧导航栏建议为：

```text
1. Dashboard
2. GA Library
3. Task Library
4. Model Builder
5. Model Library
6. Model Selector
7. Recommendation
8. Jobs & Logs
9. Settings
```

其中 Phase 3A 可以先做只读或半只读版本：

```text
1. Dashboard
2. GA Library
3. Task Library
4. Model Library
5. Model Selector Results
6. Settings
```

Phase 3B 再加入操作按钮：

```text
Upload GA
Upload Task
Start Model Build
Run Model Selection
Run Recommendation
```

---

# 6. Dashboard 页面

## 6.1 页面目的

Dashboard 是首页，用于回答：

```text
当前系统里有什么？
有多少 GA 列表？
有多少任务？
有多少模型？
最近有哪些评估结果？
```

## 6.2 页面内容

建议显示统计卡片：

```text
Projects
Runs
GA Sets
GA Versions
Models
Tasks
Evaluation Batches
Evaluation Results
```

示意：

```text
ChemDB Dashboard

┌────────────┐ ┌────────────┐ ┌────────────┐ ┌────────────┐
│ Projects   │ │ Runs       │ │ Models     │ │ GA Sets    │
│ 1          │ │ 1          │ │ 2          │ │ 1          │
└────────────┘ └────────────┘ └────────────┘ └────────────┘

┌────────────┐ ┌────────────┐ ┌────────────┐
│ Tasks      │ │ Batches    │ │ Results    │
│ 2          │ │ 3          │ │ 4          │
└────────────┘ └────────────┘ └────────────┘
```

下方显示：

```text
Recent Models
Recent Tasks
Recent Evaluation Results
```

---

# 7. GA Library 页面

## 7.1 页面目的

GA Library 负责管理 workspace-level GA resources。

它回答：

```text
系统中有哪些 GA 列表？
每个 GA 列表有哪些版本？
每个版本包含多少 GA？
用户能否上传、生成、查看、下载、编辑 GA？
```

## 7.2 页面布局

```text
GA Library

[Upload GA CSV] [Generate GA from Run] [Create New Version]

GA Sets
┌──────────────────┬────────────┬──────────┬──────────────┬────────────┐
│ GA Set           │ Version    │ Num GA   │ Created At   │ Notes      │
├──────────────────┼────────────┼──────────┼──────────────┼────────────┤
│ aromatic_basic   │ v001       │ 126      │ 2026-05-21   │ default    │
│ aromatic_basic   │ v002       │ 141      │ 2026-05-28   │ edited     │
│ crown_test       │ v001       │ 165      │ 2026-05-28   │ crown GA   │
└──────────────────┴────────────┴──────────┴──────────────┴────────────┘
```

点开某个 GA version 后显示：

```text
GA Version Detail

GA Set: aromatic_basic
Version: v002
Num GA: 141
Checksum: sha256:...
Created At: ...

GA Preview
┌───────────────┬─────────┬────────┬──────────────┐
│ GA_SMILES     │ GA_ID   │ Active │ Note         │
├───────────────┼─────────┼────────┼──────────────┤
│ c1ccccc1      │ G874135 │ yes    │ benzene ring │
│ c1ccncc1      │ G289719 │ yes    │ pyridine     │
└───────────────┴─────────┴────────┴──────────────┘

[Download CSV]
[Create Edited Version]
```

## 7.3 GA 注入方式

GA 列表应支持三种注入方式：

```text
1. 上传 GA_with_id.csv
2. 从某个 run 或 ligand dataset 自动生成 GA candidates
3. 基于已有 GA version 编辑后保存为新 version
```

用户不需要知道 `step4_5.py`。页面上应写成：

```text
Generate GA from ligand data
```

后台再调用：

```text
ga_registry_manager generate-ga-from-run
ga_registry_manager create-ga-version-from-csv
ga_registry_manager bind-ga-version-to-run
```

---

# 8. Task Library 页面

## 8.1 页面目的

Task Library 负责管理模型后评估任务。

Task 是模型后流程的输入资源，通常包含：

```text
candidate list
positive list
negative list，可选
metal
embedding type
notes
```

## 8.2 页面布局

```text
Task Library

[Create Task] [Upload Task CSV] [Import Legacy Quest]

Tasks
┌──────────────────┬───────┬────────────┬───────────┬───────────┬──────────┐
│ Task ID          │ Metal │ Candidates │ Positives │ Negatives │ Notes    │
├──────────────────┼───────┼────────────┼───────────┼───────────┼──────────┤
│ ni_lkb_p_ni      │ Ni    │ 344        │ 11        │ 0         │ Ni test  │
│ pd_lkb_p_cluster │ Pd    │ 343        │ 10        │ 1         │ Pd test  │
└──────────────────┴───────┴────────────┴───────────┴───────────┴──────────┘
```

创建 task 时：

```text
Create Task

Task ID:        [ ni_lkb_p_ni ]
Metal:          [ Ni ▼ ]
Embedding:      [ ecfp ▼ ]

Candidates CSV: [ Upload ]
Positives CSV:  [ Upload ]
Negatives CSV:  [ Upload optional ]

[Validate Task]
[Save Task]
```

点开 task 后：

```text
Task Detail

Tabs:
[Overview] [Candidates] [Positives] [Negatives] [Previous Evaluations]
```

## 8.3 Task 存储位置

后台仍使用：

```text
workspace/model_after_tasks/{task_id}/
    task.json
    candidates.csv
    positives.csv
    negatives.csv
```

用户不需要知道实际路径。

---

# 9. Model Builder 页面

## 9.1 页面目的

Model Builder 是模型前流程的主入口。

用户在这里完成：

```text
选择数据集
选择 GA 列表
启动自动模型构建
```

当前暂时不加入独立超参数库，超参数仍沿用现有默认配置或内部配置。等整体流程稳定后，再加入 Hyperparameter Library。

## 9.2 页面布局

```text
Model Builder

1. Select Input Data
Project:       [ p001 ▼ ]
Dataset:       [ PubChem 1% demo ▼ ]

2. Select GA List
GA Set:        [ aromatic_basic ▼ ]
GA Version:    [ v002 ▼ ]

3. Build Options
Run ID:        [ run_20260528_001 ]
Embedding:     [ ecfp ▼ ]
Pipeline:      [ Full model build ▼ ]

[Start Build]
```

## 9.3 用户点击 Start Build 后的后台流程

前端只显示“构建模型”，后台实际执行：

```text
init-run
bind GA version
step1
step2
apply-ga
step6_7
step8
step12
step13
build embedding
prepare training index
prepare training config
training data
training train
rebuild SQLite index
```

用户不需要看到这些内部步骤，除非进入 Jobs & Logs 页面。

## 9.4 进度显示

```text
Building Model: run_20260528_001

Current stage:
[██████████████░░░░░░] 70%

Stages:
✓ Initialize run
✓ Bind GA list
✓ Clean and repair data
✓ Apply GA and generate IRL
✓ Build training samples
✓ Build embeddings
→ Train model
○ Register model
```

## 9.5 完成后显示

```text
Build Completed

Created models:
┌──────────────────────────────┬────────────┬────────────┬────────────┐
│ Model ID                     │ Checkpoint │ GA Version │ Val MRR    │
├──────────────────────────────┼────────────┼────────────┼────────────┤
│ p001_run_20260528_001_best   │ best       │ v002       │ 0.083      │
│ p001_run_20260528_001_last   │ last       │ v002       │ 0.078      │
└──────────────────────────────┴────────────┴────────────┴────────────┘

[Go to Model Library]
[Evaluate on Tasks]
```

---

# 10. Model Library 页面

## 10.1 页面目的

Model Library 展示系统已有模型。

模型后流程主要从这里选择模型。用户不需要把 GA、run、训练细节放在第一优先级，但可以在详情里查看。

## 10.2 页面布局

```text
Model Library

Filters:
Project [all ▼]  Model status [available ▼]  Embedding [ecfp ▼]

┌──────────────────────────────┬──────────┬────────────┬──────────┬──────────┐
│ Model ID                     │ Source   │ Checkpoint │ Val MRR  │ Status   │
├──────────────────────────────┼──────────┼────────────┼──────────┼──────────┤
│ p001_run_demo_best           │ run_demo │ best       │ 0.076    │ ready    │
│ p001_run_demo_last           │ run_demo │ last       │ 0.076    │ ready    │
│ p001_run_20260528_001_best   │ run_001  │ best       │ 0.083    │ ready    │
└──────────────────────────────┴──────────┴────────────┴──────────┴──────────┘
```

点开模型后显示：

```text
Model Detail

Model:
- model_id
- checkpoint
- val_mrr
- embedding_backend

Build Context:
- input dataset
- GA set / version
- run_root
- config_path
- index_path
```

如果后续检测到 stale model，应显示：

```text
Warning: This model may be stale because GA binding changed after training.
```

---

# 11. Model Selector 页面

## 11.1 页面目的

Model Selector 是模型后流程的核心页面。

它不要求用户关心每个模型的构建细节，而是让用户选择：

```text
一组任务
一组模型范围
一个排序指标
```

系统自动评估模型并排序。

## 11.2 页面布局

```text
Model Selector

1. Select Tasks
[✓] ni_lkb_p_ni
[✓] pd_lkb_p_cluster
[ ] other_task

2. Select Model Scope
Model source:
( ) All available models
( ) Models from selected project
( ) Models from selected runs
( ) Manually selected models

Project: [ p001 ▼ ]
Runs:    [ run_demo, run_20260528_001 ▼ ]

3. Selection Metric
Primary metric: [ MRR ▼ ]
Tie-breaker:    [ Hit@5 ▼ ]

[Run Model Selection]
```

## 11.3 后台逻辑

用户点击 `Run Model Selection` 后，后台执行：

```text
for task in selected_tasks:
    for model in selected_models:
        evaluate model on task
    rank models by selected metric
save model_after_results
rebuild SQLite index
```

## 11.4 结果页面

```text
Model Selection Results

Task: ni_lkb_p_ni
┌──────┬──────────────────────────────┬──────┬───────┬────────┬──────────┐
│ Rank │ Model                        │ MRR  │ Hit@5 │ Hit@10 │ Select   │
├──────┼──────────────────────────────┼──────┼───────┼────────┼──────────┤
│ 1    │ p001_run_demo_last           │ .171 │ .187  │ .273   │ Use      │
│ 2    │ p001_run_20260528_001_best   │ .083 │ .091  │ .143   │ Use      │
│ 3    │ p001_run_demo_best           │ .016 │ 0.000 │ 0.000  │ Use      │
└──────┴──────────────────────────────┴──────┴───────┴────────┴──────────┘

Task: pd_lkb_p_cluster
┌──────┬──────────────────────────────┬──────┬───────┬────────┬──────────┐
│ Rank │ Model                        │ MRR  │ Hit@5 │ Hit@10 │ Select   │
├──────┼──────────────────────────────┼──────┼───────┼────────┼──────────┤
│ 1    │ p001_run_20260528_001_best   │ .092 │ .100  │ .200   │ Use      │
│ 2    │ p001_run_demo_best           │ .016 │ 0.000 │ 0.000  │ Use      │
└──────┴──────────────────────────────┴──────┴───────┴────────┴──────────┘
```

---

# 12. Recommendation 页面

## 12.1 页面目的

Recommendation 页面用于使用某个选定模型，对实际候选列表进行推荐。

它可以从 Model Selector 的结果跳转过来，也可以手动选择 task 和 model。

## 12.2 页面布局

```text
Recommendation

Task source:
( ) Use existing task
( ) Create temporary task

Task:  [ ni_lkb_p_ni ▼ ]
Model: [ p001_run_demo_last ▼ ]

Optional override:
Candidates CSV: [ optional upload ]
Positives CSV:  [ optional upload ]
Negatives CSV:  [ optional upload ]

[Run Recommendation]
```

## 12.3 输出结果

```text
Recommendation Result

┌──────┬──────────────┬──────────────┬────────────┐
│ Rank │ Candidate ID │ SMILES       │ Score      │
├──────┼──────────────┼──────────────┼────────────┤
│ 1    │ L3_xxx       │ ...          │ 0.982      │
│ 2    │ L3_yyy       │ ...          │ 0.954      │
└──────┴──────────────┴──────────────┴────────────┘

[Download CSV]
[Save Result]
```

---

# 13. Jobs & Logs 页面

## 13.1 页面目的

Jobs & Logs 是给高级用户或维护者使用的调试页面。

普通用户只需要知道任务是否完成；维护者需要能查看底层 step 和日志。

## 13.2 页面内容

```text
Jobs & Logs

Running / Recent Jobs
┌──────────────┬──────────┬─────────┬────────────┬──────────┐
│ Job          │ Run      │ Status  │ Started At │ Log      │
├──────────────┼──────────┼─────────┼────────────┼──────────┤
│ Build Model  │ run_001  │ running │ ...        │ View     │
│ Evaluate     │ batch_01 │ success │ ...        │ View     │
└──────────────┴──────────┴─────────┴────────────┴──────────┘
```

Run detail 中可以显示：

```text
step1
step2
apply-ga
step6_7
step8
step12
step13
embedding
training
```

但这些不应作为主导航流程。

---

# 14. Settings 页面

## 14.1 页面目的

Settings 用于配置路径和刷新数据库。

## 14.2 页面内容

```text
Settings

Workspace root:
[ /path/to/ChemDBWebVersion/workspace ]

SQLite DB:
[ /path/to/ChemDBWebVersion/workspace/chemdb.sqlite ]

ChemDB repo root:
[ /path/to/ChemDBWebVersion/ChemDB ]

[Validate Paths]
[Rebuild SQLite Index]
```

`Rebuild SQLite Index` 对应：

```bash
python -m app.storage.cli rebuild-index \
  --workspace-root workspace \
  --db-path workspace/chemdb.sqlite
```

---

# 15. 暂缓的 Hyperparameter Library

超参数库暂时不进入 Phase 3A/3B 的核心范围。

原因：

```text
1. 当前部分训练参数可能仍然写死在 src 或 training 中。
2. 若过早抽象 hyperparameter_sets，会牵涉 training config、orchestrator、model registry 的进一步改造。
3. 目前最重要的是先让 GA、Task、Model Builder、Model Selector 这条主线跑通。
```

未来可以新增：

```text
workspace/hyperparameter_sets/
```

以及 SQLite 表：

```text
hyperparameter_sets
hyperparameter_versions
```

但这应放在前端主流程稳定之后再做。

当前 Model Builder 可以先使用默认训练配置，或只提供非常有限的参数选项。

---

# 16. Phase 3 实施顺序建议

## Phase 3A：只读资源库和状态展示

目标：

```text
先让用户看见系统中有什么。
```

实现页面：

```text
Dashboard
GA Library
Task Library
Model Library
Model Selector Results
Settings
```

特点：

```text
只读为主
不启动 pipeline
不编辑 GA
不创建 task
不训练模型
```

---

## Phase 3B：资源注入与模型构建

目标：

```text
让用户可以通过页面注入 GA、注入 task，并启动模型构建。
```

实现功能：

```text
Upload GA CSV
Create GA version
Upload Task
Create Task
Select Dataset + GA
Start Model Build
View build progress
Register models into Model Library
```

---

## Phase 3C：模型选择与推荐闭环

目标：

```text
让用户可以选择任务和模型范围，自动得到排序，并进一步做推荐。
```

实现功能：

```text
Select tasks
Select model scope
Run model selection
View ranked models
Choose best model
Run recommendation
Download recommendation result
Save recommendation result
```

---

# 17. 关键设计原则

## 17.1 前端主线不是 Step

不要让用户按下面方式操作：

```text
step1
step2
step4_5
step6_7
step8
step12
step13
training
```

这些是后台细节。

## 17.2 前端主线应该是 Resource → Model → Selection

用户应看到：

```text
GA Library
Task Library
Model Builder
Model Library
Model Selector
Recommendation
```

## 17.3 Run 是后台记录，不是用户主对象

Run 仍然重要，因为它记录一次模型构建过程。

但普通用户真正关心的是：

```text
我用了哪个 GA？
我生成了哪些模型？
这些模型在哪些 task 上表现如何？
哪个模型应该用于推荐？
```

## 17.4 SQLite 是前端主要读取来源

Streamlit 应主要读取：

```text
workspace/chemdb.sqlite
```

通过：

```text
app/storage/repositories.py
```

获取结构化数据。

## 17.5 Streamlit 不直接调用 ChemDB 科学代码

Streamlit 不应直接 import 或运行：

```text
ChemDB/src/step*.py
```

而应调用：

```text
app/services/orchestrator.py
app/services/ga_registry_manager.py
app/services/model_after_manager.py
```

---

# 18. 最终一句话总结

ChemDB Web 的前端不应该是命令行脚本的网页复制版，而应该是一个资源驱动的模型构建与模型选择系统：

```text
GA 列表 + 数据集
        ↓
自动构建模型
        ↓
模型库

Task 列表 + 模型范围
        ↓
自动模型选择
        ↓
推荐结果
```

Run、Step、Artifact、Log 作为后台记录存在，用于追踪、调试和复现，但不应成为普通用户的主要操作入口。
