# ChemDB 流水线输入输出与文件依赖

本文档描述 ChemDB 项目从原始数据到推荐训练/评估的**完整文件依赖关系**。路径默认相对于项目根目录 `ChemDB_restructured/ChemDB/`。

**相关文档**：技术报告 [`report.md`](report.md)；Step13→评估细节 [`../training/results_aromatic_small_rings/training_tmp_results/PIPELINE_STEP13_TO_EVAL.md`](../training/results_aromatic_small_rings/training_tmp_results/PIPELINE_STEP13_TO_EVAL.md)。

---

## 1. 约定

| 符号 | 含义 |
|------|------|
| `[USER]` | 需用户自备的输入 |
| `[TMP]` | 写入 `tmp/`，流水线中间/最终表 |
| `[DATA]` | 写入 `data/`，长期嵌入与配置 |
| `[TRAIN]` | 写入 `training/` |
| `(opt)` | 可选步骤或可选文件 |
| `(del)` | 运行结束后会删除的临时文件 |
| `-->` | 生成 / 依赖方向 |

**工作目录**：建议在项目根执行 `python src/main.py` 或 `python -m training.*`，使 `./tmp`、`./data` 落在项目根下。若在 `src/` 内单独跑 `step*.py`，`./tmp` 会变为 `src/tmp`，与训练路径不一致。

**编排入口 `src/main.py` 覆盖范围**：`step1` → `step2` → `step4_5` → `step6_7` → `step8` → `step9` → `step11`。**不包含** `step12`、`step13` 与 `training/*`。

**ChemDBWebVersion `app.services.orchestrator`（Phase 1E+）**：`step1` → `step2` → `require_bound_ga` → `apply_ga_to_run` → `step6_7` → …  
GA 列表为 **workspace** 资源（`workspace/ga_sets/...`），绑定后物化到 `run/tmp/GA_with_id.csv`；core **不**自动跑 `step4_5` full-auto。Authoring：`step4_5 --mode generate-ga`；下游：`--mode apply-ga`。

---

## 2. 总览依赖树（ASCII）

以下为可逐行解析的依赖树：`+--` 表示子节点，`-->` 表示产出。

```
ChemDB_ROOT/
|
+-- [USER] data/pubchem/*.csv ............... cid, isosmiles
+-- [USER] data/metal_list.txt .............. 金属元素列表
+-- [USER] data/p_elements_list.txt ......... (opt) P区元素
+-- [USER] data/metal_embedding/element_features.csv
|
+-- STEP1  src/step1.py  (src/main.py)
|   IN:  data/pubchem/*.csv, data/metal_list.txt
|   OUT: [TMP] complex_data.csv
|        [TMP] ligand_data.csv
|        [TMP] processing_stats.json
|        (opt) [TMP] pubchem_data.csv
|
+-- STEP2  src/step2.py
|   IN:  [TMP] ligand_data.csv
|   OUT: [TMP] repaired_ligand_data.csv  .... ligand_new_did, new_smiles, ...
|        [TMP] repair_stats.json
|
+-- STEP4_5  src/step4_5.py
|   IN:  [TMP] repaired_ligand_data.csv
|   OUT: [TMP] ligand_data_deduplicated.csv
|        [TMP] GA_with_id.csv ............. GA_SMILES, GA_ID
|        [TMP] ligand_with_gac.csv ........ DID, SMILES, GAC
|        [TMP] IRL_filtered.csv
|        [TMP] IRL_filtered_cleaned.csv ... (step6_7 不用此文件)
|        [TMP] step4_5_stats.json
|   TMP: (del) ga_patterns_temp.csv
|
+-- STEP6_7  src/step6_7.py
|   IN:  [TMP] repaired_ligand_data.csv
|        [TMP] GA_with_id.csv
|        [TMP] IRL_filtered.csv
|   OUT: [TMP] fragments.csv
|        [TMP] step6_7_stats.json
|
+-- STEP8  src/step8.py
|   IN:  [TMP] complex_data.csv, ligand_data.csv, repaired_ligand_data.csv
|        [TMP] fragments.csv, IRL_filtered.csv
|        data/metal_list.txt
|   OUT: [TMP] neo4j_metals.csv
|        [TMP] neo4j_complexes.csv, neo4j_ligands.csv
|        [TMP] neo4j_repaired_ligands.csv, neo4j_fragments.csv, neo4j_irl.csv
|        [TMP] neo4j_m_l1_relationships.csv
|        [TMP] neo4j_l1_l2_relationships.csv, neo4j_l2_l3_relationships.csv
|        [TMP] neo4j_l1_l3_relationships.csv
|        [TMP] neo4j_l3_l4_relationships.csv, neo4j_l4_l5_relationships.csv
|        [TMP] step8_stats.json
|   (opt) tools/neo4j/import.sh --> Neo4j (bolt 3087, browser 3074)
|
+-- STEP9 (opt)  src/step9.py
|   IN:  [TMP] complex_data.csv, ligand_data.csv, repaired_ligand_data.csv
|        [TMP] fragments.csv
|        data/IRL_filtered.csv ............ 注意: 在 data/ 而非 tmp/
|   OUT: [TMP] did_index.csv
|        [TMP] step9_stats.json
|
+-- STEP11 (opt)  src/step11.py  [需 GPU / MolCLR]
|   IN:  [TMP] repaired_ligand_data.csv
|        [TMP] neo4j_m_l1_relationships.csv, neo4j_l1_l3_relationships.csv
|        [TMP] neo4j_metals.csv
|   OUT: [TMP] l3_embeddings.csv, m_l3_relationships.csv
|        [TMP] step11_stats.json
|
+-- STEP12  src/step12.py
|   IN:  [TMP] neo4j_l3_l4_relationships.csv
|        [TMP] neo4j_l4_l5_relationships.csv
|        [TMP] ligand_with_gac.csv (或回退 neo4j_repaired_ligands.csv)
|        [TMP] neo4j_m_l1_relationships.csv, neo4j_l1_l3_relationships.csv
|        [TMP] neo4j_repaired_ligands.csv
|   OUT: [TMP] l3_l5.json
|        [TMP] l5_l3.json
|        [TMP] l5_freq_weight.json
|        [TMP] l3_gac.json
|        [TMP] l5_l3_filtered_K{K}.json .... 默认 K=30
|        [TMP] m_l3_pairs.csv
|        [TMP] step12_stats.json
|
+-- STEP13  src/step_13_complete.py
|   IN:  [TMP] m_l3_pairs.csv, l3_l5.json, l5_l3.json
|        [TMP] l5_freq_weight.json, l3_gac.json
|        (opt) [TMP] l5_l3_filtered_K{K}.json
|        (opt) [TMP] metal_l3_index.csv .... 缺则 Neo4j 自动生成
|   OUT: [TMP] step13_kl_nl_samples.csv ... 训练主表 (T,M,label,candidate_did,...)
|        [TMP] step13_kl_nl_samples.stats.csv
|        [TMP] step13_stats.json
|        (opt) step13_targets.csv, step13_candidates.csv, metal_l3_index.csv
|
+-- EMBED_L3  src/tools/build_L3_embedding_index.py
|   IN:  [TMP] metal_l3_index.csv, repaired_ligand_data.csv
|   OUT: [DATA] L3_embedding/L3_embedding_gin.npz
|        [DATA] L3_embedding/L3_embedding_gcn.npz
|        [DATA] L3_embedding/L3_embedding_ecfp.npz
|
+-- EMBED_METAL  data/metal_embedding/zscore_element_features.py
|   IN:  [DATA] metal_embedding/element_features.csv
|   OUT: [DATA] metal_embedding/element_features_zscore.csv
|
+-- TRAIN_DATA  training/data.py
|   IN:  [TRAIN] index.json --> step13_kl_nl_samples.csv, m_l3_pairs.csv
|   OUT: [TRAIN] split_index.json
|        [TRAIN] train_records.pkl, val_records.pkl, test_records.pkl
|   NOTE: 仅当 split_index.json 不存在时重建 pkl
|
+-- TRAIN  training/train.py
|   IN:  [TRAIN] config.yaml, split_index.json, *_records.pkl
|        [TRAIN] index.json --> ligand *.npz, metal zscore csv
|   OUT: [TRAIN] ckpts/best.pt, last.pt, history.json
|
+-- EVAL  training/evaluation/run_all_quests.py
|   IN:  [TRAIN] models/*.pt, evaluation/quests/*.csv, index.json
|   OUT: [TRAIN] evaluation/results/*
```

---

## 3. 分层模型与数据流（ASCII）

化学分层（L 层）及主要载体文件：

```
                    PubChem CSV
                         |
                         v
              +---------------------+
              |  Step1: M/L1/L2   |
              +---------------------+
                    |         |
            complex_data    ligand_data
                    |         |
                    |         v
                    |    Step2 repair
                    |         |
                    |    repaired_ligand_data  <-- L3 配体主表
                    |         |
                    |    Step4/5 GA + IRL + GAC
                    |    GA_with_id, IRL_filtered, ligand_with_gac
                    |         |
                    |    Step6/7 L4 fragments
                    |         |
                    v         v
              +---------------------+
              |  Step8: Neo4j CSV  |  neo4j_* + 关系表
              +---------------------+
                         |
              +----------+----------+
              |                     |
         Step12 索引            (opt) Neo4j DB
    l3_l5, l5_l3, l3_gac,
    m_l3_pairs, ...
              |
              v
         Step13 采样
    step13_kl_nl_samples.csv
              |
    +---------+---------+
    |                   |
 L3 npz            metal zscore
    |                   |
    v                   v
 training.data --> *_records.pkl
    |
    v
 training.train --> ckpts/*.pt
    |
    v
 run_all_quests / 推理脚本
```

---

## 4. 脚本 IO 速查表

### 4.1 数据处理（`src/`）

| 脚本 | 命令示例 | 必需输入 | 主要输出 |
|------|----------|----------|----------|
| `step1.py` | `python src/step1.py [workers] [limit]` | `data/pubchem/*.csv`, `data/metal_list.txt` | `tmp/complex_data.csv`, `tmp/ligand_data.csv` |
| `step2.py` | `python src/step2.py` | `tmp/ligand_data.csv` | `tmp/repaired_ligand_data.csv` |
| `step4_5.py` | `python src/step4_5.py [--workers N] [--limit N]` | `tmp/repaired_ligand_data.csv` | `tmp/GA_with_id.csv`, `tmp/ligand_with_gac.csv`, `tmp/IRL_filtered.csv`, … |
| `step6_7.py` | `python src/step6_7.py` | `tmp/repaired_ligand_data.csv`, `tmp/GA_with_id.csv`, `tmp/IRL_filtered.csv` | `tmp/fragments.csv` |
| `step8.py` | `python src/step8.py` | Step1–6/7 的 `tmp/*.csv`（见上） | `tmp/neo4j_*.csv`（节点+关系） |
| `step9.py` | `python src/step9.py` | 多个 `tmp/*.csv` + `data/IRL_filtered.csv` | `tmp/did_index.csv` |
| `step11.py` | `python src/step11.py` | `tmp/repaired_ligand_data.csv`, 部分 `neo4j_*` 关系 | `tmp/l3_embeddings.csv`, `tmp/m_l3_relationships.csv` |
| `step12.py` | `cd src && python step12.py --input-dir ../tmp --output-dir ../tmp --K 30` | Step8 关系 CSV + `ligand_with_gac.csv` | `tmp/l3_l5.json`, `m_l3_pairs.csv`, … |
| `step_13_complete.py` | `cd src && python step_13_complete.py --input-dir ../tmp --output-dir ../tmp --kl-nl-only` | Step12 全部 JSON/CSV | `tmp/step13_kl_nl_samples.csv` |
| `main.py` | `python src/main.py` | 同 step1–11 所选步骤 | 同各步输出 |

### 4.2 嵌入与训练

| 脚本 | 命令示例 | 必需输入 | 主要输出 |
|------|----------|----------|----------|
| `build_L3_embedding_index.py` | `python src/tools/build_L3_embedding_index.py --tmp-dir tmp` | `tmp/metal_l3_index.csv`, `tmp/repaired_ligand_data.csv` | `data/L3_embedding/L3_embedding_{gin,gcn,ecfp}.npz` |
| `zscore_element_features.py` | `python data/metal_embedding/zscore_element_features.py` | `data/metal_embedding/element_features.csv` | `element_features_zscore.csv` |
| `training/data.py` | `python -m training.data --seed 42` | `training/index.json` 指向的 CSV + `m_l3_pairs.csv` | `training/split_index.json`, `*_records.pkl` |
| `training/train.py` | `python -m training.train --config training/config.yaml` | `config.yaml`, pkl, `index.json` 嵌入路径 | `training/ckpts/best.pt`, `last.pt`, `history.json` |
| `run_all_quests.py` | `python -m training.evaluation.run_all_quests` | `training/models/*.pt`, `quests/*.csv`, `index.json` | `training/evaluation/results/` |

### 4.3 服务与工具（只读 `tmp/` 或查询）

| 组件 | 端口 | 读取的数据 |
|------|------|------------|
| `src/server.py` | 3041 | `tmp/*.csv`, `*_stats.json` |
| `src/molclr_api.py` | 3044 | 运行时算嵌入 |
| `src/affinity_api.py` | 3045 | 模型 + SMILES |
| `src/proxy.py` | 3040 | 反代上述服务 |
| `tools/chemtool.py` | — | SMILES → DID |
| `tools/neo4j/import.sh` | — | `tmp/neo4j_*.csv` |

---

## 5. 关键文件字段说明

### 5.1 `tmp/repaired_ligand_data.csv`（枢纽表）

| 列名 | 说明 |
|------|------|
| `ligand_new_did` | 修复后 L3 DID |
| `source_did` | 修复前 DID |
| `new_smiles` | 修复后 SMILES |
| `old_smiles` | 原始 SMILES |

下游 Step4/5、6/7、embedding 均依赖 **`ligand_new_did` + `new_smiles`**。

### 5.2 `tmp/step13_kl_nl_samples.csv`（训练索引源）

| 列名 | 说明 |
|------|------|
| `T` | 目标配体 DID（L3 格式） |
| `M` | 金属符号；UL 行可为空 |
| `label` | `KL` / `NL` / `UL` |
| `candidate_did` | 候选配体 DID |
| `distance`, `d_tilde`, `weight` | 采样与权重元数据 |

### 5.3 `training/index.json`

集中配置路径（当前仓库内为**绝对路径**，迁移需修改）：

- `kl_nl_ul_index.path` → `step13_kl_nl_samples.csv`
- `ligand_embedding.variants.{gin,gcn,ecfp}.path` → `data/L3_embedding/*.npz`
- `metal_embedding.path` → `element_features_zscore.csv`

### 5.4 `training/*_records.pkl`

每条 record 字典键：`T`, `M`, `kl_dids`, `nl_dids`, `ul_dids`（由 `training/data.py` 从 Step13 CSV 流式生成）。

---

## 6. 按目录的文件清单

### 6.1 `data/` — 用户输入与嵌入产物

```
data/
├── pubchem/
│   └── *.csv                    [USER] Step1 输入 (cid, isosmiles)
├── metal_list.txt               [USER] Step1, Step8
├── p_elements_list.txt          [USER] (opt)
├── IRL_filtered.csv             [USER] 仅 Step9 读取；与 tmp 副本可能不同步
├── L3_embedding/
│   ├── L3_embedding_gin.npz     <-- build_L3_embedding_index.py
│   ├── L3_embedding_gcn.npz
│   └── L3_embedding_ecfp.npz
└── metal_embedding/
    ├── element_features.csv     [USER]
    └── element_features_zscore.csv  <-- zscore_element_features.py
```

### 6.2 `tmp/` — 流水线主工作区

```
tmp/
├── complex_data.csv             <-- step1
├── ligand_data.csv              <-- step1
├── repaired_ligand_data.csv     <-- step2  *** 枢纽 ***
├── ligand_data_deduplicated.csv <-- step4_5
├── GA_with_id.csv               <-- step4_5
├── ligand_with_gac.csv          <-- step4_5
├── IRL_filtered.csv             <-- step4_5
├── IRL_filtered_cleaned.csv     <-- step4_5 (step6_7 不用)
├── fragments.csv                <-- step6_7
├── neo4j_*.csv                  <-- step8 (节点 + 关系)
├── l3_l5.json                   <-- step12
├── l5_l3.json
├── l5_freq_weight.json
├── l3_gac.json
├── l5_l3_filtered_K30.json      <-- step12 (K 可配置)
├── m_l3_pairs.csv               <-- step12
├── metal_l3_index.csv           <-- step13 或 embedding 脚本
├── step13_kl_nl_samples.csv     <-- step13  *** 训练主表 ***
├── did_index.csv                <-- step9 (opt)
├── l3_embeddings.csv            <-- step11 (opt)
├── processing_stats.json        各 step 的 *_stats.json
└── (del) ga_patterns_temp.csv   step4_5 内部临时
```

### 6.3 `training/` — 训练与评估

```
training/
├── index.json                   路径配置
├── config.yaml                  训练超参
├── split_index.json             <-- data.py (7:1:2 划分)
├── train_records.pkl            <-- data.py (体积可达 GB 级)
├── val_records.pkl
├── test_records.pkl
├── ckpts/
│   ├── best.pt                  <-- train.py
│   ├── last.pt
│   └── history.json
├── models/                      评估用 checkpoint 副本
└── evaluation/
    ├── quests/*.csv             评估任务 (如 ni_lkb_p_ni.csv)
    └── results/                 run_all_quests 输出
```

---

## 7. 推荐执行顺序（ASCII 流程）

### 7.1 建库路径（至 Neo4j CSV）

```
[USER] pubchem + metal_list
    --> step1 --> step2 --> step4_5 --> step6_7 --> step8
                                              |
                                    (opt) neo4j import
                                    (opt) step9
```

### 7.2 推荐训练路径（完整）

```
step8 产出 neo4j_*.csv
    --> step12 --> step13 --> step13_kl_nl_samples.csv
          |                        |
          |                        +--> build_L3_embedding_index --> *.npz
          +--> m_l3_pairs.csv              |
                                           v
                              zscore_element_features --> metal zscore
                                           |
                                           v
                              training/data.py --> *_records.pkl
                                           |
                                           v
                              training/train.py --> ckpts/*.pt
                                           |
                                           v
                              run_all_quests --> evaluation/results/
```

### 7.3 `main.py` 与完整链的关系

```
main.py:  step1 -> step2 -> step4_5 -> step6_7 -> step8 -> step9 -> step11
完整推荐: 上述建库 + step12 + step13 + embedding + training/*
```

---

## 8. 临时文件与缓存行为

| 文件/机制 | 类型 | 行为 |
|-----------|------|------|
| `ga_patterns_temp.csv` | 真临时 | Step4/5 生成 `GA_with_id.csv` 后 **删除** |
| `split_index.json` + `*_records.pkl` | 缓存 | **已存在则 `training/data.py` 不重建**；更新 Step13 后需手动删除 |
| `training/index.json` | 配置 | 使用绝对路径时需随环境更新 |
| `pubchem_data.csv` | 可选大文件 | Step1 `SAVE_PUBCHEM_DATA` 控制 |
| GAC 评估 `tmp_gac_*` | 实验目录 | `run_gac_eval.sh` 可能另建临时 `tmp` 副本 |

---

## 9. 步骤间硬依赖（机器可读简表）

以下为 `REQUIRES` 关系，格式：`下游脚本 REQUIRES 输入文件`（文件缺失将导致失败或跳过）。

```
step1          REQUIRES  data/pubchem/*.csv, data/metal_list.txt
step2          REQUIRES  tmp/ligand_data.csv
step4_5        REQUIRES  tmp/repaired_ligand_data.csv
step6_7        REQUIRES  tmp/repaired_ligand_data.csv, tmp/GA_with_id.csv, tmp/IRL_filtered.csv
step8          REQUIRES  tmp/complex_data.csv, tmp/ligand_data.csv, tmp/repaired_ligand_data.csv,
                       tmp/fragments.csv, tmp/IRL_filtered.csv, data/metal_list.txt
step9          REQUIRES  tmp/complex_data.csv, tmp/ligand_data.csv, tmp/repaired_ligand_data.csv,
                       tmp/fragments.csv, data/IRL_filtered.csv
step11         REQUIRES  tmp/repaired_ligand_data.csv, tmp/neo4j_m_l1_relationships.csv,
                       tmp/neo4j_l1_l3_relationships.csv, tmp/neo4j_metals.csv
step12         REQUIRES  tmp/neo4j_l3_l4_relationships.csv, tmp/neo4j_l4_l5_relationships.csv,
                       tmp/ligand_with_gac.csv (或 neo4j_repaired_ligands.csv),
                       tmp/neo4j_m_l1_relationships.csv, tmp/neo4j_l1_l3_relationships.csv,
                       tmp/neo4j_repaired_ligands.csv
step13         REQUIRES  tmp/m_l3_pairs.csv, tmp/l3_l5.json, tmp/l5_l3.json,
                       tmp/l5_freq_weight.json, tmp/l3_gac.json
build_L3_emb   REQUIRES  tmp/metal_l3_index.csv, tmp/repaired_ligand_data.csv
training.data  REQUIRES  step13_kl_nl_samples.csv, m_l3_pairs.csv (经 index.json)
training.train REQUIRES  split_index.json, *_records.pkl, embeddings (经 index.json)
```

---

## 10. 常见问题

1. **工作目录不一致**：在 `src/` 下跑 step 会使 `tmp/` 位置与 `training/index.json` 不一致。
2. **Step9 的 IRL 路径**：读 `data/IRL_filtered.csv`，Step4/5 默认写在 `tmp/IRL_filtered.csv`；若跑 Step9 需复制或改路径。
3. **`main.py` 不含 Step12/13**：完整推荐链需单独执行。
4. **embedding 与 Step13 顺序**：`build_L3_embedding_index.py` 需要 `metal_l3_index.csv`，通常先 Step13（或 Neo4j 生成索引）再建 npz。
5. **pkl 体积**：`train_records.pkl` 可达数 GB，注意磁盘与内存。

---

## 11. 文档修订

| 版本 | 说明 |
|------|------|
| 1.0 | 初始版：覆盖 Step1–13、embedding、training、eval 全链路 IO 与 ASCII 依赖树 |

若脚本默认路径或 `main.py` 编排变更，请同步更新本节与第 2、9 节。
