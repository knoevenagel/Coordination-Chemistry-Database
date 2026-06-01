# ChemDB `src/` 与 `training/` 参数全量普查（规则/阈值/分块/采样）

## A. `ChemDB/src` 全量清单（逐脚本）

## 1) `affinity_api.py`

- `DEFAULT_K = 10`（kNN 邻居数）
  - 用于控制亲和度估计时参与投票的近邻数量。
- `DEFAULT_MODEL_TYPE = "gin"`
  - 用于指定默认加载的分子表示模型类型。
- `API_PORT = 3045`
  - 用于定义服务默认监听端口。
- `calculate_affinity(..., k=DEFAULT_K)` 默认 k=10
  - 用于在请求级别覆盖近邻数量。
- GPU 选择逻辑硬编码：选 `nvidia-smi` 报告的最大 `memory.free` GPU
  - 用于优先选择空闲显存最多的设备以减少 OOM。

## 2) `comparison.py`

- `compare_pair(..., sample_mismatches=0)`
  - 用于控制每对比较时额外抽样输出的错配数量。
- CLI:
  - `--sample-mismatches` 默认 `0`
    - 用于设置错配样本抽样数，`0` 表示不抽样。
  - `--pairs` 默认空字符串
    - 用于指定要比较的 pair 列表来源。

## 3) `fingerprint_utils.py`

- `_get_morgan_generator(radius=2, fp_size=2048)`
  - `radius=2` 用于定义 Morgan 指纹的邻域半径。
  - `fp_size=2048` 用于定义指纹位向量长度。

## 4) `main.py`

- 不是算法参数中心，主要是 step 调度映射 `STEP_DEFINITIONS/STEP_ORDER`
- 内部兼容参数：
  - `record_limit = -1`（默认不限制）
    - 用于在 legacy 流程中控制最大处理记录数。

## 5) `molclr_api.py`

- `get_model(..., device='cpu')`
  - 用于控制模型加载时使用的默认设备。
- `extract_embeddings(..., device='cpu')`
  - 用于控制嵌入提取时使用的默认设备。

## 6) `proxy.py`

- 后端地址硬编码：
  - `SIM:3042`, `QRY:3041`, `SPA:3043`, `EMB:3044`
    - 用于固定四类后端服务的转发目标端口。

## 7) `server.py`

- SQLite 批写硬编码：`batch_size = 1000`
  - 用于控制每次写入 SQLite 的提交批量大小。
- PRAGMA 固定策略：
  - `journal_mode=OFF`
    - 用于关闭日志以提升写入吞吐。
  - `synchronous=OFF`
    - 用于降低刷盘同步强度以换取速度。
  - `temp_store=MEMORY`
    - 用于将临时数据放入内存减少磁盘 IO。
  - `cache_size=-200000`
    - 用于增大 SQLite 页缓存以减少频繁读写。

## 8) `step1.py`（清洗前处理）

- 默认常量：
  - `DEFAULT_RECORD_LIMIT = -1`
    - 用于默认不限制处理记录数。
  - `DEFAULT_MAX_WORKERS = 200`
    - 用于限制自动并发的最大 worker 数量。
  - `SAVE_PUBCHEM_DATA = True`
    - 用于控制是否保留 PubChem 中间产物文件。
  - `SORT_COMPLEX_BY_SMILES_LENGTH = False`（也可由 `CHEMDB_STEP1_SORT_BY_LENGTH` 覆盖）
    - 用于控制是否按 SMILES 长度排序复合物。
- worker 规则：
  - 自动并发：`min(cpu_count, 200)`
    - 用于在机器核数与上限之间自动取较小值。
  - worker `get(timeout=1)`（轮询超时）
    - 用于控制队列轮询阻塞时长。
  - stall 日志周期 `120s`
    - 用于控制卡顿状态日志的输出频率。
- 大 complex 可疑阈值（由 `step1` 调用 `CC_split.questionable_detection`）：
  - `ligand_num_boundary = 11`
    - 用于将配体片段数过多的 complex 标记为 `inactive=3`。
- CLI:
  - `--workers`
    - 用于手动覆盖并发 worker 数。
  - `--limit`（默认 `-1`）
    - 用于手动限制处理记录上限。

## 9) `step2.py`（配体修复）

- 默认常量：
  - `DEFAULT_RECORD_LIMIT = -1`
    - 用于默认全量处理输入数据。
- worker 与队列硬编码：
  - `worker_process(..., batch_size=100)`
    - 用于控制 worker 聚合回传结果的粒度。
  - `self.batch_size = max(50, 1000 // self.num_workers)`
    - 用于按并发度自适应调整每批任务规模。
  - 文件流式分块：`chunk_size = 10000`
    - 用于控制 CSV 读取的内存占用与吞吐平衡。
  - 队列回压阈值：`while task_queue.qsize() >= 4000`
    - 用于防止生产速度过快导致队列堆积。
  - 队列容量：`Queue(maxsize=self.num_workers * 5000)`
    - 用于限制最大排队任务数避免内存暴涨。
  - worker 自动上限：`min(cpu_count, 200)`
    - 用于限制自动并发上界。
- CLI:
  - `--workers`
    - 用于手工设置并发进程数。
  - `--limit`
    - 用于手工限制处理总量。

## 10) `step4_5.py`（GA/IRL）

- 默认常量：
  - `DEFAULT_RECORD_LIMIT = -1`
    - 用于默认不截断输入数据。
- 规则阈值（函数内硬编码）：
  - IRL 预筛：`GAC >= 2`
    - 用于剔除 GAC 过低的候选复合物。
  - 额外排除：`complex_smiles == 'CC'`
    - 用于移除已知无效或噪声样本。
- 并发/分块：
  - worker 上限：`min(num_workers, 200)`
    - 用于限制并发上界以避免资源竞争。
  - `chunksize = max(1, len(args_list)//(workers*4))`
    - 用于平衡进程调度开销与任务吞吐。
  - tqdm `mininterval=5.0`
    - 用于限制进度条刷新频率减少开销。
- CLI:
  - `--mode {full-auto, generate-ga, apply-ga}`
    - 用于切换不同的 GA/IRL 执行模式。
  - `--workers`
    - 用于指定并发 worker 数。
  - `--limit`
    - 用于限制处理记录数量。

## 11) `step6_7.py`（标注+拆分）

- 默认常量：
  - `DEFAULT_RECORD_LIMIT = -1`
    - 用于默认全量处理样本。
- 并发与分块硬编码：
  - worker cap：`min(num_workers, 64)`
    - 用于限制拆分流程并发上界。
  - 默认 `batch_size = 50`
    - 用于控制单批任务处理规模。
  - producer 读 CSV `batch_size = 10000`
    - 用于控制生产者侧流式读取块大小。
  - 回压阈值：`while task_queue.qsize() >= 4000`
    - 用于防止队列积压过大。
  - 默认 `max_queue_size = workers * 5000`
    - 用于限制队列容量随并发线性扩展。
- CLI:
  - `--workers`
    - 用于手动配置并发进程数。
  - `--limit`
    - 用于手动设置记录上限。
  - `--batch-size`
    - 用于手动设置单批处理规模。

## 12) `step8.py`

- 默认常量：
  - `DEFAULT_RECORD_LIMIT = -1`
    - 用于默认不进行截断。
- 数据截断规则：
  - 对多个实体表统一采用 `head(record_limit)`（仅在 `limit > 0` 时）
    - 用于在调试或小样本模式下保持多表规模一致。
- CLI:
  - `--limit`
    - 用于设置统一截断条数。
  - 路径参数 `--input-dir/--output-dir/--metal-list`
    - 用于指定输入输出与金属名单文件位置。

## 13) `step9.py`

- 默认常量：
  - `DEFAULT_RECORD_LIMIT = -1`
    - 用于默认全量处理输入。
- 并发与分块硬编码：
  - worker cap：`min(cpu_count, 200)`
    - 用于限制自动并发上限。
  - `worker_process(..., batch_size=100)`
    - 用于控制每次回传的数据批量。
  - `self.batch_size = max(10, 500 // workers)`
    - 用于按并发动态调整任务切片大小。
  - producer `chunk_size = 10000`
    - 用于控制生产者读入分块粒度。
  - 回压阈值：`while qsize >= workers * 500`
    - 用于限制队列积压防止内存压力。
  - 队列容量：`workers * 1000`
    - 用于定义最大排队任务容量。
- CLI（legacy argv）：
  - workers / record_limit
    - 用于兼容旧入口参数传递方式。

## 14) `step10.py`

- `Step10Processor(record_limit=-1)`
  - 用于控制全局处理条数上限，`-1` 表示不限制。
- 分块读取硬编码：`chunksize = 10000`
  - 用于控制流式读取时的单块大小。
- 全局集合达到 `record_limit` 提前停止
  - 用于在达到目标规模后提前结束流程。

## 15) `step11.py`

- 默认常量：
  - `DEFAULT_BATCH_SIZE = 256`
    - 用于控制嵌入推理时的批大小。
  - `DEFAULT_RECORD_LIMIT = -1`
    - 用于默认不限制处理总量。
- GPU 策略：
  - 自动选最大空闲显存 GPU
    - 用于自动选择更不易 OOM 的设备。
- 处理细节：
  - 批次数计算 `(N + batch_size - 1) // batch_size`
    - 用于计算完整处理所需的 batch 总数。
- 入口支持参数（legacy argv）：
  - `batch_size`, `record_limit`, `model_type`, `device`
    - 用于通过命令行覆盖默认推理参数。

## 16) `step12.py`

- 默认常量：
  - `DEFAULT_RECORD_LIMIT = -1`
    - 用于默认全量构建索引。
  - `DEFAULT_K = 30`
    - 用于控制候选筛选时的邻近截断阈值。
  - `DEFAULT_MIN_GAC = 10`
    - 用于限制最小 GAC 过滤下界。
  - `DEFAULT_MAX_GAC = 15`
    - 用于限制最大 GAC 过滤上界。
- CLI:
  - `--K`
    - 用于覆盖默认 `K` 值。
  - `--min-gac`
    - 用于覆盖最小 GAC 阈值。
  - `--max-gac`
    - 用于覆盖最大 GAC 阈值。
  - `--record-limit`
    - 用于限制本次处理记录数。

## 17) `step_13_complete.py`（采样核心）

- 默认常量：
  - `DEFAULT_NUM_SAMPLES = 100`
    - 用于设置默认采样目标数量。
  - `DEFAULT_D = 1000`
    - 用于定义 rare/common 的分界计数阈值。
  - `DEFAULT_K = 30`
    - 用于控制候选筛选邻域范围。
  - `DEFAULT_SEED = 42`
    - 用于确保采样流程可复现。
  - `DEFAULT_RECORD_LIMIT = -1`
    - 用于默认不截断输入记录。
- 采样规则硬编码（关键）：
  - `seed_tu(base_seed, target_num, l5_num)`（稳定种子混合，含固定大常数）
    - 用于为不同 target 组合生成稳定且分散的随机种子。
  - UL 采样循环批大小：`batch_size = 128`
    - 用于控制 UL 候选批次遍历粒度。
  - KL 权重：`alpha=4, lambda_param=1`
    - 用于强化 KL 采样中高相关候选的权重。
  - NL 权重：`alpha=2, lambda_param=0.7`
    - 用于平衡 NL 采样中的探索与利用。
  - KL 抽样数：`k=16`
    - 用于限制每个目标输出的 KL 数量。
  - NL 抽样数：`k=128`
    - 用于限制每个目标输出的 NL 数量。
  - 默认 UL 数量：`num_unrelated = 32`
    - 用于限制每个目标的无关样本数量。
  - 多进程 `imap_unordered(..., chunksize=10)`
    - 用于控制多进程任务切分粒度。
  - 批提交常量 `BATCH_SIZE = 5000`
    - 用于控制结果批量写出与合并节奏。
- CLI:
  - `--num-samples`, `--D`, `--K`, `--seed`
    - 用于覆盖采样规模、分界阈值、候选范围和随机种子。
  - `--num-unrelated`
    - 用于覆盖 UL 数量。
  - `--workers`
    - 用于设置并发处理进程数。
  - `--generate-kl-nl / --no-generate-kl-nl`
    - 用于控制是否生成 KL/NL 数据集。
  - `--kl-nl-only / --no-kl-nl-only`
    - 用于控制是否仅执行 KL/NL 相关流程。

## 18) `utils.py`

- 未发现显式阈值/采样类参数（工具函数为主）。

---

## B. `ChemDB/src/tools` 全量清单（拆分/清洗细节）

## 1) `CC_split.py`

- 函数默认参数：
  - `questionable_detection(ligand_num_boundary=11)`
    - 用于定义疑似异常复合物的 ligand 数边界。
- 含义：
  - 超过边界则标记 inactive（状态 3）
    - 用于将超阈值样本标为不可用。

## 2) `CC_split_fast.py`

- 与 `CC_split.py` 对齐：
  - `questionable_detection(ligand_num_boundary=11)`
    - 用于保持快速版本与原版一致的异常判定阈值。

## 3) `fix_smiles.py`（修复规则最重）

- 构造参数：
  - `__init__(smiles, debug=False, draw=False)`
    - `debug` 用于开启调试日志输出而 `draw` 用于控制是否绘图辅助排错。
- 硬编码规则表：
  - `self.typical_valences`（元素典型价态全表，硬编码）
    - 用于在修复时约束每种元素允许的价态范围。
  - `self.stable_isotopes`（稳定同位素映射表，硬编码）
    - 用于校验或修复同位素标注的合法性。
  - `self.culmul_map / culmul_map2`（结构修复映射表）
    - 用于对特定可识别错误结构执行规则化替换。
  - 卤素集合：`{F, Cl, Br, I, At}`（原子序数与符号两套规则）
    - 用于触发卤素相关的专门修复分支。
- 电荷分离修复规则：
  - `charge_se` 计数阈值 `< 2` 才继续修复
    - 用于避免在电荷分离过多时继续进行高风险修复。
- 其他细节：
  - 多处 “最接近价态”/氢原子增减规则均在函数内硬编码判定。
    - 用于把不满足化学约束的结构拉回可接受状态。

## 4) `DID_calculate.py`

- 未发现典型“阈值/采样”参数；
- 为 DID 计算算法实现，主要是哈希/替换盒流程细节。

## 5) `DID_embedding.py`

- 默认参数：
  - `backend='molclr_gin'`
    - 用于指定默认嵌入后端。
  - `radius=2`
    - 用于指定 ECFP 的邻域半径。
  - `n_bits=2048`
    - 用于指定 ECFP 向量长度。
  - `format='json'`
    - 用于指定默认输出格式。
- CLI 对应同名参数。
  - 用于让命令行可覆盖默认嵌入配置。

## 6) `L4_create.py`

- 构造参数：
  - `__init__(..., ga_list=None, source_did=None)`
    - 用于限定可用 GA 列表并可指定单一来源 DID。
- 未发现数值阈值型参数（以结构处理规则为主）。

## 7) `L4_create_independent.py`

- 构造参数与 `L4_create.py` 类似：
  - `ga_list=None`, `source_did=None`
    - 用于限定独立构建流程的来源范围。
- 另外存在固定 DID 常量：`CC_DID = "D218966891838592"`。
  - 用于标识固定参照复合物 DID。

## 8) `build_L3_embedding_index.py`

- 默认参数：
  - `radius=2`
    - 用于定义分子指纹提取半径。
  - `n_bits=2048`
    - 用于定义指纹向量维度。
  - `chunk_size=50000`
    - 用于控制索引构建时的流式处理块大小。
  - `estimate const = 1_000_000`
    - 用于估算内存或索引容量时的基准常量。
  - `device='cpu'`
    - 用于指定默认执行设备。
  - `backends=['gin','gcn','ecfp']`
    - 用于指定默认启用的多种嵌入后端。
- 估算模式样本数规则：
  - `sample_n = min(200, N) else max(50, N)`
    - 用于在估算模式下控制最小与最大采样规模。

## 9) `distance.py`

- 未发现显式阈值/采样参数（数学工具函数为主）。

---

## C. `ChemDB/training` 全量清单（逐脚本）

## 1) `data.py`

- 默认常量：
  - `DEFAULT_SEED = 42`
    - 用于固定数据切分与采样的随机性以保证复现。
- 切分规则（硬编码）：
  - `train/val/test = 7:1:2`
    - 用于定义训练、验证、测试集的默认比例。
- DataLoader 默认：
  - `batch_size=32`
    - 用于控制训练时每步样本数。
  - `num_workers=0`
    - 用于控制 DataLoader 并发加载进程数。
  - `max_candidates=None`
    - 用于控制候选上限，`None` 表示不截断。
- 记录过滤硬编码：
  - `expected_candidates = 161`（不满足即丢弃）
    - 用于过滤候选数不符合预期的数据记录。
- CLI:
  - `--seed`（默认 42）
    - 用于从命令行覆盖随机种子。

## 2) `model.py`

- `RankModel(..., dropout=0.1)` 默认 dropout
  - 用于在训练中做随机失活以减轻过拟合。

## 3) `train.py`（训练超参主入口）

- `resolve_config()` 默认回填：
  - Data:
    - `embedding='gin'`
      - 用于指定默认分子嵌入特征来源。
    - `batch_size=32`
      - 用于控制每个训练 step 的样本批大小。
    - `num_workers=0`
      - 用于控制数据加载并发进程数。
    - `seed=42`
      - 用于固定训练流程随机性。
  - Model:
    - `activation='gelu'`
      - 用于指定前馈网络的默认激活函数。
    - `dropout=0.1`
      - 用于控制默认失活比例。
  - Train:
    - `epochs=50`
      - 用于控制最大训练轮数。
    - `lr=1e-3`
      - 用于控制优化器的初始学习率。
    - `weight_decay=0.0`
      - 用于控制 L2 正则强度。
    - `optimizer='adamw'`
      - 用于指定默认优化器类型。
    - `scheduler=None`
      - 用于指定是否启用学习率调度器。
    - `scheduler_patience=3`
      - 用于定义调度器等待性能停滞的轮数。
    - `scheduler_factor=0.5`
      - 用于定义调度器触发时学习率衰减倍率。
    - `early_stop_patience=8`
      - 用于定义早停可容忍的连续无提升轮数。
    - `val_every=1`
      - 用于定义验证评估间隔轮次。
    - `device='cuda'`
      - 用于设置默认训练设备。
- 调度规则：
  - `scheduler='reduce_on_plateau'` 或 `scheduler='cosine'`
    - 用于在两种学习率调度策略之间切换。
  - early-stop 使用 `epochs_no_improve` 比较 `early_stop_patience`
    - 用于在验证集长期无提升时提前终止训练。
- CLI:
  - `--config`
    - 用于指定外部 YAML 配置文件路径。
  - `--resume`
    - 用于从已有 checkpoint 恢复训练状态。

## 4) `__init__.py`

- 无参数。

---

## D. `ChemDB/training/model_after` 全量清单

## 1) `scoring.py`

- ECFP 默认：
  - `compute_ecfp_vector(radius=2, n_bits=2048)`
    - 用于定义推荐打分中的默认指纹参数。
- 推荐默认设备：
  - `recommend_single_model(..., device='cpu')`
    - 用于控制单模型推荐时的默认推理设备。

## 2) `bundle.py`

- 设备默认：
  - `resolve_device(device='cpu')`
    - 用于在未显式指定时回退到 CPU 运行。
- 模型构建默认回退：
  - `dropout = m_cfg.get('dropout', 0.1)`
    - 用于在配置缺失时保证模型仍有默认正则强度。

## 3) `eval_task.py`

- 命中率 Top-K 固定：
  - `hit@5`, `hit@10`, `hit@20`, `hit@50`
    - 用于统一报告多粒度排名命中指标。
- `max_k` 逻辑按 positives 数量与传参约束
  - 用于防止评价 K 值超过可计算范围。

## 4) `batch.py`

- 多模型排序规则：
  - `selection_score` 降序
    - 用于按综合评分由高到低排列候选模型。
  - 失败模型留空 rank
    - 用于避免失败任务污染有效排名序列。

## 5) `io.py`, `__init__.py`

- 无超参数性质阈值。

---

## E. 与 `theory.md` 的关键一致性核对

- 理论给出的：
  - `D=1000`
    - 用于定义 rare/common 分界阈值。
  - `K=30`
    - 用于定义候选筛选邻域规模。
  - KL `(alpha=4, lambda=1)`
    - 用于定义 KL 采样的权重函数形状。
  - NL `(alpha=2, lambda=0.7)`
    - 用于定义 NL 采样的权重函数形状。
  - `|NL|=128`
    - 用于定义每个目标的 NL 抽样数量。
- 当前代码：
  - 上述值基本一致
    - 表示关键超参整体与理论文档一致。
  - **差异**：`|KL|` 当前实现是 `k=16`（不是 64）
    - 表示当前 KL 抽样规模比理论设定更小。

---

## F. 备注

- 本清单聚焦“可调与硬编码参数点”，不展开所有化学规则实现细节代码。
- 若要继续，可把本文件再升级为“参数字典”格式（参数名 -> 文件 -> 行号 -> 当前值 -> 是否可配置 -> 建议暴露方式）。

---

## G. 统计目的/范围/方法（放在最后）

### 目的

- 对 `src` 与 `training` 中与行为相关的参数做系统盘点，覆盖可配项与硬编码项。
- 特别识别“拆分/清洗流程”中的规则、阈值、分块与采样细节。

### 范围

- 覆盖目录：
  - `ChemDB/src/**/*.py`
  - `ChemDB/src/tools/**/*.py`
  - `ChemDB/training/**/*.py`
  - `ChemDB/training/model_after/**/*.py`
- 覆盖类型：
  - CLI 参数（`argparse`）
  - `DEFAULT_*` 与常量赋值
  - 函数默认参数
  - 函数体内硬编码规则（阈值、分块大小、采样 k、队列阈值等）

### 方法

1. 通过静态扫描提取参数候选（`add_argument`、默认值、常量、关键字参数）。
2. 对关键脚本人工复核，补全“函数体硬编码”的业务语义（如 Step13 的 KL/NL 采样逻辑）。
3. 以“先 `src`、后 `training`”顺序组织，并在末尾做 `theory.md` 一致性核对。

