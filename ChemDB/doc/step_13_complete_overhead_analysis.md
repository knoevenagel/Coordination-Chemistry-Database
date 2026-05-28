# step_13_complete.py 相对 step13_with_kl_nl.py 的增量步骤与开销分析

## 一、增量功能概述

`step_13_complete.py` 在 `step13_with_kl_nl.py` 的基础上**仅增加**：为每个 target 生成 **num_unrelated（默认 32）个“无关 NL”**，label 为 `UL`。  
“无关”定义：L5(l) ∩ L5(target) = ∅。

因此增量步骤全部围绕 **无关 L3 的抽样与写入**，以及为支撑该逻辑所增加的**数据与调用**。

---

## 二、增量步骤拆解与开销

### 1. 索引加载阶段（主进程，一次）

| 步骤 | 说明 | 与 with_kl_nl 的差异 |
|------|------|----------------------|
| 加载 `l3_l5.json` / `l5_l3.json` 等 | 两脚本均加载相同文件 | 无 |
| **预计算全体 L3 列表** | `all_l3_list = list(get_all_l3_set(_SHARED_L5_L3_INDEX))` | **仅 complete 有** |

- **`get_all_l3_set(l5_l3_index)`**  
  - 遍历 `l5_l3_index` 的所有 value（每个 L5 下的 L3 列表），合并到 set，再转 list。  
  - 复杂度：O(Σ |L3 列表|) ≈ O(总 L3 条数)，约 10^6 量级。  
  - 时间：一次遍历 + set 去重 + list 化，通常为索引加载时间的附加尾段（与 JSON 解析、磁盘 I/O 相比通常较小）。  
- 实测中 complete 的 `index_load_time`（10.28s）比 with_kl_nl（3.27s）大，可能还受运行顺序、磁盘缓存等影响；`get_all_l3_set` 本身是单次 O(N) 的合理开销。

---

### 2. 每个 target 的 KL/NL 之后：无关 NL（UL）生成

对**每个**成功构建候选集且参与 KL/NL 的 target，在写完 KL/NL 行之后会调用：

```text
sample_unrelated_l3_for_target(target_did, k=32, l3_l5_map, l5_l3_index, all_l3_list, base_seed, target_num)
```

该函数内部步骤与复杂度如下。

#### 2.1 构建 `pool`（每 target 一次）

```python
pool = [x for x in all_l3_list if x != normalized and x != target_did]
```

- **含义**：从约 10^6 的 `all_l3_list` 中筛掉当前 target 的两种写法，得到 `pool`。  
- **复杂度**：O(|all_l3_list|) ≈ O(10^6) **每个 target**。  
- **开销**：100 个 target 即约 **10^8** 次迭代与列表构造，是**主要 CPU 开销之一**（大量内存遍历与 list 分配）。

#### 2.2 按批抽样并判定“无关”

- **循环**：`loop_num = 1, 2, ...`，每轮用 `seed_tu(base_seed, target_num, loop_num)` 固定种子，从 `pool` 中 `rng.sample(pool, 128)` 抽 128 个 L3。  
- **对每个被抽到的 L3**：  
  - 调用 `_is_unrelated_l3(l3_l5_map, l, target_did)`：  
    - 取 L5(l) 与 L5(target)，做 `set` 再求交，判断是否为空。  
  - 若无关且未在 `result` 中则加入，直到凑满 32 个。

- **`_is_unrelated_l3` 单次**：  
  - 2 次 `l3_l5_map.get` + 2 个 set 构造 + 1 次交集，O(|L5(l)| + |L5(target)|)，通常为小常数。  
- **每 target 调用次数**：最少约 32 次（第一批 128 里就凑够），最多可能多批（每批 128 次），平均约几十到一两百次。  
- 100 个 target × 若干批 × 128 ≈ **数万次** `_is_unrelated_l3`，单次成本小，但累加明显。

#### 2.3 写入 UL 行

- 对 32 个 UL 各写一行 CSV，与 KL/NL 写入类似，开销相对抽样可忽略。

---

## 三、开销量化与占比（100 targets 实测）

基于你之前的对比（同一环境、100 样本、kl-nl-only、complete 使用已有 metal_l3_index）：

| 指标 | step13_with_kl_nl | step_13_complete | 增量 |
|------|--------------------|------------------|------|
| **build_time** | 9.57 s | 97.19 s | **+87.62 s** |
| **总 processing_time** | 53.50 s | 120.44 s | **+66.94 s**（其余为 extract/index_load 等） |

- **增量几乎全部来自“构建阶段”**，即每个 target 的候选集构建 + KL/NL + **UL 抽样**。  
- 候选集构建与 KL/NL 部分两脚本一致，因此 **~87 s 中的绝大部分可归因于 UL 相关逻辑**。

粗略分解 UL 相关开销：

1. **每 target 一次 `pool = [x for x in all_l3_list if ...]`**  
   - 100 × 10^6 级迭代 → 易达 **数十秒** 量级（视 CPU/内存带宽而定）。  
2. **`sample_unrelated_l3_for_target` 内**  
   - 多轮 `rng.sample(pool, 128)` + 大量 `_is_unrelated_l3` → 估计 **数秒到十几秒**。  
3. **`get_all_l3_list` 一次**  
   - 与 1 相比可忽略。

因此可以概括：**“每 target 对百万级 all_l3_list 做一次 list 推导”是当前实现中最突出的增量开销。**

---

## 四、优化建议（针对 UL 抽样）

1. **去掉“每 target 重建 pool”**  
   - 不在 `sample_unrelated_l3_for_target` 里做 `pool = [x for x in all_l3_list if x != normalized and x != target_did]`。  
   - 直接使用 `all_l3_list`，每批 `rng.sample(all_l3_list, 128)`；在遍历 batch 时跳过 `l == normalized or l == target_did`，再调用 `_is_unrelated_l3`。  
   - 这样每 target 从 O(|all_l3_list|) 降为 O(批数 × 128)，可显著减少 CPU 与内存分配。

2. **可选：预计算“全体 L3”的 list 一次，并保证可复现**  
   - 若对顺序无要求，可对 `get_all_l3_set` 的结果做一次 `sorted` 再 `list`，保证 `all_l3_list` 固定，随机性完全由 `seed_tu` + `rng.sample` 控制。

3. **批量或并行**  
   - 若主进程仍是瓶颈，可考虑把“UL 抽样”挪到 worker 中与候选集/KL-NL 同批处理，或对多个 target 的 UL 抽样做轻量并行（需注意 CSV 写入串行化）。

实施建议 1 即可明显降低 step_13_complete 的增量开销，且不改变语义（仍为同一批 L3、同一套“无关”定义与种子逻辑）。

---

## 五、小结

| 增量步骤 | 发生位置 | 复杂度/次数 | 主要开销 |
|----------|----------|-------------|----------|
| `get_all_l3_set` → `all_l3_list` | 主进程，索引加载后一次 | O(总 L3 数) | 小 |
| 每 target 构建 `pool` | `sample_unrelated_l3_for_target` 内 | 100 × O(10^6) | **大** |
| 每 target 多批 `rng.sample` + `_is_unrelated_l3` | 同上 | 100 × 多批 × 128 | 中 |
| 写入 32 行 UL/ target | 主进程 CSV | 100 × 32 | 小 |

**结论**：step_13_complete 相对 step13_with_kl_nl 的额外时间，主要来自**每个 target 对全体 L3 列表做一次 O(N) 的 list 推导以得到 pool**。优化该处（改为直接基于 `all_l3_list` 抽样并在循环内跳过 target 自身）即可显著降低增量开销。

---

## 六、优化实施与复测结果（2026-02-02）

- **优化**：`sample_unrelated_l3_for_target` 不再每 target 构建 `pool`，改为直接对 `all_l3_list` 做 `rng.sample(all_l3_list, 128)`，在循环内跳过 `l == normalized or l == target_did`。
- **target 一致性**：两脚本均改为使用 `random.Random(self.seed).sample(unique_dids, self.num_samples)` 抽样 target，同种子（42）则同 100 个 target；复测时两脚本的 target 集合经 diff 验证一致，总候选数均为 643,894。

| 指标 | step13_with_kl_nl | step_13_complete（优化后） | 优化前 complete（参考） |
|------|-------------------|----------------------------|-------------------------|
| **build_time** | 10.66 s | **37.80 s** | 97.19 s |
| **总 processing_time** | 52.83 s | 60.12 s | 120.44 s |
| **墙钟** | 56.10 s | 70.23 s | 129.72 s |
| **平均每 target 构建** | 0.107 s/个 | 0.378 s/个 | 0.972 s/个 |

优化后 complete 的 build 由约 97 s 降至约 38 s（约 **2.6× 加速**），增量开销由约 87 s 降至约 27 s。

---

## 七、复测消耗对比（2026-02-02，KL/NL 按 (T,M) 种子一致后）

条件：100 样本、seed=42、同一 metal_l3_index（bench_kl_nl 已有索引），先跑 with_kl_nl 再跑 complete。

| 指标 | step13_with_kl_nl | step_13_complete |
|------|-------------------|------------------|
| **墙钟时间** | 56.71 s | **20.82 s** |
| **脚本 processing_time** | 53.21 s | 17.56 s |
| **extract_time** | 0.79 s | 0.58 s |
| **index_load_time** | 3.31 s | 3.77 s |
| **build_time** | 10.76 s | **9.39 s** |
| **平均每 target 构建** | 0.108 s/个 | 0.094 s/个 |
| **CPU 占用** | 183% | 302% |
| **User+System 时间** | 66.20+37.62 s | 34.91+28.02 s |
| **最大常驻内存 (RSS)** | 2190 MB | 2195 MB |

说明：本次 complete 墙钟与 build 均略优于 with_kl_nl，可能因 complete 后跑、metal_l3_index 等文件已被缓存；且 complete 并行度更高（CPU 302% vs 183%）。同一环境下多轮跑或对调顺序，数值会有波动，但优化后 complete 的 UL 增量已明显小于早期版本。
