# model_after — 改哪里

**数据流 / 路径 / I/O（一般不动）**  
`io.py`（task 读取）、`bundle.py`（run 加载）、`batch.py` 里 CSV 循环与写 summary 的骨架。

**科学逻辑（改这里）**

| 目标 | 文件 |
|------|------|
| 怎么打分、context 怎么建、候选怎么嵌入、rank 含义 | `scoring.py` |
| 单模型好不好：held-out MRR/Hit@k、margin、status | `eval_task.py` |
| 多模型里选谁：汇总指标 → `selection_score` 规则 | `batch.py` → `_selection_score()` |

单模型评估 ≈ `scoring.py` + `eval_task.py`；群中最适 ≈ 上两者产出 + `batch.py` 选模规则。
