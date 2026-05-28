# -*- coding: utf-8 -*-
"""
Data 模块：按 target 7:1:2 划分（仅从 tmp/m_l3_pairs.csv 取 T 索引），
单遍扫大 CSV 写 train/val/test records；Dataset 将 (T,M,KL,NL,UL) 转为 PyTorch 张量。
对应 .cursor/nn_structure.md 与 training 计划。
"""
from __future__ import annotations

import csv
import json
import logging
import os
import pickle
import random
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import torch
from tqdm import tqdm

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# DID 归一化：内部统一为 "D..." 形式，与 L3_embedding npz 查表一致
# ---------------------------------------------------------------------------


def normalize_did(did: str) -> str:
    """归一化 DID：去掉 L3_ 前缀，内部统一用 D... 形式。"""
    if not did or not did.strip():
        return ""
    s = did.strip()
    if s.upper().startswith("L3_"):
        return s[3:]
    return s


def to_l3_did(did: str) -> str:
    """转为 L3 格式（用于与 CSV 中 T 列一致）。"""
    if not did or not did.strip():
        return ""
    s = did.strip()
    if s.upper().startswith("L3_"):
        return s
    return f"L3_{s}"


# ---------------------------------------------------------------------------
# 从 tmp/m_l3_pairs.csv 获取唯一 T 并做 7:1:2 划分
# ---------------------------------------------------------------------------


def load_unique_t_from_m_l3_pairs(m_l3_pairs_path: str | Path) -> List[str]:
    """
    从 m_l3_pairs.csv 读取唯一 ligand_did，转为 L3 格式，返回有序列表。
    """
    path = Path(m_l3_pairs_path)
    if not path.is_file():
        raise FileNotFoundError(f"m_l3_pairs.csv not found: {path}")
    seen: set[str] = set()
    unique_t: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        lid_col = "ligand_did"
        if reader.fieldnames and lid_col not in reader.fieldnames:
            for c in reader.fieldnames or []:
                if "ligand" in c.lower() and "did" in c.lower():
                    lid_col = c
                    break
        for row in reader:
            did = (row.get(lid_col) or "").strip()
            if not did:
                continue
            t = to_l3_did(did)
            if t not in seen:
                seen.add(t)
                unique_t.append(t)
    logger.info("load_unique_t_from_m_l3_pairs: %d unique targets from %s", len(unique_t), path)
    return unique_t


def build_split_targets(unique_t: List[str], seed: int) -> Tuple[List[str], List[str], List[str]]:
    """7:1:2 划分：shuffle 后 70% train, 10% val, 20% test。"""
    rng = random.Random(seed)
    order = list(unique_t)
    rng.shuffle(order)
    n = len(order)
    n_train = int(round(0.7 * n))
    n_val = int(round(0.1 * n))
    n_test = n - n_train - n_val
    train_t = order[:n_train]
    val_t = order[n_train : n_train + n_val]
    test_t = order[n_train + n_val :]
    return train_t, val_t, test_t


# ---------------------------------------------------------------------------
# 单遍流式读大 CSV，按 T 归属写入 train/val/test_records.pkl
# ---------------------------------------------------------------------------


def _flush_t_block(
    prev_t: str,
    block: List[Tuple[str, List[str], List[str]]],
    ul_dids: List[str],
    train_set: set[str],
    val_set: set[str],
    test_set: set[str],
    train_records: List[Dict],
    val_records: List[Dict],
    test_records: List[Dict],
) -> None:
    """将上一 T 的所有 (M, kl_dids, nl_dids) 与 ul_dids 写入对应 split。"""
    for m, kl, nl in block:
        rec = {"T": prev_t, "M": m, "kl_dids": list(kl), "nl_dids": list(nl), "ul_dids": list(ul_dids)}
        if prev_t in train_set:
            train_records.append(rec)
        elif prev_t in val_set:
            val_records.append(rec)
        elif prev_t in test_set:
            test_records.append(rec)


def build_records_single_pass(
    csv_path: str | Path,
    train_targets: List[str],
    val_targets: List[str],
    test_targets: List[str],
    training_dir: str | Path,
    seed: int,
) -> Dict[str, Any]:
    """
    单遍流式读 step13_kl_nl_samples.csv，按 T 归属写入 train/val/test_records.pkl。
    返回统计信息（用于写 split_index.json）。
    """
    train_set = set(train_targets)
    val_set = set(val_targets)
    test_set = set(test_targets)

    train_records: List[Dict] = []
    val_records: List[Dict] = []
    test_records: List[Dict] = []

    path = Path(csv_path)
    if not path.is_file():
        raise FileNotFoundError(f"KL/NL/UL CSV not found: {path}")

    logger.info(
        "build_records_single_pass: reading %s (train_targets=%d, val=%d, test=%d)",
        path, len(train_targets), len(val_targets), len(test_targets),
    )

    # 当前 T 的 (M, kl_dids, nl_dids) 块缓冲；当前 T 的 UL 列表
    block: List[Tuple[str, List[str], List[str]]] = []
    current_ul: List[str] = []
    prev_t: Optional[str] = None
    current_m: Optional[str] = None
    current_kl: List[str] = []
    current_nl: List[str] = []

    def flush_current_t_m():
        nonlocal block, current_ul, prev_t, current_m, current_kl, current_nl
        if prev_t is None:
            return
        if current_m is not None and (current_kl or current_nl):
            block.append((current_m, current_kl, current_nl))
        if block:
            _flush_t_block(
                prev_t, block, current_ul,
                train_set, val_set, test_set,
                train_records, val_records, test_records,
            )
        block = []
        current_ul = []
        current_m = None
        current_kl = []
        current_nl = []

    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        row_iter = tqdm(
            reader,
            desc="Scanning KL/NL/UL CSV",
            unit=" rows",
            mininterval=1.0,
            leave=True,
        )
        for row in row_iter:
            t = (row.get("T") or "").strip()
            m = (row.get("M") or "").strip()
            label = (row.get("label") or "").strip()
            candidate_did = (row.get("candidate_did") or "").strip()
            did_norm = normalize_did(candidate_did) if candidate_did else ""

            if label == "UL":
                if did_norm:
                    current_ul.append(did_norm)
                continue

            # KL 或 NL
            if t != prev_t:
                flush_current_t_m()
                prev_t = t

            if m != current_m:
                if current_m is not None and (current_kl or current_nl):
                    block.append((current_m, current_kl, current_nl))
                current_m = m
                current_kl = []
                current_nl = []

            if label == "KL" and did_norm:
                current_kl.append(did_norm)
            elif label == "NL" and did_norm:
                current_nl.append(did_norm)

    flush_current_t_m()

    out_dir = Path(training_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Writing train_records.pkl (%d records)...", len(train_records))
    with open(out_dir / "train_records.pkl", "wb") as f:
        pickle.dump(train_records, f, protocol=4)
    logger.info("Writing val_records.pkl (%d records)...", len(val_records))
    with open(out_dir / "val_records.pkl", "wb") as f:
        pickle.dump(val_records, f, protocol=4)
    logger.info("Writing test_records.pkl (%d records)...", len(test_records))
    with open(out_dir / "test_records.pkl", "wb") as f:
        pickle.dump(test_records, f, protocol=4)
    logger.info(
        "build_records_single_pass done: train=%d val=%d test=%d",
        len(train_records), len(val_records), len(test_records),
    )
    return {
        "n_train": len(train_records),
        "n_val": len(val_records),
        "n_test": len(test_records),
    }


# ---------------------------------------------------------------------------
# 索引生成入口：仅当 split_index.json 不存在时执行
# ---------------------------------------------------------------------------

SPLIT_INDEX_NAME = "split_index.json"
DEFAULT_SEED = 42


def ensure_split_index(
    training_dir: str | Path,
    kl_nl_ul_csv_path: str | Path,
    m_l3_pairs_path: str | Path,
    seed: int = DEFAULT_SEED,
) -> Path:
    """
    若 training_dir/split_index.json 不存在，则从 m_l3_pairs 取唯一 T、7:1:2 划分，
    单遍扫大 CSV 写 train/val/test_records.pkl 并写 split_index.json。
    返回 split_index.json 的 Path。
    """
    training_dir = Path(training_dir)
    split_index_path = training_dir / SPLIT_INDEX_NAME
    if split_index_path.is_file():
        logger.info("split index already exists: %s", split_index_path)
        return split_index_path

    logger.info("Generating split index and records (seed=%s)...", seed)
    logger.info("Step 1/3: loading unique targets from m_l3_pairs: %s", m_l3_pairs_path)
    unique_t = load_unique_t_from_m_l3_pairs(m_l3_pairs_path)
    logger.info("  -> %d unique targets", len(unique_t))
    logger.info("Step 2/3: 7:1:2 split...")
    train_t, val_t, test_t = build_split_targets(unique_t, seed)
    logger.info("  -> train=%d val=%d test=%d targets", len(train_t), len(val_t), len(test_t))
    logger.info("Step 3/3: single-pass scan of KL/NL/UL CSV and write records...")
    stats = build_records_single_pass(
        kl_nl_ul_csv_path,
        train_t, val_t, test_t,
        training_dir,
        seed,
    )
    index_data = {
        "seed": seed,
        "target_list_source": "m_l3_pairs.csv",
        "train_targets": train_t,
        "val_targets": val_t,
        "test_targets": test_t,
        "n_train": stats["n_train"],
        "n_val": stats["n_val"],
        "n_test": stats["n_test"],
    }
    with open(split_index_path, "w", encoding="utf-8") as f:
        json.dump(index_data, f, indent=2, ensure_ascii=False)
    logger.info("Wrote %s (n_train=%s n_val=%s n_test=%s)", split_index_path, stats["n_train"], stats["n_val"], stats["n_test"])
    return split_index_path


# ---------------------------------------------------------------------------
# Embedding 查表
# ---------------------------------------------------------------------------


def load_ligand_embedding(npz_path: str | Path) -> Tuple[Dict[str, np.ndarray], int]:
    """
    加载 npz（dids, embeddings），构建归一化 did -> vec 查表。
    返回 (dict, d_l)。
    """
    path = Path(npz_path)
    with np.load(path, allow_pickle=True) as z:
        dids = z["dids"]
        emb = z["embeddings"]
    if dids.ndim == 0:
        dids = dids.reshape(-1)
    lookup: Dict[str, np.ndarray] = {}
    for i, did in enumerate(dids):
        d = did if isinstance(did, str) else str(did)
        key = normalize_did(d)
        if key and i < emb.shape[0]:
            lookup[key] = np.asarray(emb[i], dtype=np.float32)
    d_l = emb.shape[1] if emb.size else 0
    return lookup, d_l


def load_metal_embedding(csv_path: str | Path) -> Tuple[Dict[str, np.ndarray], int]:
    """
    按 element 列加载，构建 symbol -> vec 查表；排除 element 列本身，其余为数值特征。
    返回 (dict, d_m)。
    """
    path = Path(csv_path)
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    if not rows:
        return {}, 0
    elem_col = "element"
    if elem_col not in fieldnames and "symbol" in fieldnames:
        elem_col = "symbol"
    feature_cols = [c for c in fieldnames if c != elem_col and c.strip()]
    lookup: Dict[str, np.ndarray] = {}
    for row in rows:
        sym = (row.get(elem_col) or "").strip()
        if not sym:
            continue
        vec = []
        for c in feature_cols:
            try:
                vec.append(float(row.get(c, 0)))
            except (TypeError, ValueError):
                vec.append(0.0)
        lookup[sym] = np.array(vec, dtype=np.float32)
    d_m = len(feature_cols)
    return lookup, d_m


# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------


class RankDataset(torch.utils.data.Dataset):
    """
    条件排序 Dataset：每条样本 (T,M)，返回 metal (d_m), context (d_l), candidates (K, d_l), target_idx=0。
    候选集合仅 [T] + NL + UL；context 为 KL 的 mean pooling。
    """

    def __init__(
        self,
        split: str,
        training_dir: str | Path,
        ligand_embedding_path: str | Path,
        metal_embedding_path: str | Path,
        max_candidates: Optional[int] = None,
        kl_pooling: str = "mean",
        dtype: torch.dtype = torch.float32,
        device: torch.device = torch.device("cpu"),
        missing_embedding: str = "zero",
    ):
        assert split in ("train", "val", "test"), split
        self.split = split
        self.training_dir = Path(training_dir)
        self.dtype = dtype
        self.device = device
        self.missing_embedding = missing_embedding  # "zero" | "skip"
        self.max_candidates = max_candidates
        self.kl_pooling = kl_pooling

        records_path = self.training_dir / f"{split}_records.pkl"
        if not records_path.is_file():
            raise FileNotFoundError(f"Records not found: {records_path}. Run ensure_split_index first.")
        with open(records_path, "rb") as f:
            raw_records = pickle.load(f)

        # 只保留候选数 = 1(T)+128(NL)+32(UL) = 161 的条目，其余丢弃
        expected_candidates = 161
        self.records = [
            r for r in raw_records
            if (1 + len(r.get("nl_dids") or []) + len(r.get("ul_dids") or [])) == expected_candidates
        ]
        dropped = len(raw_records) - len(self.records)
        if dropped:
            logger.info(
                "RankDataset(%s): dropped %d records with candidates != %d (kept %d)",
                split, dropped, expected_candidates, len(self.records),
            )

        self.ligand_lookup, self.d_l = load_ligand_embedding(ligand_embedding_path)
        self.metal_lookup, self.d_m = load_metal_embedding(metal_embedding_path)
        self._zero_ligand = np.zeros(self.d_l, dtype=np.float32)
        self._zero_metal = np.zeros(self.d_m, dtype=np.float32)

    def __len__(self) -> int:
        return len(self.records)

    def _get_ligand_vec(self, did: str) -> np.ndarray:
        key = normalize_did(did) if did else ""
        if key and key in self.ligand_lookup:
            return self.ligand_lookup[key]
        if self.missing_embedding == "zero":
            return self._zero_ligand.copy()
        return None

    def _get_metal_vec(self, symbol: str) -> np.ndarray:
        if symbol and symbol in self.metal_lookup:
            return self.metal_lookup[symbol]
        if self.missing_embedding == "zero":
            return self._zero_metal.copy()
        return None

    def __getitem__(self, i: int) -> Dict[str, torch.Tensor]:
        rec = self.records[i]
        t = rec["T"]
        m = rec["M"]
        kl_dids = rec.get("kl_dids") or []
        nl_dids = rec.get("nl_dids") or []
        ul_dids = rec.get("ul_dids") or []

        metal_vec = self._get_metal_vec(m)
        if metal_vec is None:
            metal_vec = self._zero_metal.copy()

        # context: KL mean pooling
        kl_vecs = []
        for d in kl_dids:
            v = self._get_ligand_vec(d)
            if v is not None:
                kl_vecs.append(v)
        if kl_vecs:
            context = np.stack(kl_vecs, axis=0).mean(axis=0).astype(np.float32)
        else:
            context = self._zero_ligand.copy()

        # candidates = [T] + NL + UL
        cand_vecs: List[np.ndarray] = []
        t_vec = self._get_ligand_vec(t)
        if t_vec is not None:
            cand_vecs.append(t_vec)
        for d in nl_dids:
            v = self._get_ligand_vec(d)
            if v is not None:
                cand_vecs.append(v)
        for d in ul_dids:
            v = self._get_ligand_vec(d)
            if v is not None:
                cand_vecs.append(v)

        K = self.max_candidates
        if K is not None and len(cand_vecs) > K:
            cand_vecs = cand_vecs[:K]
        if K is not None and len(cand_vecs) < K:
            for _ in range(K - len(cand_vecs)):
                cand_vecs.append(self._zero_ligand.copy())
        candidates = np.stack(cand_vecs, axis=0).astype(np.float32)
        if K is not None:
            assert candidates.shape[0] == K, (candidates.shape[0], K)

        def to_t(x: np.ndarray) -> torch.Tensor:
            return torch.as_tensor(x, dtype=self.dtype, device=self.device)

        return {
            "metal": to_t(metal_vec),
            "context": to_t(context),
            "candidates": to_t(candidates),
            "target_idx": torch.tensor(0, dtype=torch.long, device=self.device),
        }


def load_index_paths(index_path: Optional[str | Path] = None) -> Dict[str, str]:
    """
    从 training/index.json 读取路径。若未传 index_path 则使用本文件同目录下的 index.json。
    返回 dict: kl_nl_ul_csv, m_l3_pairs（若 index 中无则用 kl_nl_ul 同目录的 m_l3_pairs.csv 推断）,
    ligand_embedding_gin/gcn/ecfp, metal_embedding。
    """
    if index_path is None:
        index_path = Path(__file__).resolve().parent / "index.json"
    with open(index_path, "r", encoding="utf-8") as f:
        idx = json.load(f)
    csv_path = idx["kl_nl_ul_index"]["path"]
    base = str(Path(csv_path).parent)
    m_l3_default = os.path.join(base, "m_l3_pairs.csv")
    out = {
        "kl_nl_ul_csv": csv_path,
        "m_l3_pairs": idx.get("m_l3_pairs_path") or m_l3_default,
        "metal_embedding": idx["metal_embedding"]["path"],
    }
    for k, v in idx["ligand_embedding"]["variants"].items():
        out[f"ligand_embedding_{k}"] = v["path"]
    return out


def get_dataloader(
    split: str,
    training_dir: str | Path,
    ligand_embedding_path: str | Path,
    metal_embedding_path: str | Path,
    batch_size: int = 32,
    max_candidates: Optional[int] = None,
    shuffle: Optional[bool] = None,
    num_workers: int = 0,
    **kwargs: Any,
) -> torch.utils.data.DataLoader:
    """
    构建 DataLoader。train 默认 shuffle=True（无放回随机），val/test 默认 shuffle=False。
    """
    if shuffle is None:
        shuffle = split == "train"
    dataset = RankDataset(
        split=split,
        training_dir=training_dir,
        ligand_embedding_path=ligand_embedding_path,
        metal_embedding_path=metal_embedding_path,
        max_candidates=max_candidates,
    )
    return torch.utils.data.DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
        **kwargs,
    )


# ---------------------------------------------------------------------------
# 全量运行入口：仅当无 split_index.json 时生成划分与 records
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    parser = argparse.ArgumentParser(
        description="全量运行 data：若不存在 split_index.json，则从 m_l3_pairs 做 7:1:2 划分并单遍扫大 CSV 写 train/val/test_records.pkl。"
    )
    parser.add_argument(
        "--training-dir",
        type=str,
        default=None,
        help="输出目录（split_index.json 与 *_records.pkl），默认为本文件所在目录",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="随机种子")
    parser.add_argument(
        "--index",
        type=str,
        default=None,
        help="index.json 路径，默认 training/index.json",
    )
    args = parser.parse_args()

    base = Path(__file__).resolve().parent
    training_dir = Path(args.training_dir) if args.training_dir else base
    index_path = Path(args.index) if args.index else (base / "index.json")

    paths = load_index_paths(index_path)
    ensure_split_index(
        training_dir=training_dir,
        kl_nl_ul_csv_path=paths["kl_nl_ul_csv"],
        m_l3_pairs_path=paths["m_l3_pairs"],
        seed=args.seed,
    )
    print("Done. Split index and records in:", training_dir)
