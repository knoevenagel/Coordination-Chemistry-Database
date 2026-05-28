#!/usr/bin/env python3
"""
ChemDB Step13 Complete (with KL/NL + 无关 NL)
基于 step13_with_kl_nl，在处理每个 target 的 NL 时额外生成与 target 无关的 NL（label=UL），
数量可配置，默认 32。无关 NL 抽样逻辑内联于此脚本，使用 seed_tu(base_seed, target_num, loop_num) 控制每轮 128 个随机 L3。

输入/输出同 step13_with_kl_nl。
"""

# Default configuration variables
DEFAULT_RECORD_LIMIT = -1  # -1 for unlimited
DEFAULT_NUM_SAMPLES = 100  # Default number of samples to process
OUTPUT_DIR = "./tmp"
INPUT_DIR = "./tmp"
DEFAULT_D = 1000  # Default rare/common L5 threshold
DEFAULT_K = 30    # Default GAC cutoff
DEFAULT_SEED = 42 # Default random seed

import os
import csv
import json
import random
import re
import time
import logging
import gc
import math
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict
import pandas as pd
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool, cpu_count

# Neo4j import (optional, only needed for index generation)
try:
    from neo4j import GraphDatabase
    NEO4J_AVAILABLE = True
except ImportError:
    NEO4J_AVAILABLE = False


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def get_optimal_worker_count() -> int:
    """Get optimal worker count based on CPU cores"""
    cpu_cores = cpu_count()
    optimal_workers = min(cpu_cores, 200)
    logger.info(f"CPU核心数: {cpu_cores}, 使用工作进程数: {optimal_workers}")
    return optimal_workers


def extract_did_number(did: str) -> int:
    """
    从DID字符串中提取15位数字部分并转换为整数
    
    Args:
        did: DID字符串，格式为 "D" + 15位数字，如 "D000001209752003" 或 "L3_D000001209752003"
    
    Returns:
        提取的数字部分作为整数，如 1209752003
    """
    # 匹配 D 后面跟着的数字
    match = re.search(r'D(\d+)', did)
    if match:
        return int(match.group(1))
    else:
        # 如果没有匹配到，尝试直接提取所有数字
        digits = ''.join(filter(str.isdigit, did))
        if digits:
            return int(digits)
        raise ValueError(f"无法从DID中提取数字: {did}")


def seed_tu(base_seed: int, target_num: int, l5_num: int) -> int:
    """
    为每个(target_did, L5_did)对生成稳定、可复现、分散性好的随机种子
    
    使用纯整数混合：两次大质数乘法 + 异或 + 轻量avalanche
    
    Args:
        base_seed: 基础随机种子（通常来自命令行参数）
        target_num: target_did的数字部分（15位数字）
        l5_num: l5_did的数字部分（15位数字）
    
    Returns:
        31位非负整数种子（范围：0 到 2^31-1）
    """
    # 大质数用于混合（选择两个较大的质数）
    P1 = 2654435761  # 2^32 * 0.618 (黄金比例相关)
    P2 = 2246822519  # 另一个大质数
    
    # 第一次混合：base_seed 和 target_num
    h1 = (base_seed * P1) ^ (target_num * P2)
    
    # 第二次混合：加入 l5_num
    h2 = (h1 * P1) ^ (l5_num * P2)
    
    # 轻量avalanche效果：位移和乘法扰动
    h3 = h2 ^ (h2 >> 16)
    h3 = h3 * P1
    h3 = h3 ^ (h3 >> 15)
    
    # Mask到31位（0x7FFFFFFF = 2^31 - 1）
    seed = h3 & 0x7FFFFFFF
    
    return seed


# Module-level globals for shared data (inherited via fork, NOT pickled)
# These are set by the main process BEFORE Pool creation
_SHARED_L3_L5_MAP = None
_SHARED_L5_L3_INDEX = None
_SHARED_L5_L3_INDEX_FILTERED = None
_SHARED_L5_STATS = None
_SHARED_L3_GAC_MAP = None
_SHARED_L5_WEIGHT_MAP = None

# Worker-specific parameters (small, OK to pickle)
_worker_D = None
_worker_K = None
_worker_seed = None
_worker_compute_distances = None


def _init_worker(D: int, K: int, seed: int, compute_distances: bool):
    """Initialize worker process - only sets small parameters.
    Large data is already inherited via fork (copy-on-write).
    """
    global _worker_D, _worker_K, _worker_seed, _worker_compute_distances
    
    # Suppress logging in worker processes
    import logging
    logging.getLogger().setLevel(logging.ERROR)
    
    # Only set small parameters - large data already inherited via fork
    _worker_D = D
    _worker_K = K
    _worker_seed = seed
    _worker_compute_distances = compute_distances


def _process_single_did(did: str) -> Tuple[bool, str, Optional[Dict], float]:
    """Worker function to process a single DID using inherited shared data"""
    global _worker_D, _worker_K, _worker_seed, _worker_compute_distances
    
    start_time = time.time()
    try:
        result = _build_cand_with_shared_data(
            target_did=did,
            l3_l5_map=_SHARED_L3_L5_MAP,
            l5_l3_index=_SHARED_L5_L3_INDEX,
            l5_l3_index_filtered=_SHARED_L5_L3_INDEX_FILTERED,
            l5_stats=_SHARED_L5_STATS,
            l3_gac_map=_SHARED_L3_GAC_MAP,
            l5_weight_map=_SHARED_L5_WEIGHT_MAP,
            D=_worker_D,
            K=_worker_K,
            seed=_worker_seed,
            compute_distances=_worker_compute_distances
        )
        elapsed_time = time.time() - start_time
        return (True, did, result, elapsed_time)
    except Exception as e:
        elapsed_time = time.time() - start_time
        return (False, did, str(e), elapsed_time)


def _build_cand_with_shared_data(
    target_did: str,
    l3_l5_map: dict,
    l5_l3_index: dict,
    l5_l3_index_filtered: Optional[dict],
    l5_stats: dict,
    l3_gac_map: dict,
    l5_weight_map: dict,
    D: int,
    K: int,
    seed: Optional[int],
    compute_distances: bool
) -> Dict:
    """
    Build candidate set using shared data (no CandidateBuilder instance needed).
    This is a standalone function that operates on shared data references.
    
    为每个(target_did, L5_did)对使用独立的随机种子，确保采样结果独立且可复现。
    """
    import random
    
    # Step 1: Get G_T = L5(T)
    if target_did not in l3_l5_map:
        raise ValueError(f"目标配体 {target_did} 不存在于L3-L5映射中")
    
    g_t_list = sorted(l3_l5_map[target_did])
    g_t = set(g_t_list)
    
    # 提取target_did的数字部分（用于生成种子）
    try:
        target_num = extract_did_number(target_did)
    except (ValueError, AttributeError) as e:
        # 如果提取失败，使用hash作为fallback（但这不是理想情况）
        target_num = hash(target_did) & 0x7FFFFFFF
        logger.warning(f"无法从target_did提取数字，使用hash: {target_did}, {e}")
    
    # Step 2: Partition rare/common L5 units
    g_t_rare = []
    g_t_common = []
    
    for u in g_t_list:
        if u not in l5_stats:
            continue
        freq = l5_stats[u]["freq"]
        if freq < D:
            g_t_rare.append(u)
        else:
            g_t_common.append(u)
    
    # 按 theory.md：B(u)=S(u)\{T}，cand(T) 中不包含 T 本身
    target_norm = normalize_did(target_did)
    exclude_t = {target_did, target_norm}
    
    # Step 3a: Rare L5 full expansion (B(u)=S(u)\setminus{T})
    candidates = set()
    rare_expansion_count = 0
    
    for u in g_t_rare:
        if u in l5_l3_index:
            l3_list = l5_l3_index[u]
            added = [l for l in l3_list if l not in exclude_t]
            candidates.update(added)
            rare_expansion_count += len(added)
    
    # Step 3b: Common L5 restricted expansion
    # 为每个common L5单元使用独立的随机种子
    common_expansion_count = 0
    use_filtered_index = (l5_l3_index_filtered is not None)
    index_to_use = l5_l3_index_filtered if use_filtered_index else l5_l3_index
    
    # base_seed用于seed_tu函数（如果seed为None，使用默认值0）
    base_seed = seed if seed is not None else 0
    
    for u in g_t_common:
        if u not in index_to_use:
            continue
        
        if use_filtered_index:
            b_u = [l for l in index_to_use[u] if l not in exclude_t]
        else:
            all_l3_with_u = l5_l3_index[u]
            b_u = []
            for l3_did in all_l3_with_u:
                if l3_did in exclude_t:
                    continue
                gac = l3_gac_map.get(l3_did)
                if gac is not None and gac <= K:
                    b_u.append(l3_did)
        b_u.sort()
        
        if len(b_u) > 0:
            # 为当前(target_did, L5_did)对生成独立的随机种子
            try:
                l5_num = extract_did_number(u)
            except (ValueError, AttributeError) as e:
                # 如果提取失败，使用hash作为fallback
                l5_num = hash(u) & 0x7FFFFFFF
                logger.warning(f"无法从L5_did提取数字，使用hash: {u}, {e}")
            
            # 生成独立的随机种子
            tu_seed = seed_tu(base_seed, target_num, l5_num)
            
            # 使用独立的Random实例进行采样
            rng = random.Random(tu_seed)
            sample_size = min(D, len(b_u))
            sampled = rng.sample(b_u, sample_size)
            
            candidates.update(sampled)
            common_expansion_count += len(sampled)
    
    # Step 4: Compute distances (if enabled)
    candidates_list = sorted(candidates)
    candidates_with_distances = None
    
    if compute_distances:
        distances = {}
        for cand_did in candidates_list:
            try:
                # Compute weighted Jaccard distance
                l5_1 = set(l3_l5_map.get(cand_did, []))
                l5_2 = set(l3_l5_map.get(target_did, []))
                
                if not l5_1 and not l5_2:
                    similarity = 0.0
                else:
                    intersection = l5_1 & l5_2
                    union = l5_1 | l5_2
                    intersection_weight = sum(l5_weight_map.get(u, 0.0) for u in sorted(intersection))
                    union_weight = sum(l5_weight_map.get(u, 0.0) for u in sorted(union))
                    similarity = intersection_weight / union_weight if union_weight > 0 else 0.0
                
                distances[cand_did] = 1.0 - similarity
            except Exception:
                distances[cand_did] = 1.0
        
        candidates_with_distances = [
            {"did": cand_did, "distance": float(distances.get(cand_did, 1.0))}
            for cand_did in candidates_list
        ]
        candidates_with_distances.sort(key=lambda x: (-x["distance"], x["did"]))
    
    # Build result
    result = {
        "target_did": target_did,
        "l5_generating_set": g_t_list,
        "rare_l5": g_t_rare,
        "common_l5": g_t_common,
        "candidates": candidates_list,
        "statistics": {
            "total_candidates": len(candidates),
            "rare_expansion_count": rare_expansion_count,
            "common_expansion_count": common_expansion_count,
            "rare_l5_count": len(g_t_rare),
            "common_l5_count": len(g_t_common),
        }
    }
    
    if compute_distances and candidates_with_distances:
        result["candidates_with_distances"] = candidates_with_distances
        if candidates_with_distances:
            distances_only = [item["distance"] for item in candidates_with_distances]
            result["statistics"]["distance_stats"] = {
                "min_distance": float(min(distances_only)),
                "max_distance": float(max(distances_only)),
                "mean_distance": float(sum(distances_only) / len(distances_only))
            }
    
    return result


# KL/NL相关函数
def calculate_normalized_distance(distance: float, q10: float, q90: float) -> float:
    """计算归一化距离 d_tilde = clip((d - q10) / (q90 - q10), 0, 1)"""
    if q90 == q10:
        return 0.0
    normalized = (distance - q10) / (q90 - q10)
    return max(0.0, min(1.0, normalized))


def weight_function(d_tilde: float, alpha: float, lambda_param: float) -> float:
    """权重函数 w(d;α,λ) = λexp(-αd_tilde) + (1-λ)"""
    return lambda_param * math.exp(-alpha * d_tilde) + (1 - lambda_param)


def weighted_sampling_without_replacement(
    items: List[Tuple[str, float, float]],
    weights: List[float],
    k: int,
    seed: Optional[int] = None,
) -> List[Tuple[str, float, float]]:
    """使用 Efraimidis–Spirakis 方法实现加权不放回抽样。seed 不为 None 时使用该种子保证可复现。"""
    if len(items) <= k:
        return items
    
    rng = random.Random(seed) if seed is not None else random
    keys = []
    for i, weight in enumerate(weights):
        if weight > 0:
            r = rng.random()
            key = math.pow(r, 1.0 / weight)
        else:
            key = 0.0
        keys.append((key, i))
    
    keys.sort(reverse=True)
    selected_indices = [idx for _, idx in keys[:k]]
    return [items[i] for i in selected_indices]


def normalize_did(did: str) -> str:
    """归一化DID：去掉L3_前缀"""
    return did.replace('L3_', '') if did.startswith('L3_') else did


# ---------- 无关 L3 抽样（内联，不引用 tools） ----------
def get_all_l3_set(l5_l3_index: dict) -> Set[str]:
    """从 l5_l3_index 预计算全体 L3 集合。"""
    all_l3: Set[str] = set()
    for l3_list in l5_l3_index.values():
        all_l3.update(l3_list)
    return all_l3


def _is_unrelated_l3(l3_l5_map: dict, l3_did_1: str, l3_did_2: str) -> bool:
    """判断两个 L3 是否在 L5 空间下无关：L5(l1) ∩ L5(l2) = ∅。"""
    n1 = normalize_did(l3_did_1)
    n2 = normalize_did(l3_did_2)
    if n1 == n2:
        return False
    l5_1 = set(l3_l5_map.get(n1, []) or l3_l5_map.get(l3_did_1, []))
    l5_2 = set(l3_l5_map.get(n2, []) or l3_l5_map.get(l3_did_2, []))
    return len(l5_1 & l5_2) == 0


def sample_unrelated_l3_for_target(
    target_did: str,
    k: int,
    l3_l5_map: dict,
    l5_l3_index: dict,
    all_l3_list: List[str],
    base_seed: int,
    target_num: int,
) -> List[str]:
    """
    为给定 target 随机抽取 k 个与其在 L5 空间下无关的 L3。
    每次从全体 L3 中随机抽 128 个（种子 seed_tu(base_seed, target_num, loop_num)，loop_num 从 1 递增），
    从第一个起逐个检查是否无关（并跳过 target 自身），无关则加入，直到凑够 k 个；不够则 loop_num+1 再抽 128 个重复。
    假定总能凑够 k，不处理凑不齐的情况。
    优化：不再每 target 构建 pool，直接对 all_l3_list 抽样，在循环内跳过 target 自身。
    """
    normalized = normalize_did(target_did)
    if not all_l3_list:
        return []
    result: List[str] = []
    loop_num = 1
    batch_size = 128
    n_total = len(all_l3_list)
    while len(result) < k:
        tu_seed = seed_tu(base_seed, target_num, loop_num)
        rng = random.Random(tu_seed)
        batch = rng.sample(all_l3_list, min(batch_size, n_total))
        for l in batch:
            if l == normalized or l == target_did:
                continue
            if _is_unrelated_l3(l3_l5_map, l, target_did) and l not in result:
                result.append(l)
                if len(result) == k:
                    return result
        loop_num += 1
    return result


def build_metal_l3_index_from_neo4j(
    neo4j_uri: str,
    neo4j_database: str,
    neo4j_user: str,
    neo4j_password: str,
    output_path: Path
) -> Dict[str, Set[str]]:
    """
    从Neo4j构建Metal到RepairedLigand的索引并保存为CSV
    
    Args:
        neo4j_uri: Neo4j连接URI
        neo4j_database: 数据库名
        neo4j_user: 用户名
        neo4j_password: 密码
        output_path: 输出CSV文件路径
    
    Returns:
        Dict[metal_symbol, Set[l3_did]]
    """
    if not NEO4J_AVAILABLE:
        raise ImportError("需要安装 neo4j 库来生成索引。安装命令: pip install neo4j")
    
    logger.info("=" * 60)
    logger.info("从Neo4j构建Metal到RepairedLigand索引")
    logger.info("=" * 60)
    
    # 步骤1: 获取所有Metal
    logger.info("步骤1: 获取所有Metal...")
    driver = GraphDatabase.driver(neo4j_uri, auth=(neo4j_user, neo4j_password))
    
    try:
        metals = []
        with driver.session(database=neo4j_database) as session:
            query = "MATCH (m:Metal) RETURN DISTINCT m.symbol AS symbol ORDER BY m.symbol"
            result = session.run(query)
            metals = [record["symbol"] for record in result]
        
        logger.info(f"  找到 {len(metals):,} 个不同的Metal")
        
        # 步骤2: 对每个Metal查询其连接的RepairedLigand
        logger.info("步骤2: 查询每个Metal的RepairedLigand...")
        metal_l3_index = defaultdict(set)
        
        start_time = time.time()
        total_l3_count = 0
        
        with driver.session(database=neo4j_database) as session:
            for i, metal_symbol in enumerate(metals, 1):
                query = """
                MATCH (:Metal {symbol: $metal_symbol})-[:M_L1]->(:Complex)-[:L1_L3]->(l3:RepairedLigand)
                RETURN DISTINCT l3.did AS l3_did
                """
                
                result = session.run(query, metal_symbol=metal_symbol)
                l3_dids = [record["l3_did"] for record in result]
                
                metal_l3_index[metal_symbol] = set(l3_dids)
                total_l3_count += len(l3_dids)
                
                if i % 10 == 0 or i == len(metals):
                    elapsed = time.time() - start_time
                    avg_time = elapsed / i
                    remaining = (len(metals) - i) * avg_time
                    logger.info(f"  进度: {i}/{len(metals)} ({i/len(metals)*100:.1f}%) | "
                              f"当前Metal: {metal_symbol} ({len(l3_dids):,} 个L3) | "
                              f"已用时: {elapsed:.1f}s | 预计剩余: {remaining:.1f}s")
        
        total_time = time.time() - start_time
        logger.info(f"  完成！总耗时: {total_time:.2f} 秒")
        logger.info(f"  总L3连接数: {total_l3_count:,}")
        logger.info(f"  平均每个Metal: {total_l3_count/len(metals):.1f} 个L3")
        
        # 保存为CSV
        logger.info(f"保存索引到: {output_path}")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w', encoding='utf-8', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['metal_symbol', 'l3_did'])
            
            for metal_symbol, l3_dids in sorted(metal_l3_index.items()):
                for l3_did in sorted(l3_dids):
                    writer.writerow([metal_symbol, l3_did])
        
        logger.info(f"  ✓ 已保存 {sum(len(dids) for dids in metal_l3_index.values()):,} 行")
        
        return dict(metal_l3_index)
    
    finally:
        driver.close()


def load_metal_l3_index(index_path: Path) -> Dict[str, Set[str]]:
    """
    加载Metal到RepairedLigand的索引，并归一化所有DID（去掉L3_前缀）
    返回的索引中，所有DID都是归一化后的格式（无L3_前缀）
    """
    suffix = index_path.suffix.lower()
    
    if suffix == '.csv':
        metal_l3_index = defaultdict(set)
        with open(index_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                metal_symbol = row['metal_symbol']
                l3_did = row['l3_did']
                normalized_did = normalize_did(l3_did)
                metal_l3_index[metal_symbol].add(normalized_did)
        return dict(metal_l3_index)
    elif suffix == '.json':
        with open(index_path, 'r', encoding='utf-8') as f:
            json_data = json.load(f)
        return {
            metal_symbol: {normalize_did(did) for did in l3_dids}
            for metal_symbol, l3_dids in json_data.items()
        }
    else:
        raise ValueError(f"不支持的索引文件格式: {suffix}")


def batch_query_kl_proto_from_index(
    metal_l3_index: Dict[str, Set[str]],
    metal_symbol: str,
    candidate_dids: List[str]
) -> Set[str]:
    """
    从索引中批量查询cand(T)中哪些candidate属于A(M)
    
    注意：metal_l3_index 和 candidate_dids 都应该是归一化的（无L3_前缀）
    """
    if metal_symbol not in metal_l3_index:
        return set()
    
    a_m_set = metal_l3_index[metal_symbol]
    
    # candidate_dids 应该已经是归一化的，直接做 O(1) membership 检查
    return {did for did in candidate_dids if did in a_m_set}


def build_target_metals_map(m_l3_pairs_csv: Path) -> Dict[str, List[str]]:
    """
    预处理m_l3_pairs.csv，建立target_did_normalized -> metals(list)的映射
    返回的字典中，key是归一化的target_did（无L3_前缀），value是metal列表
    """
    target_metals_map = defaultdict(set)
    
    with open(m_l3_pairs_csv, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)  # 跳过表头
        
        for row in reader:
            if len(row) >= 4:
                ligand_did = row[0].strip()
                metal_symbol = row[3].strip()
                
                if metal_symbol:
                    ligand_did_normalized = normalize_did(ligand_did)
                    target_metals_map[ligand_did_normalized].add(metal_symbol)
    
    # 转换为列表并排序
    return {
        target_did: sorted(list(metals))
        for target_did, metals in target_metals_map.items()
    }


def process_kl_nl_for_target_online(
    target_did: str,
    candidates_with_distances: List[Dict],  # [{"did": "...", "distance": ...}, ...]
    metals: List[str],
    metal_l3_index: Dict[str, Set[str]],
    kl_nl_writer: csv.DictWriter,
    kl_nl_stats_writer: csv.DictWriter,
    l3_l5_map: Optional[dict] = None,
    l5_l3_index: Optional[dict] = None,
    all_l3_list: Optional[List[str]] = None,
    num_unrelated: int = 32,
    base_seed: int = 42,
    target_num: Optional[int] = None,
):
    """
    在线处理单个target的KL/NL生成，直接写入文件。
    若传入 l3_l5_map / l5_l3_index / all_l3_list，则额外为该 target 生成 num_unrelated 个无关 NL（label='UL'）。

    Args:
        target_did: target的DID（原始格式，可能带L3_前缀）
        candidates_with_distances: candidates列表，每个元素包含did和distance
        metals: 该target对应的metals列表
        metal_l3_index: 归一化的Metal到RepairedLigand索引
        kl_nl_writer: KL/NL结果CSV写入器
        kl_nl_stats_writer: KL/NL统计CSV写入器
        l3_l5_map, l5_l3_index, all_l3_list: 用于生成无关 NL
        num_unrelated: 每个 target 写入的无关 NL 数量，默认 32
        base_seed: 基础随机种子（用于 seed_tu）
        target_num: target_did 的数字部分，用于 seed_tu；不传则从 target_did 解析
    """
    if not candidates_with_distances or not metals:
        return
    
    # 归一化target_did用于查找metals
    target_did_normalized = normalize_did(target_did)
    
    # 提取距离并计算q10和q90
    distances = [item.get('distance', 0.0) for item in candidates_with_distances]
    if not distances or all(d == 0.0 for d in distances):
        # 如果没有有效距离，跳过KL/NL生成
        return
    
    q10 = np.percentile(distances, 10)
    q90 = np.percentile(distances, 90)
    
    # 归一化所有candidate_dids
    candidates_data_normalized = []
    for item in candidates_with_distances:
        candidate_did = item.get('did', '')
        distance = item.get('distance', 0.0)
        candidate_did_normalized = normalize_did(candidate_did)
        d_tilde = calculate_normalized_distance(distance, q10, q90)
        candidates_data_normalized.append((candidate_did_normalized, distance, d_tilde))
    
    # 处理每个metal
    for metal_symbol in metals:
        candidate_dids_normalized = [did for did, _, _ in candidates_data_normalized]
        
        # 从索引查询KL_proto（O(1) membership检查），并显式排除 target 自身
        kl_proto_set = batch_query_kl_proto_from_index(
            metal_l3_index, metal_symbol, candidate_dids_normalized
        )
        kl_proto_set.discard(target_did_normalized)  # 排除 T，避免信息泄露

        # 划分KL_proto和NL_proto（KL 中已不含 T）
        kl_proto = []
        nl_proto = []
        
        for candidate_did_norm, distance, d_tilde in candidates_data_normalized:
            if candidate_did_norm in kl_proto_set:
                kl_proto.append((candidate_did_norm, distance, d_tilde))
            else:
                nl_proto.append((candidate_did_norm, distance, d_tilde))
        
        # 计算权重并抽样（按 (T,M) 确定性种子，与 step13_with_kl_nl 一致）
        _target_num = target_num if target_num is not None else (hash(target_did) & 0x7FFFFFFF)
        metal_hash = sum(ord(c) for c in metal_symbol) & 0x7FFFFFFF
        seed_kl = seed_tu(base_seed, _target_num, metal_hash)
        seed_nl = seed_tu(base_seed, _target_num, metal_hash + 10000)
        if len(kl_proto) > 0:
            kl_weights = [weight_function(d_tilde, alpha=4, lambda_param=1) 
                         for _, _, d_tilde in kl_proto]
            kl_sampled = weighted_sampling_without_replacement(kl_proto, kl_weights, k=16, seed=seed_kl)
        else:
            kl_sampled = []
        
        if len(nl_proto) > 0:
            nl_weights = [weight_function(d_tilde, alpha=2, lambda_param=0.7) 
                         for _, _, d_tilde in nl_proto]
            nl_sampled = weighted_sampling_without_replacement(nl_proto, nl_weights, k=128, seed=seed_nl)
        else:
            nl_sampled = []
        
        # 统计信息
        if kl_sampled:
            kl_distances = [d for _, d, _ in kl_sampled]
            kl_stats = {
                'min': np.min(kl_distances),
                'median': np.median(kl_distances),
                'mean': np.mean(kl_distances),
                'q10': np.percentile(kl_distances, 10),
                'q90': np.percentile(kl_distances, 90),
            }
        else:
            kl_stats = {}
        
        if nl_sampled:
            nl_distances = [d for _, d, _ in nl_sampled]
            nl_stats = {
                'min': np.min(nl_distances),
                'median': np.median(nl_distances),
                'mean': np.mean(nl_distances),
                'q10': np.percentile(nl_distances, 10),
                'q90': np.percentile(nl_distances, 90),
            }
        else:
            nl_stats = {}
        
        # 直接写入KL/NL结果（使用原始target_did格式）
        for candidate_did_norm, distance, d_tilde in kl_sampled:
            weight = weight_function(d_tilde, alpha=4, lambda_param=1)
            kl_nl_writer.writerow({
                'T': target_did,  # 使用原始格式
                'M': metal_symbol,
                'label': 'KL',
                'candidate_did': candidate_did_norm,  # 归一化后的格式
                'distance': distance,
                'd_tilde': d_tilde,
                'weight': weight
            })
        
        for candidate_did_norm, distance, d_tilde in nl_sampled:
            weight = weight_function(d_tilde, alpha=2, lambda_param=0.7)
            kl_nl_writer.writerow({
                'T': target_did,  # 使用原始格式
                'M': metal_symbol,
                'label': 'NL',
                'candidate_did': candidate_did_norm,  # 归一化后的格式
                'distance': distance,
                'd_tilde': d_tilde,
                'weight': weight
            })
        
        # 写入统计摘要
        kl_nl_stats_writer.writerow({
            'T': target_did,
            'M': metal_symbol,
            '|cand(T)|': len(candidates_with_distances),
            'q10': q10,
            'q90': q90,
            '|KL_proto|': len(kl_proto),
            '|KL_proto|/|cand(T)|': len(kl_proto) / len(candidates_with_distances) if candidates_with_distances else 0,
            '|NL_proto|': len(nl_proto),
            '|NL_proto|/|cand(T)|': len(nl_proto) / len(candidates_with_distances) if candidates_with_distances else 0,
            'KL_min': kl_stats.get('min', None),
            'KL_median': kl_stats.get('median', None),
            'KL_mean': kl_stats.get('mean', None),
            'KL_q10': kl_stats.get('q10', None),
            'KL_q90': kl_stats.get('q90', None),
            'NL_min': nl_stats.get('min', None),
            'NL_median': nl_stats.get('median', None),
            'NL_mean': nl_stats.get('mean', None),
            'NL_q10': nl_stats.get('q10', None),
            'NL_q90': nl_stats.get('q90', None),
        })

    # 无关 NL（UL）：每个 target 随机选取 num_unrelated 个 L5 无关 L3 并写入
    if num_unrelated > 0 and l3_l5_map is not None and l5_l3_index is not None and all_l3_list is not None:
        try:
            _target_num = target_num if target_num is not None else extract_did_number(target_did)
            ul_sampled = sample_unrelated_l3_for_target(
                target_did=target_did,
                k=num_unrelated,
                l3_l5_map=l3_l5_map,
                l5_l3_index=l5_l3_index,
                all_l3_list=all_l3_list,
                base_seed=base_seed,
                target_num=_target_num,
            )
            for ul_did in ul_sampled:
                kl_nl_writer.writerow({
                    'T': target_did,
                    'M': '',
                    'label': 'UL',
                    'candidate_did': ul_did,
                    'distance': '',
                    'd_tilde': '',
                    'weight': '',
                })
        except Exception as e:
            logger.debug("无关NL抽样跳过 target %s: %s", target_did, e)


class Step13Processor:
    def __init__(self, record_limit: int = DEFAULT_RECORD_LIMIT,
                 num_samples: int = DEFAULT_NUM_SAMPLES,
                 D: int = DEFAULT_D, K: int = DEFAULT_K,
                 seed: int = DEFAULT_SEED,
                 compute_distances: bool = True,
                 num_workers: Optional[int] = None,
                 generate_kl_nl: bool = False,
                 metal_l3_index_path: Optional[Path] = None,
                 skip_candidates_output: bool = False,
                 num_unrelated: int = 32):
        self.input_dir = Path(INPUT_DIR)
        self.output_dir = Path(OUTPUT_DIR)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.record_limit = record_limit
        self.num_samples = num_samples
        self.D = D
        self.K = K
        self.seed = seed
        self.compute_distances = compute_distances
        self.num_workers = num_workers if num_workers else get_optimal_worker_count()
        self.skip_candidates_output = skip_candidates_output
        
        # KL/NL生成：只有在compute_distances=True时才启用
        self.generate_kl_nl = generate_kl_nl and compute_distances
        if generate_kl_nl and not compute_distances:
            logger.warning("KL/NL生成需要距离计算，但--no-distances已启用，KL/NL生成将被禁用")
        
        # 如果跳过候选集输出，必须启用KL/NL生成
        if skip_candidates_output and not self.generate_kl_nl:
            raise ValueError("--kl-nl-only 选项需要同时启用 --generate-kl-nl")
        
        self.metal_l3_index_path = metal_l3_index_path
        self.num_unrelated = num_unrelated  # 每个 target 生成的无关 NL 数量
        
        # KL/NL相关数据（在初始化时预处理）
        self.metal_l3_index: Optional[Dict[str, Set[str]]] = None
        self.target_metals_map: Optional[Dict[str, List[str]]] = None
        self.auto_generate_index = False
        self.neo4j_uri = None
        self.neo4j_database = None
        self.neo4j_user = None
        self.neo4j_password = None
        
        # 如果启用KL/NL，预处理索引和target-metals映射
        # 注意：索引加载在process_all中延迟进行，以便可以使用命令行参数设置的Neo4j配置
        
        self.stats = {
            'total_targets': 0,
            'total_candidates': 0,
            'success_count': 0,
            'fail_count': 0,
            'processing_time': 0,
            'extract_time': 0,
            'index_load_time': 0,
            'build_time': 0,
            'record_limit': self.record_limit,
            'num_samples': self.num_samples,
            'D': self.D,
            'K': self.K,
            'seed': self.seed
        }
        
        # Input files
        self.m_l3_pairs_input = self.input_dir / "m_l3_pairs.csv"
        
        # Output files
        self.targets_output = self.output_dir / "step13_targets.csv"
        self.candidates_output = self.output_dir / "step13_candidates.csv"
        self.stats_output = self.output_dir / "step13_stats.json"
        
        # KL/NL输出文件
        if self.generate_kl_nl:
            self.kl_nl_output = self.output_dir / "step13_kl_nl_samples.csv"
            self.kl_nl_stats_output = self.output_dir / "step13_kl_nl_samples.stats.csv"
        
        logger.info(f"Step13处理器初始化完成，输出目录: {self.output_dir}")
        logger.info(f"样本数量: {self.num_samples}")
        logger.info(f"参数: D={self.D}, K={self.K}, seed={self.seed}")
        logger.info(f"计算距离: {self.compute_distances}")
        logger.info(f"并行工作进程数: {self.num_workers}")
        if self.skip_candidates_output:
            logger.info(f"候选集输出: 已禁用（仅输出KL/NL数据）")
        if self.generate_kl_nl:
            logger.info(f"KL/NL生成: 已启用（每个 target 无关 NL 数: {self.num_unrelated}）")
            if self.metal_l3_index_path:
                logger.info(f"  索引文件: {self.metal_l3_index_path}")
    
    def _save_stats(self):
        """Save statistics to JSON file"""
        stats_data = {
            'timestamp': datetime.now().isoformat(),
            'processing_time_seconds': self.stats['processing_time'],
            'statistics': self.stats,
            'output_files': {
                'targets': str(self.targets_output),
                'candidates': str(self.candidates_output)
            }
        }
        
        with open(self.stats_output, 'w', encoding='utf-8') as f:
            json.dump(stats_data, f, indent=2, ensure_ascii=False)
        
        logger.info(f"统计信息已保存到: {self.stats_output}")
    
    def extract_unique_dids(self) -> List[str]:
        """
        从m_l3_pairs.csv提取唯一的ligand_did并转换为L3格式
        """
        start_time = time.time()
        logger.info(f"读取M-L3对文件: {self.m_l3_pairs_input}")
        
        if not self.m_l3_pairs_input.exists():
            raise FileNotFoundError(f"M-L3对文件不存在: {self.m_l3_pairs_input}")
        
        unique_dids = set()
        
        with open(self.m_l3_pairs_input, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in tqdm(reader, desc="提取DID", leave=False):
                ligand_did = row.get('ligand_did', '').strip()
                if ligand_did:
                    # Convert to L3 format: D... -> L3_D...
                    if not ligand_did.startswith('L3_'):
                        ligand_did = f"L3_{ligand_did}"
                    unique_dids.add(ligand_did)
        
        extract_time = time.time() - start_time
        self.stats['extract_time'] = extract_time
        
        # Sort for determinism
        sorted_dids = sorted(unique_dids)
        
        logger.info(f"找到 {len(sorted_dids)} 个唯一DID，耗时 {extract_time:.2f}秒")
        
        return sorted_dids
    
    def _load_indices(self) -> dict:
        """Load all indices directly from JSON files and return as shared data dict"""
        logger.info("加载预计算索引...")
        load_start = time.time()
        
        try:
            data_dir = self.input_dir
            
            # Load L3-L5 map
            l3_l5_file = data_dir / "l3_l5.json"
            with open(l3_l5_file, 'r', encoding='utf-8') as f:
                l3_l5_map = json.load(f)
            logger.info(f"加载L3-L5映射: {len(l3_l5_map)} 个L3节点")
            
            # Load L5-L3 inverted index
            l5_l3_file = data_dir / "l5_l3.json"
            with open(l5_l3_file, 'r', encoding='utf-8') as f:
                l5_l3_index = json.load(f)
            logger.info(f"加载L5-L3倒排索引: {len(l5_l3_index)} 个L5单元")
            
            # Load L5 stats
            l5_stats_file = data_dir / "l5_freq_weight.json"
            with open(l5_stats_file, 'r', encoding='utf-8') as f:
                l5_stats = json.load(f)
            logger.info(f"加载L5统计信息: {len(l5_stats)} 个L5单元")
            
            # Load L3 GAC map
            l3_gac_file = data_dir / "l3_gac.json"
            with open(l3_gac_file, 'r', encoding='utf-8') as f:
                l3_gac_map = json.load(f)
            logger.info(f"加载L3 GAC映射: {len(l3_gac_map)} 个L3节点")
            
            # Try to load filtered index
            l5_l3_index_filtered = None
            filtered_file = data_dir / f"l5_l3_filtered_K{self.K}.json"
            if filtered_file.exists():
                with open(filtered_file, 'r', encoding='utf-8') as f:
                    l5_l3_index_filtered = json.load(f)
                logger.info(f"加载GAC过滤索引 (K={self.K}): {len(l5_l3_index_filtered)} 个L5单元")
            
            # Extract weight map from stats
            l5_weight_map = {
                l5_did: stats["weight"]
                for l5_did, stats in l5_stats.items()
            }
            
            load_time = time.time() - load_start
            self.stats['index_load_time'] = load_time
            logger.info(f"索引加载完成，耗时 {load_time:.2f}秒")
            
            if l5_l3_index_filtered is not None:
                logger.info(f"✓ 已加载预过滤索引 (K={self.K})")
            else:
                logger.warning(f"⚠ 未加载预过滤索引，将使用运行时过滤")
            
            return {
                'l3_l5_map': l3_l5_map,
                'l5_l3_index': l5_l3_index,
                'l5_l3_index_filtered': l5_l3_index_filtered,
                'l5_stats': l5_stats,
                'l3_gac_map': l3_gac_map,
                'l5_weight_map': l5_weight_map,
            }
        except Exception as e:
            logger.error(f"加载索引失败: {e}")
            raise
    
    def _check_input_files(self) -> bool:
        """检查所需的输入文件是否存在"""
        missing_files = []
        
        # 必需文件
        required_files = [
            ("m_l3_pairs.csv", self.m_l3_pairs_input),
            ("l3_l5.json", self.input_dir / "l3_l5.json"),
            ("l5_l3.json", self.input_dir / "l5_l3.json"),
            ("l5_freq_weight.json", self.input_dir / "l5_freq_weight.json"),
            ("l3_gac.json", self.input_dir / "l3_gac.json"),
        ]
        
        for name, path in required_files:
            if not path.exists():
                missing_files.append((name, path))
        
        if missing_files:
            logger.error("=" * 60)
            logger.error("缺少必需的输入文件:")
            for name, path in missing_files:
                logger.error(f"  - {name}: {path}")
            logger.error("=" * 60)
            return False
        
        # 可选文件
        optional_file = self.input_dir / f"l5_l3_filtered_K{self.K}.json"
        if optional_file.exists():
            logger.info(f"找到可选文件: l5_l3_filtered_K{self.K}.json")
        else:
            logger.info(f"未找到可选文件: l5_l3_filtered_K{self.K}.json（将使用运行时过滤）")
        
        return True
    
    def process_all(self) -> Dict:
        """执行所有处理步骤"""
        start_time = time.time()
        logger.info("开始Step13候选集构建")
        logger.info("=" * 60)
        
        # 检查输入文件
        logger.info("检查输入文件...")
        if not self._check_input_files():
            raise FileNotFoundError("缺少必需的输入文件，请检查输入目录")
        logger.info("  ✓ 所有必需文件存在")
        
        # 如果启用KL/NL，先处理索引加载或生成
        if self.generate_kl_nl:
            # 如果索引文件不存在，尝试自动生成
            if not self.metal_l3_index_path or not self.metal_l3_index_path.exists():
                logger.warning(f"索引文件不存在: {self.metal_l3_index_path}")
                logger.info("尝试自动生成索引...")
                
                if not NEO4J_AVAILABLE:
                    logger.error("无法自动生成索引：neo4j库未安装")
                    logger.error("请先运行 build_metal_l3_index.py 生成索引，或安装neo4j库")
                    logger.error("安装命令: pip install neo4j")
                    logger.error("KL/NL生成被禁用")
                    self.generate_kl_nl = False
                elif not self.auto_generate_index:
                    logger.error("无法自动生成索引：未提供Neo4j配置")
                    logger.error("请使用 --metal-l3-index 指定索引文件，或提供Neo4j配置参数")
                    logger.error("KL/NL生成被禁用")
                    self.generate_kl_nl = False
                else:
                    # 使用默认路径（output_dir）
                    auto_index_path = self.output_dir / "metal_l3_index.csv"
                    logger.info(f"将在 {auto_index_path} 生成索引")
                    
                    try:
                        self.metal_l3_index = build_metal_l3_index_from_neo4j(
                            neo4j_uri=self.neo4j_uri,
                            neo4j_database=self.neo4j_database,
                            neo4j_user=self.neo4j_user,
                            neo4j_password=self.neo4j_password,
                            output_path=auto_index_path
                        )
                        self.metal_l3_index_path = auto_index_path
                        logger.info("  ✓ 索引生成成功")
                    except Exception as e:
                        logger.error(f"索引生成失败: {e}")
                        logger.error("KL/NL生成被禁用")
                        self.generate_kl_nl = False
            
            if self.generate_kl_nl:
                # 加载并归一化索引
                logger.info(f"加载Metal到RepairedLigand索引: {self.metal_l3_index_path}")
                load_start = time.time()
                self.metal_l3_index = load_metal_l3_index(self.metal_l3_index_path)
                load_time = time.time() - load_start
                logger.info(f"  加载了 {len(self.metal_l3_index):,} 个Metal的索引，耗时 {load_time:.2f} 秒")
                
                # 预处理m_l3_pairs.csv
                logger.info("预处理m_l3_pairs.csv，建立target-metals映射...")
                preprocess_start = time.time()
                self.target_metals_map = build_target_metals_map(self.m_l3_pairs_input)
                preprocess_time = time.time() - preprocess_start
                logger.info(f"  建立了 {len(self.target_metals_map):,} 个target的metals映射，耗时 {preprocess_time:.2f} 秒")
        
        # Note: Random seed is now handled per (target_did, L5_did) pair
        # using seed_tu() function for independent and reproducible sampling
        logger.info(f"基础随机种子: {self.seed} (将用于seed_tu函数生成每个(target, L5)对的独立种子)")
        
        # Step 1: Extract unique DIDs
        unique_dids = self.extract_unique_dids()
        
        # Step 2: Sample DIDs（使用 self.seed 保证可复现，与 step13_with_kl_nl 同种子则同 target 集合）
        if self.num_samples > 0 and len(unique_dids) > self.num_samples:
            rng = random.Random(self.seed)
            selected_dids = rng.sample(unique_dids, self.num_samples)
            logger.info(f"随机抽取 {self.num_samples} 个DID (seed={self.seed})")
        else:
            selected_dids = unique_dids
            logger.info(f"使用所有 {len(selected_dids)} 个DID")
        
        self.stats['total_targets'] = len(selected_dids)
        
        # Step 3: Load data into module-level globals BEFORE Pool creation
        # This allows workers to inherit via fork's copy-on-write (no pickle!)
        global _SHARED_L3_L5_MAP, _SHARED_L5_L3_INDEX, _SHARED_L5_L3_INDEX_FILTERED
        global _SHARED_L5_STATS, _SHARED_L3_GAC_MAP, _SHARED_L5_WEIGHT_MAP
        
        shared_data = self._load_indices()
        _SHARED_L3_L5_MAP = shared_data['l3_l5_map']
        _SHARED_L5_L3_INDEX = shared_data['l5_l3_index']
        _SHARED_L5_L3_INDEX_FILTERED = shared_data['l5_l3_index_filtered']
        _SHARED_L5_STATS = shared_data['l5_stats']
        _SHARED_L3_GAC_MAP = shared_data['l3_gac_map']
        _SHARED_L5_WEIGHT_MAP = shared_data['l5_weight_map']
        del shared_data  # Free the dict wrapper
        
        # 预计算全体 L3 列表（仅 KL/NL 时），供每个 target 生成无关 NL 用
        all_l3_list = None
        if self.generate_kl_nl:
            all_l3_list = list(get_all_l3_set(_SHARED_L5_L3_INDEX))
        
        # Step 4: Open output CSV files
        logger.info("=" * 60)
        logger.info("开始批量构建候选集...")
        logger.info(f"D={self.D}, K={self.K}, seed={self.seed}")
        logger.info("=" * 60)
        
        # 打开候选集输出文件（如果未禁用）
        targets_file = None
        candidates_file = None
        targets_writer = None
        candidates_writer = None
        
        if not self.skip_candidates_output:
            targets_file = open(self.targets_output, 'w', newline='', encoding='utf-8')
            candidates_file = open(self.candidates_output, 'w', newline='', encoding='utf-8')
            targets_writer = csv.writer(targets_file)
            candidates_writer = csv.writer(candidates_file)
            
            # Write headers
            targets_writer.writerow(['target_did', 'l5_generating_set', 'rare_l5', 'common_l5', 'total_candidates'])
            candidates_writer.writerow(['target_did', 'candidate_did', 'distance'])
        
        # KL/NL输出文件（如果启用）
        kl_nl_file = None
        kl_nl_stats_file = None
        kl_nl_writer = None
        kl_nl_stats_writer = None
        if self.generate_kl_nl:
            kl_nl_file = open(self.kl_nl_output, 'w', newline='', encoding='utf-8')
            kl_nl_stats_file = open(self.kl_nl_stats_output, 'w', newline='', encoding='utf-8')
            kl_nl_writer = csv.DictWriter(kl_nl_file, fieldnames=['T', 'M', 'label', 'candidate_did', 'distance', 'd_tilde', 'weight'])
            kl_nl_stats_writer = csv.DictWriter(kl_nl_stats_file, fieldnames=[
                'T', 'M', '|cand(T)|', 'q10', 'q90', 
                '|KL_proto|', '|KL_proto|/|cand(T)|', '|NL_proto|', '|NL_proto|/|cand(T)|',
                'KL_min', 'KL_median', 'KL_mean', 'KL_q10', 'KL_q90',
                'NL_min', 'NL_median', 'NL_mean', 'NL_q10', 'NL_q90'
            ])
            kl_nl_writer.writeheader()
            kl_nl_stats_writer.writeheader()
        
        build_start = time.time()
        total_candidates = 0
        failed_dids = []  # Collect failures to show at end
        
        try:
            # Process in batches to prevent memory accumulation
            # Restart Pool between batches to release IPC buffers
            BATCH_SIZE = 5000
            num_batches = (len(selected_dids) + BATCH_SIZE - 1) // BATCH_SIZE
            
            with tqdm(total=len(selected_dids), desc="构建候选集", unit="个") as pbar:
                for batch_idx in range(num_batches):
                    batch_start = batch_idx * BATCH_SIZE
                    batch_end = min(batch_start + BATCH_SIZE, len(selected_dids))
                    batch_dids = selected_dids[batch_start:batch_end]
                    
                    # Create new Pool for each batch (releases memory on close)
                    with Pool(
                        processes=self.num_workers,
                        initializer=_init_worker,
                        initargs=(self.D, self.K, self.seed, self.compute_distances)
                    ) as pool:
                        for success, did, result, elapsed in pool.imap_unordered(
                            _process_single_did, batch_dids, chunksize=10
                        ):
                            if success and result:
                                self.stats['success_count'] += 1
                                
                                # Write target row (如果未禁用候选集输出)
                                target_did = result.get('target_did', '')
                                l5_generating_set = ';'.join(result.get('l5_generating_set', []))
                                rare_l5 = ';'.join(result.get('rare_l5', []))
                                common_l5 = ';'.join(result.get('common_l5', []))
                                candidate_count = result.get('statistics', {}).get('total_candidates', 0)
                                
                                if targets_writer is not None:
                                    targets_writer.writerow([
                                        target_did,
                                        l5_generating_set,
                                        rare_l5,
                                        common_l5,
                                        candidate_count
                                    ])
                                
                                # Write candidate rows (如果未禁用候选集输出)
                                candidates_with_distances = result.get('candidates_with_distances', [])
                                if not candidates_with_distances and 'candidates' in result:
                                    candidates_with_distances = [
                                        {'did': cand, 'distance': ''}
                                        for cand in result['candidates']
                                    ]
                                
                                if candidates_writer is not None:
                                    for cand in candidates_with_distances:
                                        candidates_writer.writerow([
                                            target_did,
                                            cand.get('did', ''),
                                            cand.get('distance', '')
                                        ])
                                
                                total_candidates += len(candidates_with_distances)
                                
                                # 在线处理KL/NL（如果启用）
                                if self.generate_kl_nl and candidates_with_distances:
                                    # 获取该target的metals（使用预处理好的映射）
                                    target_did_normalized = normalize_did(target_did)
                                    metals = self.target_metals_map.get(target_did_normalized, [])
                                    
                                    if metals:
                                        target_num_for_ul = None
                                        try:
                                            target_num_for_ul = extract_did_number(target_did)
                                        except (ValueError, AttributeError):
                                            target_num_for_ul = hash(target_did) & 0x7FFFFFFF
                                        process_kl_nl_for_target_online(
                                            target_did=target_did,
                                            candidates_with_distances=candidates_with_distances,
                                            metals=metals,
                                            metal_l3_index=self.metal_l3_index,
                                            kl_nl_writer=kl_nl_writer,
                                            kl_nl_stats_writer=kl_nl_stats_writer,
                                            l3_l5_map=_SHARED_L3_L5_MAP,
                                            l5_l3_index=_SHARED_L5_L3_INDEX,
                                            all_l3_list=all_l3_list,
                                            num_unrelated=self.num_unrelated,
                                            base_seed=self.seed,
                                            target_num=target_num_for_ul,
                                        )
                                
                                del candidates_with_distances
                                del result
                            else:
                                self.stats['fail_count'] += 1
                                if len(failed_dids) < 100:
                                    failed_dids.append((did, str(result)[:200]))
                            
                            pbar.update(1)
                    
                    # Force garbage collection after each batch
                    if targets_file:
                        targets_file.flush()
                    if candidates_file:
                        candidates_file.flush()
                    if kl_nl_file:
                        kl_nl_file.flush()
                    if kl_nl_stats_file:
                        kl_nl_stats_file.flush()
                    gc.collect()
        
        finally:
            if targets_file:
                targets_file.close()
            if candidates_file:
                candidates_file.close()
            if kl_nl_file:
                kl_nl_file.close()
            if kl_nl_stats_file:
                kl_nl_stats_file.close()
        
        build_time = time.time() - build_start
        self.stats['build_time'] = build_time
        self.stats['total_candidates'] = total_candidates
        
        # Calculate total time
        total_time = time.time() - start_time
        self.stats['processing_time'] = total_time
        
        # Show collected warnings at end
        if failed_dids:
            logger.warning(f"以下 {len(failed_dids)} 个DID构建失败:")
            # Show first 10 failures
            for did, error_msg in failed_dids[:10]:
                logger.warning(f"  - {did}: {error_msg}")
            if len(failed_dids) > 10:
                logger.warning(f"  ... 还有 {len(failed_dids) - 10} 个失败")
        
        # Save statistics
        self._save_stats()
        
        logger.info("=" * 60)
        logger.info("Step13候选集构建完成！统计信息:")
        logger.info(f"  目标数: {self.stats['total_targets']}")
        logger.info(f"  成功: {self.stats['success_count']}")
        logger.info(f"  失败: {self.stats['fail_count']}")
        logger.info(f"  总候选数: {self.stats['total_candidates']}")
        logger.info(f"  提取时间: {self.stats['extract_time']:.2f}秒")
        logger.info(f"  索引加载时间: {self.stats['index_load_time']:.2f}秒")
        logger.info(f"  构建时间: {self.stats['build_time']:.2f}秒")
        logger.info(f"  总处理时间: {total_time:.2f}秒")
        if self.stats['success_count'] > 0:
            avg_time = build_time / self.stats['success_count']
            logger.info(f"  平均构建时间: {avg_time:.3f}秒/个")
        logger.info(f"输出文件:")
        if not self.skip_candidates_output:
            logger.info(f"  Targets: {self.targets_output}")
            logger.info(f"  Candidates: {self.candidates_output}")
        if self.generate_kl_nl:
            logger.info(f"  KL/NL Samples: {self.kl_nl_output}")
            logger.info(f"  KL/NL Stats: {self.kl_nl_stats_output}")
        logger.info("=" * 60)
        
        return self.stats


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Step13: 批量构建候选集")
    parser.add_argument('--num-samples', type=int, default=DEFAULT_NUM_SAMPLES,
                       help=f'处理的样本数量 (默认: {DEFAULT_NUM_SAMPLES}, -1 = 全部)')
    parser.add_argument('--D', type=int, default=DEFAULT_D,
                       help=f'稀有/常见L5阈值 (默认: {DEFAULT_D})')
    parser.add_argument('--K', type=int, default=DEFAULT_K,
                       help=f'GAC截断值 (默认: {DEFAULT_K})')
    parser.add_argument('--seed', type=int, default=DEFAULT_SEED, help='随机种子')
    parser.add_argument('--no-distances', action='store_true',
                       help='不计算候选集距离（加快速度）；启用时会同时关闭 KL/NL 生成')
    parser.add_argument('--workers', type=int, default=None, help='并行工作进程数（默认使用所有CPU核心）')
    parser.add_argument('--quiet', action='store_true', help='静默模式，只输出关键信息')
    parser.add_argument('--input-dir', type=str, default=INPUT_DIR,
                       help=f'输入目录 (默认: {INPUT_DIR})')
    parser.add_argument('--output-dir', type=str, default=OUTPUT_DIR,
                       help=f'输出目录 (默认: {OUTPUT_DIR})')
    
    # KL/NL相关参数（默认开启）
    parser.add_argument('--generate-kl-nl', dest='generate_kl_nl', action='store_true', default=True,
                       help='生成 KL/NL 样本（默认开启，需要 metal_l3_index）')
    parser.add_argument('--no-generate-kl-nl', dest='generate_kl_nl', action='store_false',
                       help='关闭 KL/NL 样本生成')
    parser.add_argument('--metal-l3-index', type=str, default=None,
                       help='Metal到RepairedLigand索引文件路径（默认: output_dir/metal_l3_index.csv，如果不存在会自动生成）')
    parser.add_argument('--kl-nl-only', dest='kl_nl_only', action='store_true', default=True,
                       help='仅输出 KL/NL 数据，不输出候选集文件（默认开启；需 --generate-kl-nl 开启时有效）')
    parser.add_argument('--no-kl-nl-only', dest='kl_nl_only', action='store_false',
                       help='同时输出候选集文件与 KL/NL 数据')
    parser.add_argument('--num-unrelated', type=int, default=32,
                       help='每个 target 生成的无关 NL（label=UL）数量（默认 32）')
    
    # Neo4j配置参数（用于自动生成索引）
    parser.add_argument('--neo4j-uri', type=str, default=None,
                       help='Neo4j连接URI（默认: 从环境变量NEO4J_URI获取，或bolt://localhost:3087）')
    parser.add_argument('--neo4j-database', type=str, default=None,
                       help='Neo4j数据库名称（默认: 从环境变量NEO4J_DATABASE获取，或neo4j3）')
    parser.add_argument('--neo4j-user', type=str, default=None,
                       help='Neo4j用户名（默认: 从环境变量NEO4J_USER获取，或neo4j）')
    parser.add_argument('--neo4j-password', type=str, default=None,
                       help='Neo4j密码（默认: 从环境变量NEO4J_PASSWORD获取）')
    
    args = parser.parse_args()
    
    # Configure logging level based on quiet mode
    if args.quiet:
        logging.getLogger().setLevel(logging.WARNING)
    
    # Handle num_samples = -1 as unlimited
    num_samples = args.num_samples if args.num_samples > 0 else -1
    
    # kl-nl-only 仅在 generate-kl-nl 开启时有效；关闭 KL/NL 时自动改为输出候选集
    if not args.generate_kl_nl and args.kl_nl_only:
        args.kl_nl_only = False
    
    # 确定索引文件路径
    metal_l3_index_path = None
    if args.generate_kl_nl:
        if args.metal_l3_index:
            metal_l3_index_path = Path(args.metal_l3_index)
        else:
            # 默认路径：优先使用output_dir，其次使用result_analysis目录
            output_dir = Path(args.output_dir)
            default_index = output_dir / "metal_l3_index.csv"
            
            if default_index.exists():
                metal_l3_index_path = default_index
            else:
                metal_l3_index_path = default_index
                logger.info(f"索引文件不存在，将在运行时自动生成: {metal_l3_index_path}")
    
    # Create processor and run
    processor = Step13Processor(
        num_samples=num_samples,
        D=args.D,
        K=args.K,
        seed=args.seed,
        compute_distances=not args.no_distances,
        num_workers=args.workers,
        generate_kl_nl=args.generate_kl_nl,
        metal_l3_index_path=metal_l3_index_path,
        skip_candidates_output=args.kl_nl_only,  # 默认 True：仅输出 KL/NL
        num_unrelated=args.num_unrelated,
    )
    
    # Override directories if specified
    processor.input_dir = Path(args.input_dir)
    processor.output_dir = Path(args.output_dir)
    processor.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Update file paths
    processor.m_l3_pairs_input = processor.input_dir / "m_l3_pairs.csv"
    processor.targets_output = processor.output_dir / "step13_targets.csv"
    processor.candidates_output = processor.output_dir / "step13_candidates.csv"
    processor.stats_output = processor.output_dir / "step13_stats.json"
    
    if args.generate_kl_nl:
        processor.kl_nl_output = processor.output_dir / "step13_kl_nl_samples.csv"
        processor.kl_nl_stats_output = processor.output_dir / "step13_kl_nl_samples.stats.csv"
        
        # 如果索引路径在output_dir中且不存在，设置Neo4j配置用于自动生成
        if metal_l3_index_path and metal_l3_index_path.parent == processor.output_dir and not metal_l3_index_path.exists():
            # 设置Neo4j配置到processor（用于自动生成）
            processor.auto_generate_index = True
            processor.neo4j_uri = args.neo4j_uri or os.getenv("NEO4J_URI", "bolt://localhost:3087")
            processor.neo4j_database = args.neo4j_database or os.getenv("NEO4J_DATABASE", "neo4j3")
            processor.neo4j_user = args.neo4j_user or os.getenv("NEO4J_USER", "neo4j")
            processor.neo4j_password = args.neo4j_password or os.getenv("NEO4J_PASSWORD", "testtest123")
        else:
            processor.auto_generate_index = False
    
    processor.process_all()


if __name__ == "__main__":
    main()
