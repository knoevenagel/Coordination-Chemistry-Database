#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
距离计算模块：实现加权Jaccard相似度和距离
"""

import logging
from typing import Dict, List, Set, Optional

logger = logging.getLogger(__name__)


class DistanceCalculator:
    """距离计算器"""
    
    def __init__(self, l3_l5_map: Dict[str, List[str]], 
                 l5_weight_map: Dict[str, float]):
        """
        初始化距离计算器
        
        Args:
            l3_l5_map: L3 -> L5映射字典
            l5_weight_map: L5 -> 权重映射字典（从l5_stats提取）
        """
        self.l3_l5_map = l3_l5_map
        self.l5_weight_map = l5_weight_map
    
    def compute_weighted_jaccard_similarity(self, l1_did: str, l2_did: str) -> float:
        """
        计算两个配体的加权Jaccard相似度
        
        s_L5(l1, l2) = Σ_{u ∈ L5(l1)∩L5(l2)} w(u) / Σ_{u ∈ L5(l1)∪L5(l2)} w(u)
        
        Args:
            l1_did: 配体1的did
            l2_did: 配体2的did
            
        Returns:
            相似度值（0到1之间）
        """
        # 获取L5集合
        l5_1 = set(self.l3_l5_map.get(l1_did, []))
        l5_2 = set(self.l3_l5_map.get(l2_did, []))
        
        if not l5_1 and not l5_2:
            return 0.0
        
        # 计算交集和并集
        intersection = l5_1 & l5_2
        union = l5_1 | l5_2
        
        # 🔧 FIX: 排序确保确定性的求和顺序（避免浮点数累加顺序导致的微小差异）
        intersection_weight = sum(self.l5_weight_map.get(u, 0.0) for u in sorted(intersection))
        union_weight = sum(self.l5_weight_map.get(u, 0.0) for u in sorted(union))
        
        # 如果分母为0，返回0
        if union_weight == 0:
            return 0.0
        
        similarity = intersection_weight / union_weight
        return similarity
    
    def compute_prior_distance(self, l1_did: str, l2_did: str) -> float:
        """
        计算先验距离
        
        d_L5(l1, l2) = 1 - s_L5(l1, l2)
        
        Args:
            l1_did: 配体1的did
            l2_did: 配体2的did
            
        Returns:
            距离值（0到1之间）
        """
        similarity = self.compute_weighted_jaccard_similarity(l1_did, l2_did)
        distance = 1.0 - similarity
        return distance
    
    def batch_compute_distances(self, cand_list: List[str], target_did: str) -> Dict[str, float]:
        """
        批量计算候选配体与目标配体的距离
        
        Args:
            cand_list: 候选配体did列表
            target_did: 目标配体did
            
        Returns:
            距离字典: {cand_did: distance, ...}
        """
        distances = {}
        for cand_did in cand_list:
            try:
                dist = self.compute_prior_distance(cand_did, target_did)
                distances[cand_did] = dist
            except Exception as e:
                logger.warning(f"计算距离失败 {cand_did} -> {target_did}: {e}")
                distances[cand_did] = 1.0  # 默认最大距离
        
        return distances


def load_distance_calculator(data_dir: str) -> DistanceCalculator:
    """
    从数据目录加载距离计算器
    
    Args:
        data_dir: 数据目录路径
        
    Returns:
        DistanceCalculator实例
    """
    import json
    from pathlib import Path
    
    data_path = Path(data_dir)
    
    # 加载L3-L5映射
    l3_l5_file = data_path / "l3_l5.json"
    with open(l3_l5_file, 'r', encoding='utf-8') as f:
        l3_l5_map = json.load(f)
    
    # 加载L5统计信息并提取权重
    l5_stats_file = data_path / "l5_freq_weight.json"
    with open(l5_stats_file, 'r', encoding='utf-8') as f:
        l5_stats = json.load(f)
    
    l5_weight_map = {l5_did: stats["weight"] for l5_did, stats in l5_stats.items()}
    
    return DistanceCalculator(l3_l5_map, l5_weight_map)
