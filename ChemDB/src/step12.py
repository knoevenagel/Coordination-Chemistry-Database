#!/usr/bin/env python3
"""
ChemDB Step12 Processor - Precomputed Index Builder (No Neo4j Required)
从Step8的CSV输出构建预计算索引，无需Neo4j数据库连接
同时生成m_l3_pairs.csv用于后续候选集构建

输出文件:
- l3_l5.json: L3 → L5 映射
- l5_l3.json: L5 → L3 倒排索引
- l5_freq_weight.json: L5 频率和权重统计
- l3_gac.json: L3 GAC 映射
- l5_l3_filtered_K{K}.json: GAC过滤后的L5→L3索引
- m_l3_pairs.csv: 金属-配体对（用于候选集构建）
"""

# Default configuration variables
DEFAULT_RECORD_LIMIT = -1  # -1 for unlimited
OUTPUT_DIR = "./tmp"
INPUT_DIR = "./tmp"
DEFAULT_K = 30  # Default GAC cutoff
DEFAULT_MIN_GAC = 10  # Minimum GAC for m_l3_pairs
DEFAULT_MAX_GAC = 15  # Maximum GAC for m_l3_pairs (legacy default)

import os
import json
import math
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import pandas as pd
from tqdm import tqdm
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.handlers = []
logger.propagate = False
logger.addHandler(logging.StreamHandler())


class Step12Processor:
    def __init__(self, record_limit: int = DEFAULT_RECORD_LIMIT, K: int = DEFAULT_K,
                 min_gac: int = DEFAULT_MIN_GAC, max_gac: int = DEFAULT_MAX_GAC):
        self.input_dir = Path(INPUT_DIR)
        self.output_dir = Path(OUTPUT_DIR)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.record_limit = record_limit
        self.K = K
        self.min_gac = min_gac
        self.max_gac = max_gac
        
        self.stats = {
            'total_l3_nodes': 0,
            'total_l5_nodes': 0,
            'total_l3_l5_mappings': 0,
            'total_l5_l3_mappings': 0,
            'total_filtered_l5_l3_mappings': 0,
            'total_m_l3_pairs': 0,
            'processing_time': 0,
            'record_limit': self.record_limit,
            'K': self.K
        }
        
        # Output files
        self.l3_l5_output = self.output_dir / "l3_l5.json"
        self.l5_l3_output = self.output_dir / "l5_l3.json"
        self.l5_stats_output = self.output_dir / "l5_freq_weight.json"
        self.l3_gac_output = self.output_dir / "l3_gac.json"
        self.l5_l3_filtered_output = self.output_dir / f"l5_l3_filtered_K{self.K}.json"
        self.m_l3_pairs_output = self.output_dir / "m_l3_pairs.csv"
        self.stats_output = self.output_dir / "step12_stats.json"
        
        logger.info(f"Step12处理器初始化完成，输出目录: {self.output_dir}")
        logger.info(f"记录限制: {self.record_limit if self.record_limit > 0 else '无限制'}")
        logger.info(f"GAC截断值 K: {self.K}")
        logger.info(f"M-L3对GAC范围: [{self.min_gac}, {self.max_gac}]")
    
    def _save_json(self, data: dict, output_file: Path):
        """Save dictionary to JSON file"""
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
        logger.info(f"数据已保存到: {output_file}")
    
    def _save_stats(self):
        """Save statistics to JSON file"""
        stats_data = {
            'timestamp': datetime.now().isoformat(),
            'processing_time_seconds': self.stats['processing_time'],
            'statistics': self.stats,
            'output_files': {
                'l3_l5': str(self.l3_l5_output),
                'l5_l3': str(self.l5_l3_output),
                'l5_stats': str(self.l5_stats_output),
                'l3_gac': str(self.l3_gac_output),
                'l5_l3_filtered': str(self.l5_l3_filtered_output),
                'm_l3_pairs': str(self.m_l3_pairs_output)
            }
        }
        
        with open(self.stats_output, 'w', encoding='utf-8') as f:
            json.dump(stats_data, f, indent=2, ensure_ascii=False)
        
        logger.info(f"统计信息已保存到: {self.stats_output}")
    
    def build_l3_l5_mapping(self) -> Dict[str, List[str]]:
        """
        构建L3→L5映射
        
        通过L3→L4→L5路径从CSV文件构建（复现legacy Neo4j行为）
        """
        logger.info("构建L3→L5映射...")
        
        # Load relationship CSVs
        l3_l4_file = self.input_dir / "neo4j_l3_l4_relationships.csv"
        l4_l5_file = self.input_dir / "neo4j_l4_l5_relationships.csv"
        
        if not l3_l4_file.exists():
            raise FileNotFoundError(f"L3-L4关系文件不存在: {l3_l4_file}")
        if not l4_l5_file.exists():
            raise FileNotFoundError(f"L4-L5关系文件不存在: {l4_l5_file}")
        
        logger.info(f"加载L3-L4关系: {l3_l4_file}")
        l3_l4_df = pd.read_csv(l3_l4_file, dtype=str)
        logger.info(f"  加载了 {len(l3_l4_df)} 条L3-L4关系")
        
        logger.info(f"加载L4-L5关系: {l4_l5_file}")
        l4_l5_df = pd.read_csv(l4_l5_file, dtype=str)
        logger.info(f"  加载了 {len(l4_l5_df)} 条L4-L5关系")
        
        # Get column names
        l3_col = [c for c in l3_l4_df.columns if 'START_ID' in c or 'repaired_ligand' in c.lower()][0]
        l4_col_l3l4 = [c for c in l3_l4_df.columns if 'END_ID' in c or 'fragment' in c.lower()][0]
        l4_col_l4l5 = [c for c in l4_l5_df.columns if 'START_ID' in c or 'fragment' in c.lower()][0]
        l5_col = [c for c in l4_l5_df.columns if 'END_ID' in c or 'irl' in c.lower()][0]
        
        # Build L4→L5 mapping
        logger.info("构建L4→L5映射...")
        l4_to_l5: Dict[str, Set[str]] = {}
        for _, row in tqdm(l4_l5_df.iterrows(), total=len(l4_l5_df), desc="L4→L5"):
            l4_id = str(row[l4_col_l4l5]).strip()
            l5_id = str(row[l5_col]).strip()
            if l4_id and l5_id:
                if l4_id not in l4_to_l5:
                    l4_to_l5[l4_id] = set()
                l4_to_l5[l4_id].add(l5_id)
        
        logger.info(f"L4→L5映射: {len(l4_to_l5)} 个L4节点")
        
        # Build L3→L5 mapping by joining L3→L4 and L4→L5
        logger.info("构建L3→L5映射...")
        l3_l5_map: Dict[str, Set[str]] = {}
        for _, row in tqdm(l3_l4_df.iterrows(), total=len(l3_l4_df), desc="L3→L5"):
            l3_id = str(row[l3_col]).strip()
            l4_id = str(row[l4_col_l3l4]).strip()
            if l3_id and l4_id and l4_id in l4_to_l5:
                if l3_id not in l3_l5_map:
                    l3_l5_map[l3_id] = set()
                l3_l5_map[l3_id].update(l4_to_l5[l4_id])
        
        # Convert sets to sorted lists for JSON serialization
        l3_l5_map_list = {k: sorted(v) for k, v in l3_l5_map.items()}
        
        self.stats['total_l3_nodes'] = len(l3_l5_map_list)
        self.stats['total_l3_l5_mappings'] = sum(len(v) for v in l3_l5_map_list.values())
        
        logger.info(f"L3→L5映射构建完成: {len(l3_l5_map_list)} 个L3节点, {self.stats['total_l3_l5_mappings']} 条映射")
        
        # Save to file
        self._save_json(l3_l5_map_list, self.l3_l5_output)
        
        return l3_l5_map_list
    
    def build_l5_l3_index(self, l3_l5_map: Dict[str, List[str]]) -> Dict[str, List[str]]:
        """
        从L3→L5映射构建L5→L3倒排索引
        """
        logger.info("构建L5→L3倒排索引...")
        
        l5_l3_index: Dict[str, List[str]] = {}
        
        for l3_did, l5_list in tqdm(l3_l5_map.items(), desc="构建倒排索引"):
            for l5_did in l5_list:
                if l5_did not in l5_l3_index:
                    l5_l3_index[l5_did] = []
                l5_l3_index[l5_did].append(l3_did)
        
        # Sort each list for determinism
        for l5_did in l5_l3_index:
            l5_l3_index[l5_did] = sorted(l5_l3_index[l5_did])
        
        self.stats['total_l5_nodes'] = len(l5_l3_index)
        self.stats['total_l5_l3_mappings'] = sum(len(v) for v in l5_l3_index.values())
        
        logger.info(f"L5→L3倒排索引构建完成: {len(l5_l3_index)} 个L5单元")
        
        # Save to file
        self._save_json(l5_l3_index, self.l5_l3_output)
        
        return l5_l3_index
    
    def compute_l5_statistics(self, l5_l3_index: Dict[str, List[str]]) -> Dict[str, Dict[str, float]]:
        """
        计算L5频率和权重
        
        freq(u) = |{l ∈ L3 | u ∈ L5(l)}|
        w(u) = ln(|L3| / freq(u))
        """
        logger.info("计算L5统计信息...")
        
        # Calculate total L3 count
        all_l3_set = set()
        for l3_list in l5_l3_index.values():
            all_l3_set.update(l3_list)
        total_l3_count = len(all_l3_set)
        logger.info(f"总L3数量: {total_l3_count}")
        
        # Calculate frequency and weight
        l5_stats: Dict[str, Dict[str, float]] = {}
        
        for l5_did, l3_list in tqdm(l5_l3_index.items(), desc="计算统计信息"):
            freq = len(l3_list)
            if freq > 0:
                weight = math.log(total_l3_count / freq)
            else:
                weight = 0.0
            l5_stats[l5_did] = {
                "freq": freq,
                "weight": weight
            }
        
        logger.info(f"L5统计信息计算完成: {len(l5_stats)} 个L5单元")
        
        # Save to file
        self._save_json(l5_stats, self.l5_stats_output)
        
        return l5_stats
    
    def extract_l3_gac(self) -> Dict[str, str]:
        """
        从ligand_with_gac.csv提取L3的GAC值
        如果该文件不存在，尝试从neo4j_repaired_ligands.csv提取
        """
        logger.info("提取L3 GAC值...")
        
        # Try ligand_with_gac.csv first (has GAC values)
        gac_file = self.input_dir / "ligand_with_gac.csv"
        if gac_file.exists():
            logger.info(f"加载GAC数据: {gac_file}")
            df = pd.read_csv(gac_file, dtype=str)
            logger.info(f"  加载了 {len(df)} 条记录")
            
            l3_gac_map: Dict[str, str] = {}
            
            for _, row in tqdm(df.iterrows(), total=len(df), desc="提取GAC"):
                did = str(row.get('DID', '')).strip()
                gac = str(row.get('GAC', '')).strip() if pd.notna(row.get('GAC')) else "None"
                if did:
                    # Handle nan and empty values
                    if gac == '' or gac.lower() == 'nan' or gac.lower() == 'none':
                        gac = "None"
                    # Add L3_ prefix for consistency
                    l3_did = f"L3_{did}" if not did.startswith('L3_') else did
                    l3_gac_map[l3_did] = gac
            
            logger.info(f"L3 GAC映射提取完成: {len(l3_gac_map)} 个L3节点")
            
            # Save to file
            self._save_json(l3_gac_map, self.l3_gac_output)
            
            return l3_gac_map
        
        # Fallback to neo4j_repaired_ligands.csv
        repaired_ligands_file = self.input_dir / "neo4j_repaired_ligands.csv"
        if not repaired_ligands_file.exists():
            raise FileNotFoundError(f"GAC文件不存在: {gac_file} 和 {repaired_ligands_file}")
        
        logger.info(f"加载Repaired Ligands: {repaired_ligands_file}")
        df = pd.read_csv(repaired_ligands_file, dtype=str)
        logger.info(f"  加载了 {len(df)} 条记录")
        
        # Get column names
        did_col = [c for c in df.columns if 'did' in c.lower() or 'ID' in c][0]
        gac_col = [c for c in df.columns if 'gac' in c.lower()][0] if any('gac' in c.lower() for c in df.columns) else None
        
        l3_gac_map: Dict[str, str] = {}
        
        for _, row in tqdm(df.iterrows(), total=len(df), desc="提取GAC"):
            did = str(row[did_col]).strip()
            gac = str(row[gac_col]).strip() if gac_col and pd.notna(row[gac_col]) else "None"
            if did:
                # Handle nan and empty values
                if gac == '' or gac.lower() == 'nan' or gac.lower() == 'none':
                    gac = "None"
                l3_gac_map[did] = gac
        
        logger.info(f"L3 GAC映射提取完成: {len(l3_gac_map)} 个L3节点")
        
        # Save to file
        self._save_json(l3_gac_map, self.l3_gac_output)
        
        return l3_gac_map
    
    def build_filtered_l5_l3_index(self, l5_l3_index: Dict[str, List[str]], 
                                    l3_gac_map: Dict[str, str]) -> Dict[str, List[str]]:
        """
        构建按GAC过滤的L5→L3索引
        只保留GAC <= K的L3配体
        """
        logger.info(f"构建GAC过滤索引 (K={self.K})...")
        
        l5_l3_filtered: Dict[str, List[str]] = {}
        total_before = 0
        total_after = 0
        
        for l5_did, l3_list in tqdm(l5_l3_index.items(), desc=f"过滤索引 (K={self.K})"):
            total_before += len(l3_list)
            filtered_l3 = []
            
            for l3_did in l3_list:
                gac_str = l3_gac_map.get(l3_did, "None")
                if gac_str != "None":
                    try:
                        gac = int(float(gac_str))
                        if gac <= self.K:
                            filtered_l3.append(l3_did)
                    except (ValueError, TypeError):
                        pass
            
            if filtered_l3:
                l5_l3_filtered[l5_did] = sorted(filtered_l3)
                total_after += len(filtered_l3)
        
        self.stats['total_filtered_l5_l3_mappings'] = total_after
        
        if total_before > 0:
            logger.info(f"过滤前总L3数: {total_before:,}, 过滤后总L3数: {total_after:,} ({total_after/total_before*100:.1f}%)")
        logger.info(f"GAC过滤索引 (K={self.K}) 构建完成: {len(l5_l3_filtered)} 个L5单元")
        
        # Save to file
        self._save_json(l5_l3_filtered, self.l5_l3_filtered_output)
        
        return l5_l3_filtered
    
    def build_m_l3_pairs(self, l3_gac_map: Dict[str, str]) -> int:
        """
        构建M-L3对（金属-配体对）
        
        从M-L1关系和L1-L3关系推导M-L3关系，并筛选GAC范围
        """
        logger.info(f"构建M-L3对 (GAC范围: [{self.min_gac}, {self.max_gac}])...")
        
        # Load relationship files
        m_l1_file = self.input_dir / "neo4j_m_l1_relationships.csv"
        l1_l3_file = self.input_dir / "neo4j_l1_l3_relationships.csv"
        repaired_ligands_file = self.input_dir / "neo4j_repaired_ligands.csv"
        
        if not m_l1_file.exists():
            logger.warning(f"M-L1关系文件不存在: {m_l1_file}, 跳过M-L3对构建")
            return 0
        if not l1_l3_file.exists():
            logger.warning(f"L1-L3关系文件不存在: {l1_l3_file}, 跳过M-L3对构建")
            return 0
        
        # Load M-L1 relationships
        logger.info(f"加载M-L1关系: {m_l1_file}")
        m_l1_df = pd.read_csv(m_l1_file, dtype=str)
        logger.info(f"  加载了 {len(m_l1_df)} 条M-L1关系")
        
        # Load L1-L3 relationships
        logger.info(f"加载L1-L3关系: {l1_l3_file}")
        l1_l3_df = pd.read_csv(l1_l3_file, dtype=str)
        logger.info(f"  加载了 {len(l1_l3_df)} 条L1-L3关系")
        
        # Load repaired ligands for SMILES
        logger.info(f"加载Repaired Ligands: {repaired_ligands_file}")
        repaired_df = pd.read_csv(repaired_ligands_file, dtype=str)
        
        # Get column names for M-L1
        m_col = [c for c in m_l1_df.columns if 'START_ID' in c or 'metal' in c.lower()][0]
        l1_col_ml1 = [c for c in m_l1_df.columns if 'END_ID' in c or 'complex' in c.lower()][0]
        
        # Get column names for L1-L3
        l1_col_l1l3 = [c for c in l1_l3_df.columns if 'START_ID' in c or 'complex' in c.lower()][0]
        l3_col = [c for c in l1_l3_df.columns if 'END_ID' in c or 'repaired' in c.lower()][0]
        
        # Get column names for repaired ligands
        did_col = [c for c in repaired_df.columns if 'did' in c.lower() or 'ID' in c][0]
        smiles_col = [c for c in repaired_df.columns if 'smiles' in c.lower()][0] if any('smiles' in c.lower() for c in repaired_df.columns) else None
        
        # Build L3 DID to SMILES mapping
        l3_smiles_map: Dict[str, str] = {}
        if smiles_col:
            for _, row in repaired_df.iterrows():
                did = str(row[did_col]).strip()
                smiles = str(row[smiles_col]).strip() if pd.notna(row[smiles_col]) else ""
                if did and smiles:
                    l3_smiles_map[did] = smiles
        
        # Build L1→Metal mapping
        logger.info("构建L1→Metal映射...")
        l1_to_metals: Dict[str, Set[str]] = {}
        for _, row in tqdm(m_l1_df.iterrows(), total=len(m_l1_df), desc="L1→Metal"):
            metal_id = str(row[m_col]).strip()
            l1_id = str(row[l1_col_ml1]).strip()
            # Extract metal symbol from ID (e.g., "M_Cu" -> "Cu")
            metal_symbol = metal_id.replace("M_", "") if metal_id.startswith("M_") else metal_id
            if l1_id and metal_symbol:
                if l1_id not in l1_to_metals:
                    l1_to_metals[l1_id] = set()
                l1_to_metals[l1_id].add(metal_symbol)
        
        # Build M-L3 pairs
        logger.info("构建M-L3对...")
        m_l3_pairs: List[Dict] = []
        seen_pairs: Set[Tuple[str, str]] = set()
        
        for _, row in tqdm(l1_l3_df.iterrows(), total=len(l1_l3_df), desc="M-L3对"):
            l1_id = str(row[l1_col_l1l3]).strip()
            l3_id = str(row[l3_col]).strip()
            
            if not l1_id or not l3_id:
                continue
            
            # Check GAC range
            gac_str = l3_gac_map.get(l3_id, "None")
            if gac_str == "None":
                continue
            
            try:
                gac = int(float(gac_str))
            except (ValueError, TypeError):
                continue
            
            if gac < self.min_gac or gac > self.max_gac:
                continue
            
            # Get metals for this L1
            metals = l1_to_metals.get(l1_id, set())
            if not metals:
                continue
            
            # Get SMILES
            smiles = l3_smiles_map.get(l3_id, "")
            
            # Remove L3_ prefix for ligand_did
            ligand_did = l3_id.replace("L3_", "") if l3_id.startswith("L3_") else l3_id
            
            for metal in metals:
                pair_key = (ligand_did, metal)
                if pair_key not in seen_pairs:
                    seen_pairs.add(pair_key)
                    m_l3_pairs.append({
                        'ligand_did': ligand_did,
                        'ligand_smiles': smiles,
                        'ligand_gac': gac,
                        'metal_symbol': metal
                    })
        
        # Save to CSV
        if m_l3_pairs:
            pairs_df = pd.DataFrame(m_l3_pairs)
            pairs_df.to_csv(self.m_l3_pairs_output, index=False)
            logger.info(f"M-L3对已保存到: {self.m_l3_pairs_output} ({len(m_l3_pairs)} 条记录)")
        else:
            logger.warning("没有生成M-L3对")
        
        self.stats['total_m_l3_pairs'] = len(m_l3_pairs)
        
        return len(m_l3_pairs)
    
    def process_all(self) -> Dict:
        """执行所有预处理步骤"""
        start_time = time.time()
        logger.info("开始Step12预处理")
        logger.info("=" * 60)
        
        # Step 1: Build L3→L5 mapping
        l3_l5_map = self.build_l3_l5_mapping()
        
        # Step 2: Build L5→L3 inverted index
        l5_l3_index = self.build_l5_l3_index(l3_l5_map)
        
        # Step 3: Compute L5 statistics
        l5_stats = self.compute_l5_statistics(l5_l3_index)
        
        # Step 4: Extract L3 GAC values
        l3_gac_map = self.extract_l3_gac()
        
        # Step 5: Build GAC-filtered L5→L3 index
        l5_l3_filtered = self.build_filtered_l5_l3_index(l5_l3_index, l3_gac_map)
        
        # Step 6: Build M-L3 pairs
        m_l3_count = self.build_m_l3_pairs(l3_gac_map)
        
        # Calculate total time
        total_time = time.time() - start_time
        self.stats['processing_time'] = total_time
        
        # Save statistics
        self._save_stats()
        
        logger.info("=" * 60)
        logger.info("Step12预处理完成！统计信息:")
        logger.info(f"  L3节点数: {self.stats['total_l3_nodes']}")
        logger.info(f"  L5节点数: {self.stats['total_l5_nodes']}")
        logger.info(f"  L3→L5映射数: {self.stats['total_l3_l5_mappings']}")
        logger.info(f"  L5→L3映射数: {self.stats['total_l5_l3_mappings']}")
        logger.info(f"  过滤后L5→L3映射数 (K={self.K}): {self.stats['total_filtered_l5_l3_mappings']}")
        logger.info(f"  M-L3对数: {self.stats['total_m_l3_pairs']}")
        logger.info(f"  处理时间: {total_time:.2f}秒")
        logger.info("=" * 60)
        
        return self.stats


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Step12: 从Step8输出构建预计算索引")
    parser.add_argument('--record-limit', type=int, default=DEFAULT_RECORD_LIMIT,
                       help=f'记录限制 (-1 = 无限制，默认: {DEFAULT_RECORD_LIMIT})')
    parser.add_argument('--K', type=int, default=DEFAULT_K,
                       help=f'GAC截断值 (默认: {DEFAULT_K})')
    parser.add_argument('--min-gac', type=int, default=DEFAULT_MIN_GAC,
                       help=f'M-L3对最小GAC (默认: {DEFAULT_MIN_GAC})')
    parser.add_argument('--max-gac', type=int, default=DEFAULT_MAX_GAC,
                       help=f'M-L3对最大GAC (默认: {DEFAULT_MAX_GAC})')
    parser.add_argument('--input-dir', type=str, default=INPUT_DIR,
                       help=f'输入目录 (默认: {INPUT_DIR})')
    parser.add_argument('--output-dir', type=str, default=OUTPUT_DIR,
                       help=f'输出目录 (默认: {OUTPUT_DIR})')
    
    args = parser.parse_args()
    
    # Create processor and run with specified directories
    processor = Step12Processor(
        record_limit=args.record_limit,
        K=args.K,
        min_gac=args.min_gac,
        max_gac=args.max_gac
    )
    
    # Override directories if specified
    processor.input_dir = Path(args.input_dir)
    processor.output_dir = Path(args.output_dir)
    processor.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Update output file paths
    processor.l3_l5_output = processor.output_dir / "l3_l5.json"
    processor.l5_l3_output = processor.output_dir / "l5_l3.json"
    processor.l5_stats_output = processor.output_dir / "l5_freq_weight.json"
    processor.l3_gac_output = processor.output_dir / "l3_gac.json"
    processor.l5_l3_filtered_output = processor.output_dir / f"l5_l3_filtered_K{processor.K}.json"
    processor.m_l3_pairs_output = processor.output_dir / "m_l3_pairs.csv"
    processor.stats_output = processor.output_dir / "step12_stats.json"
    
    processor.process_all()


if __name__ == "__main__":
    main()
