#!/usr/bin/env python3
"""
ChemDB Step4_5 Processor - GA和IRL提取Pipeline
基于step2的输出(repaired_ligand_data.csv)，提取GA模式并生成IRL列表

Pipeline步骤:
- Step 0: 对new_smiles进行去重
- Step 1: 提取GA模式（芳香环、三元环、四元环）
- Step 2: 为GA模式生成唯一ID
- Step 3: 计算每个分子的GAC
- Step 4: 提取IRL（基于GAC和子结构匹配）
- Step 5: 清洗IRL数据（去重）

输出文件:
- ligand_data_deduplicated.csv: 去重后的输入文件
- GA_with_id.csv: GA模式及其ID
- ligand_with_gac.csv: 包含GAC的分子数据
- IRL_filtered.csv: 提取的IRL列表
- IRL_filtered_cleaned.csv: 清洗后的IRL列表
"""

# Default configuration variables
OUTPUT_DIR = "./tmp"
INPUT_DIR = "./tmp"
DEFAULT_RECORD_LIMIT = -1  # -1 for unlimited
DEFAULT_WORKERS = None  # Will use CPU count if None

import os
import sys
import csv
import json
import time
import hashlib
import logging
import threading
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from multiprocessing import Pool, Process, Queue, Manager, cpu_count
from queue import Empty

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
RDLogger.DisableLog('rdMol.*')
RDLogger.DisableLog('rdDecomposition.*')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def get_optimal_worker_count() -> int:
    """Get optimal worker count based on CPU cores"""
    cpu_cores = cpu_count()
    optimal_workers = min(cpu_cores, 200)
    logger.info(f"CPU核心数: {cpu_cores}, 使用CPU核心数作为工作进程数")
    return optimal_workers


# ============================================================
# GA提取相关函数 (来自 FindGA.py)
# ============================================================

class FindGA:
    """GA模式查找器"""
    def __init__(self, smiles: str, object_DID: str):
        self.smiles = smiles
        self.DID = object_DID
        self.mol = Chem.MolFromSmiles(smiles)
        if not self.mol:
            raise ValueError(f"无效的 SMILES 字符串: {smiles}")
        self.marked_atoms = []

    def find_all_rings(self):
        """找到所有环结构"""
        return Chem.GetSymmSSSR(self.mol)

    def find_aromatic_rings(self):
        """查找所有芳香环"""
        rings = self.find_all_rings()
        aromatic_rings = [ring for ring in rings if all(self.mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)]
        return aromatic_rings

    def get_ring_from_indices(self, ring):
        """根据环的原子索引提取环的完整子分子"""
        new_mol = Chem.RWMol()
        atom_mapping = {}
        
        for idx in ring:
            atom = self.mol.GetAtomWithIdx(idx)
            new_idx = new_mol.AddAtom(atom)
            atom_mapping[idx] = new_idx
        
        for bond in self.mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if begin_idx in atom_mapping and end_idx in atom_mapping:
                new_mol.AddBond(atom_mapping[begin_idx], atom_mapping[end_idx], bond.GetBondType())
        
        new_mol.UpdatePropertyCache()
        return new_mol.GetMol()

    def get_ring_smiles(self, ring_size: int) -> List[str]:
        """查找所有指定大小的环并返回其 SMILES 结构式"""
        rings = self.find_all_rings()
        specified_size_rings = [ring for ring in rings if len(ring) == ring_size]

        ring_smiles = []
        for ring in specified_size_rings:
            submol = self.get_ring_from_indices(ring)
            smiles = Chem.MolToSmiles(submol, canonical=True)
            ring_smiles.append(smiles)

        return ring_smiles

    def get_three_and_four_membered_rings(self) -> Tuple[List[str], List[str]]:
        """同时查找三元环和四元环并返回其 SMILES 结构式"""
        three_membered_smiles = self.get_ring_smiles(3)
        four_membered_smiles = self.get_ring_smiles(4)
        return three_membered_smiles, four_membered_smiles

    def get_aromatic_ring_smiles(self) -> List[str]:
        """查找所有芳香环并返回其 SMILES 结构式"""
        rings = self.find_all_rings()
        aromatic_rings = [
            ring for ring in rings if all(self.mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        ]

        aromatic_ring_smiles = []
        for ring in aromatic_rings:
            submol = self.get_ring_from_indices(ring)
            if submol is None:
                continue
            smiles = Chem.MolToSmiles(submol, canonical=True)
            ring_mol = Chem.MolFromSmiles(smiles)
            if ring_mol is None:
                continue
            ring_smiles = Chem.MolToSmiles(ring_mol, canonical=True)
            if ring_smiles:
                aromatic_ring_smiles.append(ring_smiles)

        return list(set(aromatic_ring_smiles))


def process_single_molecule_for_GA(DID_smiles_pair: Tuple[str, str]) -> Tuple[Dict[str, Set[str]], List[str]]:
    """
    处理单个分子，提取GA模式
    返回: (GA_dict, non_kekulizable_smiles)
    """
    GA_dict = {}
    non_kekulizable_smiles = []
    
    mol_DID, mol_smiles = DID_smiles_pair
    
    try:
        ga_finder = FindGA(mol_smiles, mol_DID)
        
        # 获取芳香环的 SMILES
        aromatic_ring_smiles = ga_finder.get_aromatic_ring_smiles()
        
        # 处理芳香环
        for ring_smiles in aromatic_ring_smiles:
            ring_mol = Chem.MolFromSmiles(ring_smiles)
            if ring_mol is None:
                non_kekulizable_smiles.append(ring_smiles)
                continue
            
            heavy_atom_count = ring_mol.GetNumHeavyAtoms()
            if heavy_atom_count > 7:
                continue

            try:
                Chem.Kekulize(ring_mol, clearAromaticFlags=True)
                if ring_smiles not in GA_dict:
                    GA_dict[ring_smiles] = set()
                GA_dict[ring_smiles].add(mol_DID)
            except Exception:
                non_kekulizable_smiles.append(ring_smiles)

        # 获取三元环和四元环的 SMILES
        three_membered_smiles, four_membered_smiles = ga_finder.get_three_and_four_membered_rings()

        # 处理三元环
        for ring_smiles in three_membered_smiles:
            ring_mol = Chem.MolFromSmiles(ring_smiles)
            if ring_mol is None:
                non_kekulizable_smiles.append(ring_smiles)
                continue
            
            heavy_atom_count = ring_mol.GetNumHeavyAtoms()
            if heavy_atom_count > 7:
                continue

            if ring_smiles not in GA_dict:
                GA_dict[ring_smiles] = set()
            GA_dict[ring_smiles].add(mol_DID)

        # 处理四元环
        for ring_smiles in four_membered_smiles:
            ring_mol = Chem.MolFromSmiles(ring_smiles)
            if ring_mol is None:
                non_kekulizable_smiles.append(ring_smiles)
                continue
            
            heavy_atom_count = ring_mol.GetNumHeavyAtoms()
            if heavy_atom_count > 7:
                continue

            if ring_smiles not in GA_dict:
                GA_dict[ring_smiles] = set()
            GA_dict[ring_smiles].add(mol_DID)
    
    except Exception:
        pass
    
    return GA_dict, non_kekulizable_smiles


# ============================================================
# GAC计算相关函数 (来自 step4.py)
# ============================================================

def calculate_gac(smiles: str, ga_list: List[str]) -> int:
    """计算分子的GAC值"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0
        
        ga_count = 0
        matched_atoms = set()
        
        for ga_smiles in ga_list:
            ga_mol = Chem.MolFromSmiles(ga_smiles)
            if ga_mol is None:
                continue
            
            matches = mol.GetSubstructMatches(ga_mol)
            if matches:
                ga_count += len(matches)
                for match in matches:
                    matched_atoms.update(match)
        
        remaining_heavy_atoms = mol.GetNumHeavyAtoms()
        remaining_atoms = remaining_heavy_atoms - len(matched_atoms)
        
        gac_count = ga_count + remaining_atoms
        return gac_count
        
    except Exception as e:
        logger.error(f"计算GAC失败: {smiles}, 错误: {e}")
        return 0


def generate_ga_id(smiles: str) -> Optional[str]:
    """为GA模式生成唯一ID"""
    try:
        hash_object = hashlib.md5(smiles.encode())
        hash_hex = hash_object.hexdigest()
        first_six_hex = hash_hex[:12]
        hash_decimal = int(first_six_hex, 16)
        final_ga_id = f"G{hash_decimal % 1000000:06d}"
        return final_ga_id
    except Exception as e:
        logger.error(f"生成GA ID失败: {smiles}, 错误: {e}")
        return None


def _process_gac_worker(args):
    """计算单个分子的GAC（用于并行处理）"""
    did_smiles_pair, ga_list = args
    did, smiles = did_smiles_pair
    if not smiles:
        return None
    try:
        gac = calculate_gac(smiles, ga_list)
        return {'DID': did, 'SMILES': smiles, 'GAC': gac}
    except Exception:
        return None


# ============================================================
# IRL提取相关函数 (来自 extract_IRL_via_GAC.py)
# ============================================================

def _process_single_molecule_for_heavy_elements(args):
    """处理单个分子，计算非氢元素"""
    index, item = args
    smiles = item.get('complex_smiles', '')
    
    if not isinstance(smiles, str) or pd.isna(smiles):
        return (index, [])
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (index, [])
    
    heavy_elements = {atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetSymbol() != 'H'}
    return (index, sorted(heavy_elements))


# 全局变量用于在worker进程间共享数据（避免重复序列化）
_worker_molecule_list = None
_worker_mol_cache = None
_worker_shared_flags = None
_worker_shared_lock = None
_worker_shared_stats = None


def _init_worker(molecule_list, mol_cache, shared_flags, shared_lock, shared_stats):
    """初始化worker进程的全局变量"""
    global _worker_molecule_list, _worker_mol_cache, _worker_shared_flags, _worker_shared_lock, _worker_shared_stats
    _worker_molecule_list = molecule_list
    _worker_mol_cache = mol_cache
    _worker_shared_flags = shared_flags
    _worker_shared_lock = shared_lock
    _worker_shared_stats = shared_stats


def _process_single_molecule_for_IRL(args):
    """处理单个分子进行IRL提取（使用全局共享数据）"""
    i, current_gac = args
    
    global _worker_molecule_list, _worker_mol_cache, _worker_shared_flags, _worker_shared_lock, _worker_shared_stats
    
    molecule_list = _worker_molecule_list
    mol_cache = _worker_mol_cache
    shared_flags = _worker_shared_flags
    shared_lock = _worker_shared_lock
    shared_stats = _worker_shared_stats
    
    current_item = molecule_list[i]
    current_GAC = current_item['GAC']
    
    with shared_lock:
        if shared_flags.get(i, False):
            return None
    
    # 使用缓存的Mol对象
    current_mol = mol_cache.get(i)
    if current_mol is None:
        return None
    
    current_elements = current_item.get('heavy_elements', [])
    current_elements_set = set(current_elements)
    
    to_mark = []
    
    for j in range(i + 1, len(molecule_list)):
        target_item = molecule_list[j]
        
        if target_item['GAC'] <= current_GAC:
            continue
        
        # 使用缓存的Mol对象
        target_mol = mol_cache.get(j)
        if target_mol is None:
            continue
        
        target_elements = target_item.get('heavy_elements', [])
        
        # 元素检查（使用预计算的set）
        if not current_elements_set.issubset(set(target_elements)):
            continue
        
        # 子结构匹配
        if target_mol.HasSubstructMatch(current_mol):
            to_mark.append(j)
    
    marked_count = 0
    if to_mark:
        with shared_lock:
            to_mark_filtered = [j for j in to_mark if not shared_flags.get(j, False)]
            for j in to_mark_filtered:
                shared_flags[j] = True
            marked_count = len(to_mark_filtered)
    
    with shared_lock:
        shared_stats['processed'] = shared_stats.get('processed', 0) + 1
        shared_stats['marked'] = shared_stats.get('marked', 0) + marked_count
    
    return i


class ExtractIRLObject:
    """IRL提取器"""
    def __init__(self, molecule_list: List[Dict], num_workers: int = None):
        self.molecule_list = molecule_list
        self.num_workers = num_workers or cpu_count()
        self.compute_heavy_elements_parallel()

    def compute_heavy_elements_parallel(self):
        """并行计算每个分子包含的非氢元素"""
        args_list = [(i, item) for i, item in enumerate(self.molecule_list)]
        
        logger.info(f"使用 {self.num_workers} 个并行进程计算 {len(self.molecule_list)} 个分子的非氢元素...")
        with Pool(processes=self.num_workers) as pool:
            results = list(tqdm(
                pool.imap(_process_single_molecule_for_heavy_elements, args_list),
                total=len(args_list),
                desc="计算非氢元素"
            ))
        
        for index, heavy_elements in results:
            self.molecule_list[index]['heavy_elements'] = heavy_elements

    def filter_by_GAC(self):
        """筛选出 GAC >= 2 的分子，并排除SMILES为"CC"的分子"""
        original_count = len(self.molecule_list)
        self.molecule_list = [
            item for item in self.molecule_list 
            if item['GAC'] >= 2 and item.get('complex_smiles', '').strip() != 'CC'
        ]
        filtered_count = len(self.molecule_list)
        removed_count = original_count - filtered_count
        if removed_count > 0:
            logger.info(f"过滤: 移除了 {removed_count} 个分子 (GAC<2 或 SMILES='CC')")

    def extract_TIRL_parallel(self) -> List[Dict]:
        """并行版本的IRL提取（优化版：预缓存Mol对象，使用initializer共享数据）"""
        self.filter_by_GAC()
        
        logger.info(f"使用 {self.num_workers} 个并行进程处理 {len(self.molecule_list)} 个分子...")
        
        # 预先解析所有SMILES为Mol对象（避免重复解析）
        logger.info("预解析SMILES为Mol对象...")
        mol_cache = {}
        for i, item in enumerate(self.molecule_list):
            smiles = item.get('complex_smiles', '')
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    mol_cache[i] = mol
        logger.info(f"成功解析 {len(mol_cache)}/{len(self.molecule_list)} 个分子")
        
        # 按GAC分组
        gac_groups = {}
        for i, item in enumerate(self.molecule_list):
            gac = item['GAC']
            if gac not in gac_groups:
                gac_groups[gac] = []
            gac_groups[gac].append(i)
        
        sorted_gacs = sorted(gac_groups.keys())
        logger.info(f"分子按GAC分组: {len(sorted_gacs)} 个GAC组")
        
        manager = Manager()
        shared_flags = manager.dict()
        shared_lock = manager.Lock()
        shared_stats = manager.dict({'processed': 0, 'marked': 0})
        
        # 使用单个Pool处理所有GAC组（避免重复创建Pool的开销）
        # 但仍然按GAC顺序处理以保证正确性
        for gac in sorted_gacs:
            indices = gac_groups[gac]
            
            # 过滤有效索引（已在mol_cache中且未被标记）
            valid_indices = []
            for i in indices:
                if i not in mol_cache:
                    continue
                with shared_lock:
                    if shared_flags.get(i, False):
                        continue
                valid_indices.append(i)
            
            if not valid_indices:
                continue
            
            # 简化的参数列表（不再传递molecule_list）
            args_list = [(i, gac) for i in valid_indices]
            
            # 使用initializer传递共享数据（只序列化一次）
            with Pool(
                processes=self.num_workers,
                initializer=_init_worker,
                initargs=(self.molecule_list, mol_cache, shared_flags, shared_lock, shared_stats)
            ) as pool:
                list(tqdm(
                    pool.imap(_process_single_molecule_for_IRL, args_list, 
                             chunksize=max(1, len(args_list) // (self.num_workers * 4))),
                    total=len(args_list),
                    desc=f"GAC={gac}",
                    mininterval=5.0
                ))
        
        for i, item in enumerate(self.molecule_list):
            if i in shared_flags and shared_flags[i]:
                item['flag'] = True
        
        logger.info(f"处理完成: 处理了 {shared_stats['processed']} 个分子，标记了 {shared_stats['marked']} 个复杂分子")
        
        return self.molecule_list


# ============================================================
# 主处理器类
# ============================================================

class Step4_5Processor:
    """Step4_5处理器：GA和IRL提取Pipeline"""
    
    def __init__(self, num_workers: int = None, record_limit: int = DEFAULT_RECORD_LIMIT,
                 input_dir: str = None, output_dir: str = None):
        self.output_dir = Path(output_dir or OUTPUT_DIR)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.input_dir = Path(input_dir or INPUT_DIR)
        
        if num_workers is None:
            self.num_workers = get_optimal_worker_count()
        else:
            self.num_workers = min(num_workers, 200)
        
        self.record_limit = record_limit
        
        # 输出文件路径
        self.deduplicated_output = self.output_dir / "ligand_data_deduplicated.csv"
        self.ga_temp_output = self.output_dir / "ga_patterns_temp.csv"
        self.ga_output = self.output_dir / "GA_with_id.csv"
        self.gac_output = self.output_dir / "ligand_with_gac.csv"
        self.irl_output = self.output_dir / "IRL_filtered.csv"
        self.irl_cleaned_output = self.output_dir / "IRL_filtered_cleaned.csv"
        self.stats_output = self.output_dir / "step4_5_stats.json"
        
        self.stats = {
            'start_time': time.time(),
            'num_workers': self.num_workers,
            'record_limit': self.record_limit,
            'steps': {}
        }
        
        logger.info(f"Step4_5处理器初始化完成，输出目录: {self.output_dir}")
        logger.info(f"并行处理工作进程数: {self.num_workers}")
        logger.info(f"记录限制: {self.record_limit if self.record_limit > 0 else '无限制'}")

    def _get_smiles_DID_from_csv(self, csv_file_path: Path) -> List[Tuple[str, str]]:
        """从CSV文件读取DID和SMILES"""
        smiles_with_did = []
        with open(csv_file_path, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                did = row.get("ligand_new_did") or row.get("source_did", "")
                smiles = row.get("new_smiles") or row.get("old_smiles", "")
                if did and smiles:
                    smiles_with_did.append((did, smiles))
        return smiles_with_did

    def step0_deduplicate(self) -> bool:
        """步骤0: 对new_smiles进行去重"""
        logger.info("=" * 60)
        logger.info("[步骤 0] 去重new_smiles")
        logger.info("=" * 60)
        
        input_csv = self.input_dir / "repaired_ligand_data.csv"
        if not input_csv.exists():
            logger.error(f"输入文件不存在: {input_csv}")
            return False
        
        df = pd.read_csv(input_csv)
        original_count = len(df)
        logger.info(f"读取了 {original_count} 行数据")
        
        if 'new_smiles' not in df.columns:
            logger.warning("输入文件没有'new_smiles'列，跳过去重步骤")
            df.to_csv(self.deduplicated_output, index=False)
            self.stats['steps']['step0'] = {
                'original_count': original_count,
                'deduplicated_count': original_count,
                'removed_count': 0
            }
            return True
        
        df_deduplicated = df.drop_duplicates(subset=['new_smiles'], keep='first')
        deduplicated_count = len(df_deduplicated)
        removed_count = original_count - deduplicated_count
        
        df_deduplicated.to_csv(self.deduplicated_output, index=False)
        
        logger.info(f"去重后: {deduplicated_count} 行，移除了 {removed_count} 个重复记录")
        logger.info(f"去重率: {removed_count/original_count*100:.2f}%")
        
        self.stats['steps']['step0'] = {
            'original_count': original_count,
            'deduplicated_count': deduplicated_count,
            'removed_count': removed_count
        }
        
        return True

    def step1_extract_GA(self) -> bool:
        """步骤1: 提取GA模式"""
        logger.info("=" * 60)
        logger.info("[步骤 1] 提取GA模式")
        logger.info("=" * 60)
        
        input_file = self.deduplicated_output if self.deduplicated_output.exists() else (self.input_dir / "repaired_ligand_data.csv")
        
        smiles_and_did_list = self._get_smiles_DID_from_csv(input_file)
        logger.info(f"读取了 {len(smiles_and_did_list)} 个分子")
        
        if self.record_limit > 0:
            smiles_and_did_list = smiles_and_did_list[:self.record_limit]
            logger.info(f"应用记录限制: {self.record_limit}")
        
        # 并行提取GA模式
        logger.info(f"使用 {self.num_workers} 个并行进程处理...")
        
        merged_GA_dict = {}
        
        with Pool(processes=self.num_workers) as pool:
            results = list(tqdm(
                pool.imap(process_single_molecule_for_GA, smiles_and_did_list),
                total=len(smiles_and_did_list),
                desc="提取GA模式"
            ))
        
        for GA_dict, _ in results:
            for ring_smiles, dids in GA_dict.items():
                if ring_smiles not in merged_GA_dict:
                    merged_GA_dict[ring_smiles] = set()
                merged_GA_dict[ring_smiles].update(dids)
        
        # 保存GA模式
        GA_list = []
        for ring_smiles, dids in merged_GA_dict.items():
            GA_list.append({"GA_SMILES": ring_smiles, "GA_source": ', '.join(dids)})
        
        if GA_list:
            df = pd.DataFrame(GA_list)
            df.to_csv(self.ga_temp_output, index=False)
            logger.info(f"已保存 {len(GA_list)} 个GA模式到 {self.ga_temp_output}")
        
        self.stats['steps']['step1'] = {
            'molecules_processed': len(smiles_and_did_list),
            'ga_patterns_found': len(GA_list)
        }
        
        return True

    def step2_generate_GA_ID(self) -> bool:
        """步骤2: 为GA模式生成ID"""
        logger.info("=" * 60)
        logger.info("[步骤 2] 为GA模式生成ID")
        logger.info("=" * 60)
        
        if not self.ga_temp_output.exists():
            logger.error(f"步骤1的输出文件不存在: {self.ga_temp_output}")
            return False
        
        ga_df = pd.read_csv(self.ga_temp_output)
        logger.info(f"读取了 {len(ga_df)} 个GA模式")
        
        ga_with_id = []
        for _, row in tqdm(ga_df.iterrows(), total=len(ga_df), desc="生成GA_ID"):
            ga_smiles = row['GA_SMILES']
            ga_id = generate_ga_id(ga_smiles)
            if ga_id:
                ga_with_id.append({'GA_SMILES': ga_smiles, 'GA_ID': ga_id})
        
        ga_id_df = pd.DataFrame(ga_with_id)
        ga_id_df.to_csv(self.ga_output, index=False)
        logger.info(f"已保存 {len(ga_with_id)} 个GA模式到 {self.ga_output}")
        
        # 删除临时文件
        if self.ga_temp_output.exists():
            self.ga_temp_output.unlink()
            logger.info(f"已删除临时文件: {self.ga_temp_output.name}")
        
        self.stats['steps']['step2'] = {
            'ga_patterns_with_id': len(ga_with_id)
        }
        
        return True

    def step3_calculate_GAC(self) -> bool:
        """步骤3: 计算GAC"""
        logger.info("=" * 60)
        logger.info("[步骤 3] 计算GAC")
        logger.info("=" * 60)
        
        if not self.ga_output.exists():
            logger.error(f"步骤2的输出文件不存在: {self.ga_output}")
            return False
        
        ga_id_df = pd.read_csv(self.ga_output)
        ga_list = ga_id_df['GA_SMILES'].tolist()
        logger.info(f"读取了 {len(ga_list)} 个GA模式用于计算GAC")
        
        input_file = self.deduplicated_output if self.deduplicated_output.exists() else (self.input_dir / "repaired_ligand_data.csv")
        smiles_and_did_list = self._get_smiles_DID_from_csv(input_file)
        logger.info(f"读取了 {len(smiles_and_did_list)} 个分子")
        
        if self.record_limit > 0:
            smiles_and_did_list = smiles_and_did_list[:self.record_limit]
        
        args_list = [(pair, ga_list) for pair in smiles_and_did_list]
        
        logger.info(f"使用 {self.num_workers} 个并行进程计算GAC...")
        with Pool(processes=self.num_workers) as pool:
            results = list(tqdm(
                pool.imap(_process_gac_worker, args_list),
                total=len(args_list),
                desc="计算GAC"
            ))
        
        gac_results = [r for r in results if r is not None]
        
        gac_df = pd.DataFrame(gac_results)
        gac_df.to_csv(self.gac_output, index=False)
        logger.info(f"已保存 {len(gac_results)} 个分子的GAC到 {self.gac_output}")
        
        self.stats['steps']['step3'] = {
            'molecules_with_gac': len(gac_results)
        }
        
        return True

    def step4_extract_IRL(self) -> bool:
        """步骤4: 提取IRL"""
        logger.info("=" * 60)
        logger.info("[步骤 4] 提取IRL")
        logger.info("=" * 60)
        
        if not self.gac_output.exists():
            logger.error(f"步骤3的输出文件不存在: {self.gac_output}")
            return False
        
        gac_df = pd.read_csv(self.gac_output)
        logger.info(f"读取了 {len(gac_df)} 个分子的GAC数据")
        
        # 准备分子列表 - 使用to_dict比iterrows快很多
        gac_df = gac_df.rename(columns={'SMILES': 'complex_smiles'})
        molecule_list = gac_df.to_dict(orient='records')
        
        # 按GAC排序
        molecule_list.sort(key=lambda x: x['GAC'])
        
        # 提取IRL
        extractor = ExtractIRLObject(molecule_list, num_workers=self.num_workers)
        result = extractor.extract_TIRL_parallel()
        
        # 只保存flag=False的IRL
        irl_results = [item for item in result if not item.get('flag', False)]
        
        if irl_results:
            irl_df = pd.DataFrame(irl_results)
            irl_output_df = irl_df[['DID', 'complex_smiles', 'GAC']].copy()
            irl_output_df = irl_output_df.rename(columns={'complex_smiles': 'SMILES'})
            irl_output_df.to_csv(self.irl_output, index=False)
            logger.info(f"已保存 {len(irl_results)} 个IRL到 {self.irl_output}")
        else:
            logger.warning("没有找到IRL")
        
        self.stats['steps']['step4'] = {
            'irl_count': len(irl_results)
        }
        
        return True

    def step5_deduplicate_IRL(self) -> bool:
        """步骤5: 清洗IRL数据"""
        logger.info("=" * 60)
        logger.info("[步骤 5] 清洗IRL数据（去重）")
        logger.info("=" * 60)
        
        if not self.irl_output.exists():
            logger.error(f"步骤4的输出文件不存在: {self.irl_output}")
            return False
        
        irl_df = pd.read_csv(self.irl_output)
        original_count = len(irl_df)
        logger.info(f"读取了 {original_count} 行IRL数据")
        
        irl_df_cleaned = irl_df.drop_duplicates(subset=['DID', 'SMILES'], keep='first')
        cleaned_count = len(irl_df_cleaned)
        removed_count = original_count - cleaned_count
        
        logger.info(f"去重前: {original_count} 行")
        logger.info(f"去重后: {cleaned_count} 行")
        logger.info(f"移除重复: {removed_count} 行")
        
        irl_df_cleaned.to_csv(self.irl_cleaned_output, index=False)
        logger.info(f"已保存清洗后的 {cleaned_count} 行IRL数据到 {self.irl_cleaned_output}")
        
        self.stats['steps']['step5'] = {
            'original_count': original_count,
            'cleaned_count': cleaned_count,
            'removed_count': removed_count
        }
        
        return True

    def _save_stats(self):
        """保存统计信息"""
        self.stats['processing_time'] = time.time() - self.stats['start_time']
        self.stats['timestamp'] = datetime.now().isoformat()
        
        with open(self.stats_output, 'w', encoding='utf-8') as f:
            json.dump(self.stats, f, indent=2, ensure_ascii=False)
        
        logger.info(f"统计信息已保存到: {self.stats_output}")

    def process_all(self) -> Dict:
        """运行完整Pipeline"""
        logger.info("=" * 60)
        logger.info("GA和IRL提取Pipeline")
        logger.info("=" * 60)
        
        success = True
        
        # Step 0: 去重
        if success:
            success = self.step0_deduplicate()
        
        # Step 1: 提取GA模式
        if success:
            success = self.step1_extract_GA()
        
        # Step 2: 生成GA ID
        if success:
            success = self.step2_generate_GA_ID()
        
        # Step 3: 计算GAC
        if success:
            success = self.step3_calculate_GAC()
        
        # Step 4: 提取IRL
        if success:
            success = self.step4_extract_IRL()
        
        # Step 5: 清洗IRL
        if success:
            success = self.step5_deduplicate_IRL()
        
        self._save_stats()
        
        # 总结
        logger.info("=" * 60)
        if success:
            logger.info("Pipeline执行完成！")
            logger.info("=" * 60)
            
            if self.deduplicated_output.exists():
                dedup_df = pd.read_csv(self.deduplicated_output)
                logger.info(f"✓ 去重后的输入文件: {self.deduplicated_output} ({len(dedup_df)} 行)")
            
            if self.ga_output.exists():
                ga_id_df = pd.read_csv(self.ga_output)
                logger.info(f"✓ GA列表(含ID): {self.ga_output} ({len(ga_id_df)} 个)")
            
            if self.gac_output.exists():
                gac_df = pd.read_csv(self.gac_output)
                logger.info(f"✓ GAC数据: {self.gac_output} ({len(gac_df)} 个分子)")
            
            if self.irl_output.exists():
                irl_df = pd.read_csv(self.irl_output)
                logger.info(f"✓ IRL列表: {self.irl_output} ({len(irl_df)} 个)")
            
            if self.irl_cleaned_output.exists():
                irl_cleaned_df = pd.read_csv(self.irl_cleaned_output)
                logger.info(f"✓ IRL列表(已清洗): {self.irl_cleaned_output} ({len(irl_cleaned_df)} 个)")
        else:
            logger.error("Pipeline执行失败！")
        
        logger.info("=" * 60)
        
        return self.stats

    def process_generate_ga(self) -> Dict:
        """Generate GA_with_id.csv only (no GAC/IRL)."""
        logger.info("=" * 60)
        logger.info("Mode: generate-ga (GA extraction + GA_ID only)")
        logger.info("=" * 60)
        success = True
        if success:
            success = self.step0_deduplicate()
        if success:
            success = self.step1_extract_GA()
        if success:
            success = self.step2_generate_GA_ID()
        self._save_stats()
        self.stats["mode"] = "generate-ga"
        self.stats["success"] = success
        if not success:
            logger.error("generate-ga failed")
        elif not self.ga_output.is_file():
            logger.error(f"generate-ga did not produce {self.ga_output}")
            success = False
            self.stats["success"] = False
        return self.stats

    def process_apply_ga(self) -> Dict:
        """Apply existing GA_with_id.csv: GAC + IRL only (do not overwrite GA list)."""
        logger.info("=" * 60)
        logger.info("Mode: apply-ga (GAC + IRL using existing GA_with_id.csv)")
        logger.info("=" * 60)
        if not self.ga_output.is_file():
            msg = (
                f"apply-ga requires existing GA_with_id.csv at {self.ga_output}; "
                "bind a GA version to the run first (ga_registry_manager bind-ga-version-to-run)."
            )
            logger.error(msg)
            self.stats["mode"] = "apply-ga"
            self.stats["success"] = False
            self.stats["error"] = msg
            self._save_stats()
            return self.stats
        success = True
        if success:
            success = self.step0_deduplicate()
        if success:
            success = self.step3_calculate_GAC()
        if success:
            success = self.step4_extract_IRL()
        if success:
            success = self.step5_deduplicate_IRL()
        self._save_stats()
        self.stats["mode"] = "apply-ga"
        self.stats["success"] = success
        if not success:
            logger.error("apply-ga failed")
        return self.stats


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='GA和IRL提取Pipeline')
    parser.add_argument(
        '--mode',
        type=str,
        default='full-auto',
        choices=('full-auto', 'generate-ga', 'apply-ga'),
        help='full-auto: legacy all steps; generate-ga: GA_with_id only; apply-ga: GAC/IRL from existing GA_with_id.csv',
    )
    parser.add_argument('--workers', '-w', type=int, default=None,
                        help='并行工作进程数 (默认: CPU核心数)')
    parser.add_argument('--limit', '-l', type=int, default=DEFAULT_RECORD_LIMIT,
                        help='记录限制 (默认: -1 无限制)')
    parser.add_argument('--input-dir', '-i', type=str, default=INPUT_DIR,
                        help=f'输入目录 (默认: {INPUT_DIR})')
    parser.add_argument('--output-dir', '-o', type=str, default=OUTPUT_DIR,
                        help=f'输出目录 (默认: {OUTPUT_DIR})')
    
    args = parser.parse_args()
    
    try:
        print(f"mode: {args.mode}")
        print(f"输入目录: {args.input_dir}")
        print(f"输出目录: {args.output_dir}")
        if args.workers:
            print(f"工作进程数: {args.workers}")
        if args.limit > 0:
            print(f"记录限制: {args.limit}")
        
        processor = Step4_5Processor(
            num_workers=args.workers, 
            record_limit=args.limit,
            input_dir=args.input_dir,
            output_dir=args.output_dir
        )
        if args.mode == 'full-auto':
            stats = processor.process_all()
        elif args.mode == 'generate-ga':
            stats = processor.process_generate_ga()
        else:
            stats = processor.process_apply_ga()
        
        print(f"\n处理完成！")
        print(f"处理时间: {stats['processing_time']:.2f} 秒")
        if stats.get("success") is False:
            return 1
        
        return 0
        
    except KeyboardInterrupt:
        logger.info("用户中断处理")
        print("\n处理被用户中断")
        return 1
        
    except Exception as e:
        logger.error(f"处理过程中发生错误: {e}")
        import traceback
        logger.error(f"错误详情: {traceback.format_exc()}")
        return 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
