#!/usr/bin/env python3
"""
ChemDB Step8 Processor - Neo4j CSV Export
基于前面步骤的输出，生成Neo4j兼容的CSV文件用于图数据库导入
使用单处理器处理，适合磁盘密集型操作
"""

# Default configuration variables
METAL_LIST_PATH = "./data/metal_list.txt"
DEFAULT_RECORD_LIMIT = -1  # -1 for unlimited
OUTPUT_DIR = "./tmp"
INPUT_DIR = "./tmp"

import os
import json
import time
import logging
import ast
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import pandas as pd
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def extract_metals_from_metal_info(metal_info_str: str) -> List[str]:
    """Extract metal symbols from metal_info JSON string"""
    try:
        if not metal_info_str or metal_info_str == 'nan':
            return []
        
        metal_info = ast.literal_eval(metal_info_str)
        metals = []
        
        for metal_dict in metal_info:
            if isinstance(metal_dict, dict):
                for metal_symbol in metal_dict.keys():
                    if metal_symbol and metal_symbol not in metals:
                        metals.append(metal_symbol)
        
        return metals
    except Exception as e:
        logger.warning(f"解析metal_info失败: {metal_info_str}, 错误: {e}")
        return []

class Step8Processor:
    def __init__(
        self,
        record_limit: int = DEFAULT_RECORD_LIMIT,
        input_dir: str = None,
        output_dir: str = None,
        metal_list_path: str = None,
    ):
        self.input_dir = Path(input_dir or INPUT_DIR)
        self.output_dir = Path(output_dir or OUTPUT_DIR)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.metal_list_path = metal_list_path or METAL_LIST_PATH
        
        self.record_limit = record_limit
        
        self.stats = {
            'total_metals': 0,
            'total_complexes': 0,
            'total_ligands': 0,
            'total_repaired_ligands': 0,
            'total_fragments': 0,
            'total_irl': 0,
            'total_m_l1_relationships': 0,
            'total_l1_l2_relationships': 0,
            'total_l2_l3_relationships': 0,
            'total_l1_l3_relationships': 0,
            'total_l3_l4_relationships': 0,
            'total_l4_l5_relationships': 0,
            'processing_time': 0,
            'record_limit': self.record_limit
        }
        
        # Output files for Neo4j import
        self.metals_output = self.output_dir / "neo4j_metals.csv"
        self.complexes_output = self.output_dir / "neo4j_complexes.csv"
        self.ligands_output = self.output_dir / "neo4j_ligands.csv"
        self.repaired_ligands_output = self.output_dir / "neo4j_repaired_ligands.csv"
        self.fragments_output = self.output_dir / "neo4j_fragments.csv"
        self.irl_output = self.output_dir / "neo4j_irl.csv"
        self.m_l1_relationships_output = self.output_dir / "neo4j_m_l1_relationships.csv"
        self.l1_l2_relationships_output = self.output_dir / "neo4j_l1_l2_relationships.csv"
        self.l2_l3_relationships_output = self.output_dir / "neo4j_l2_l3_relationships.csv"
        self.l1_l3_relationships_output = self.output_dir / "neo4j_l1_l3_relationships.csv"
        self.l3_l4_relationships_output = self.output_dir / "neo4j_l3_l4_relationships.csv"
        self.l4_l5_relationships_output = self.output_dir / "neo4j_l4_l5_relationships.csv"
        self.stats_output = self.output_dir / "step8_stats.json"
        
        logger.info(f"Step8处理器初始化完成，输出目录: {self.output_dir}")
        logger.info(f"记录限制: {self.record_limit if self.record_limit > 0 else '无限制'}")
    
    def _load_metal_list(self) -> List[str]:
        """Load metal list from file"""
        metal_file = Path(self.metal_list_path)
        if not metal_file.exists():
            raise FileNotFoundError(f"金属列表文件不存在: {metal_file}")
        
        with open(metal_file, 'r') as f:
            metals = [line.strip() for line in f if line.strip()]
        
        logger.info(f"加载了 {len(metals)} 个金属")
        return metals
    
    def _save_data_to_csv(self, data: List[Dict], output_file: Path):
        if not data:
            logger.info(f"没有数据需要保存到: {output_file}")
            return
            
        df = pd.DataFrame(data)
        df.to_csv(output_file, index=False, encoding='utf-8')
        logger.info(f"数据已保存到: {output_file} ({len(data)} 条记录)")
    
    def _save_stats(self):
        stats_data = {
            'timestamp': datetime.now().isoformat(),
            'processing_time_seconds': self.stats['processing_time'],
            'statistics': self.stats,
            'output_files': {
                'metals': str(self.metals_output),
                'complexes': str(self.complexes_output),
                'ligands': str(self.ligands_output),
                'repaired_ligands': str(self.repaired_ligands_output),
                'fragments': str(self.fragments_output),
                'irl': str(self.irl_output),
                'm_l1_relationships': str(self.m_l1_relationships_output),
                'l1_l2_relationships': str(self.l1_l2_relationships_output),
                'l2_l3_relationships': str(self.l2_l3_relationships_output),
                'l1_l3_relationships': str(self.l1_l3_relationships_output),
                'l3_l4_relationships': str(self.l3_l4_relationships_output),
                'l4_l5_relationships': str(self.l4_l5_relationships_output)
            }
        }
        
        with open(self.stats_output, 'w', encoding='utf-8') as f:
            json.dump(stats_data, f, indent=2, ensure_ascii=False)
        
        logger.info(f"统计信息已保存到: {self.stats_output}")
    
    def _deduplicate_nodes(self, records: List[Dict], id_field: str, unique_field: str = None) -> List[Dict]:
        """Deduplicate nodes based on unique field (usually SMILES)"""
        if not records:
            return records
        
        # If no unique field specified, use the ID field
        if unique_field is None:
            unique_field = id_field
        
        # Create a dictionary to keep only unique records
        unique_records = {}
        for record in records:
            unique_key = record.get(unique_field, '')
            if unique_key and unique_key not in unique_records:
                unique_records[unique_key] = record
        
        deduplicated = list(unique_records.values())
        logger.info(f"去重: {len(records)} -> {len(deduplicated)} 条记录 (基于 {unique_field})")
        return deduplicated
    
    def _deduplicate_relationships(self, relationships: List[Dict]) -> List[Dict]:
        """Deduplicate relationships based on start_id, end_id, and type"""
        if not relationships:
            return relationships
        
        # Create a set to track unique relationships
        unique_relationships = set()
        deduplicated = []
        
        for rel in relationships:
            # Create a unique key from start_id, end_id, and type
            start_id = rel.get(':START_ID', '') or list(rel.values())[0]  # Get first value if no :START_ID
            end_id = rel.get(':END_ID', '') or list(rel.values())[1]      # Get second value if no :END_ID
            rel_type = rel.get(':TYPE', '') or list(rel.values())[2]      # Get third value if no :TYPE
            
            unique_key = f"{start_id}|{end_id}|{rel_type}"
            
            if unique_key not in unique_relationships:
                unique_relationships.add(unique_key)
                deduplicated.append(rel)
        
        logger.info(f"关系去重: {len(relationships)} -> {len(deduplicated)} 条记录")
        return deduplicated

    def process_all_data(self) -> Dict:
        start_time = time.time()
        logger.info("开始Neo4j CSV数据生成处理")
        
        # Load metal list
        metals = self._load_metal_list()
        
        # Create metal records with M prefix
        metal_records = []
        for metal in metals:
            metal_record = {
                'id:ID(Metal)': f"M_{metal}",
                'symbol': metal,
                # 'type': 'Metal',
                ':LABEL': 'Metal'
            }
            metal_records.append(metal_record)
        
        # Deduplicate metals (should be unique by symbol)
        metal_records = self._deduplicate_nodes(metal_records, 'id:ID(Metal)')
        
        self.stats['total_metals'] = len(metal_records)
        self._save_data_to_csv(metal_records, self.metals_output)
        
        # Load all CSV files once into memory for efficient processing
        logger.info("加载CSV文件到内存...")
        complex_file = self.input_dir / "complex_data.csv"
        ligand_file = self.input_dir / "ligand_data.csv"
        repaired_file = self.input_dir / "repaired_ligand_data.csv"
        fragment_file = self.input_dir / "fragments.csv"
        irl_file = self.input_dir / "IRL_filtered.csv"
        
        # Read CSV files with progress bar
        csv_files = [
            (complex_file, "Complex data"),
            (ligand_file, "Ligand data"), 
            (repaired_file, "Repaired ligand data"),
            (fragment_file, "Fragment data"),
            (irl_file, "IRL data")
        ]
        
        complex_df, ligand_df, repaired_df, fragment_df, irl_df = pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        
        for file_path, desc in tqdm(csv_files, desc="Loading CSV files"):
            if file_path.exists():
                if desc == "Complex data":
                    complex_df = pd.read_csv(file_path)
                    record_count = len(complex_df)
                elif desc == "Ligand data":
                    ligand_df = pd.read_csv(file_path)
                    record_count = len(ligand_df)
                elif desc == "Repaired ligand data":
                    repaired_df = pd.read_csv(file_path)
                    record_count = len(repaired_df)
                elif desc == "Fragment data":
                    fragment_df = pd.read_csv(file_path)
                    record_count = len(fragment_df)
                elif desc == "IRL data":
                    irl_df = pd.read_csv(file_path)
                    record_count = len(irl_df)
            else:
                logger.info(f"{desc} file not found, using empty DataFrame")
        
        # Process complex data
        logger.info("处理复合物数据...")
        complex_records = []
        m_l1_relationships = []
        
        # Apply record limit for complex processing
        complex_df_limited = complex_df.head(self.record_limit) if self.record_limit > 0 else complex_df
        
        for row in tqdm(complex_df_limited.itertuples(), total=len(complex_df_limited), desc="处理复合物"):
            # Create complex record (L1) with L1 prefix
            complex_record = {
                'did:ID(Complex)': f"L1_{getattr(row, 'did', '')}",
                'cid': getattr(row, 'cid', '') if hasattr(row, 'cid') and getattr(row, 'cid') and getattr(row, 'cid') != 'nan' else '',
                'smiles': getattr(row, 'complex_smiles', ''),
                # 'metal_info': getattr(row, 'metal_info', ''),
                # 'ligand_count': getattr(row, 'ligand_count', 0),
                # 'from_pubchem': getattr(row, 'from_pubchem', True),
                'inactive': getattr(row, 'inactive', 0),
                ':LABEL': 'Complex'
            }
            complex_records.append(complex_record)
            
            # Extract metals and create M-L1 relationships
            metals = extract_metals_from_metal_info(getattr(row, 'metal_info', ''))
            for metal in metals:
                metal_rel = {
                    'metal_symbol:START_ID(Metal)': f"M_{metal}",
                    'complex_did:END_ID(Complex)': f"L1_{getattr(row, 'did', '')}",
                    ':TYPE': 'M_L1'
                }
                m_l1_relationships.append(metal_rel)
        
        # Deduplicate complexes based on DID
        complex_records = self._deduplicate_nodes(complex_records, 'did:ID(Complex)')
        
        self.stats['total_complexes'] = len(complex_records)
        # Deduplicate M-L1 relationships
        m_l1_relationships = self._deduplicate_relationships(m_l1_relationships)
        
        self.stats['total_m_l1_relationships'] = len(m_l1_relationships)
        self._save_data_to_csv(complex_records, self.complexes_output)
        self._save_data_to_csv(m_l1_relationships, self.m_l1_relationships_output)
        
        # Process ligand data
        logger.info("处理配体数据...")
        if not ligand_df.empty:
            ligand_records = []
            
            # Apply record limit for ligand processing
            ligand_df_limited = ligand_df.head(self.record_limit) if self.record_limit > 0 else ligand_df
            
            for row in tqdm(ligand_df_limited.itertuples(), total=len(ligand_df_limited), desc="处理配体"):
                # Create ligand record (L2) with L2 prefix, remove relationship fields
                ligand_record = {
                    'did:ID(Ligand)': f"L2_{getattr(row, 'ligand_did', '')}",
                    'smiles': getattr(row, 'ligand_smiles', ''),
                    ':LABEL': 'Ligand'
                }
                ligand_records.append(ligand_record)
            
            # Deduplicate ligands based on DID
            ligand_records = self._deduplicate_nodes(ligand_records, 'did:ID(Ligand)')
            
            self.stats['total_ligands'] = len(ligand_records)
            self._save_data_to_csv(ligand_records, self.ligands_output)
        
        # Process repaired ligand data
        logger.info("处理修复配体数据...")
        if not repaired_df.empty:
            repaired_records = []
            
            # Apply record limit for repaired ligand processing
            repaired_df_limited = repaired_df.head(self.record_limit) if self.record_limit > 0 else repaired_df
            
            for row in tqdm(repaired_df_limited.itertuples(), total=len(repaired_df_limited), desc="处理修复配体"):
                # Create repaired ligand record (L3) with L3 prefix, remove relationship fields
                repaired_record = {
                    'did:ID(RepairedLigand)': f"L3_{getattr(row, 'ligand_new_did', '')}",
                    'smiles': getattr(row, 'new_smiles', ''),
                    # 'old_smiles': getattr(row, 'old_smiles', ''),
                    'is_repaired': getattr(row, 'is_repaired', 0),
                    ':LABEL': 'RepairedLigand'
                }
                repaired_records.append(repaired_record)
            
            # Deduplicate repaired ligands based on DID
            repaired_records = self._deduplicate_nodes(repaired_records, 'did:ID(RepairedLigand)')
            
            self.stats['total_repaired_ligands'] = len(repaired_records)
            self._save_data_to_csv(repaired_records, self.repaired_ligands_output)
        
        # Process fragment data
        logger.info("处理片段数据...")
        if not fragment_df.empty:
            fragment_records = []
            
            # Apply record limit for fragment processing
            fragment_df_limited = fragment_df.head(self.record_limit) if self.record_limit > 0 else fragment_df
            
            for row in tqdm(fragment_df_limited.itertuples(), total=len(fragment_df_limited), desc="处理片段"):
                # Create fragment record (L4) with L4 prefix, remove relationship fields
                fragment_record = {
                    'did:ID(Fragment)': f"L4_{getattr(row, 'fragment_DID', '')}",
                    'smiles': getattr(row, 'fragment_smiles', ''),
                    # 'fragment_IRL_did': getattr(row, 'fragment_IRL_did', ''),
                    # 'fragment_IRL_smiles': getattr(row, 'fragment_IRL_smiles', ''),
                    ':LABEL': 'Fragment'
                }
                fragment_records.append(fragment_record)
            
            # Deduplicate fragments based on DID
            fragment_records = self._deduplicate_nodes(fragment_records, 'did:ID(Fragment)')
            
            self.stats['total_fragments'] = len(fragment_records)
            self._save_data_to_csv(fragment_records, self.fragments_output)
        
        # Process IRL data
        logger.info("处理IRL数据...")
        if not irl_df.empty:
            irl_records = []
            
            # Apply record limit for IRL processing
            irl_df_limited = irl_df.head(self.record_limit) if self.record_limit > 0 else irl_df
            
            for row in tqdm(irl_df_limited.itertuples(), total=len(irl_df_limited), desc="处理IRL"):
                # Create IRL record (L5) with L5 prefix
                irl_record = {
                    'did:ID(IRL)': f"L5_{getattr(row, 'DID', '')}",
                    'smiles': getattr(row, 'SMILES', ''),
                    'GAC': getattr(row, 'GAC', 0),
                    ':LABEL': 'IRL'
                }
                irl_records.append(irl_record)
            
            # Deduplicate IRL based on DID
            irl_records = self._deduplicate_nodes(irl_records, 'did:ID(IRL)')
            
            self.stats['total_irl'] = len(irl_records)
            self._save_data_to_csv(irl_records, self.irl_output)
        
        # Generate all relationships
        logger.info("生成关系数据...")
        
        # Create mappings for deduplicated nodes
        complex_id_mapping = {record['did:ID(Complex)'].replace('L1_', ''): record['did:ID(Complex)'] for record in complex_records}
        ligand_id_mapping = {record['did:ID(Ligand)'].replace('L2_', ''): record['did:ID(Ligand)'] for record in ligand_records}
        repaired_ligand_id_mapping = {record['did:ID(RepairedLigand)'].replace('L3_', ''): record['did:ID(RepairedLigand)'] for record in repaired_records}
        fragment_id_mapping = {record['did:ID(Fragment)'].replace('L4_', ''): record['did:ID(Fragment)'] for record in fragment_records}
        irl_id_mapping = {record['did:ID(IRL)'].replace('L5_', ''): record['did:ID(IRL)'] for record in irl_records}
        
        # Generate L1-L2 relationships (Complex to Ligand) - from ligand data
        logger.info("生成L1-L2关系...")
        l1_l2_relationships = []
        if not ligand_df.empty:
            # Use the same limited DataFrame for relationship generation
            ligand_df_for_rel = ligand_df.head(self.record_limit) if self.record_limit > 0 else ligand_df
            logger.info(f"处理 {len(ligand_df_for_rel)} 条配体记录...")
            for row in tqdm(ligand_df_for_rel.itertuples(), total=len(ligand_df_for_rel), desc="生成L1-L2关系"):
                complex_original_id = getattr(row, 'source_complex_did', '')
                ligand_original_id = getattr(row, 'ligand_did', '')
                
                # Only create relationship if both nodes exist after deduplication
                if complex_original_id in complex_id_mapping and ligand_original_id in ligand_id_mapping:
                    relationship = {
                        'complex_did:START_ID(Complex)': complex_id_mapping[complex_original_id],
                        'ligand_did:END_ID(Ligand)': ligand_id_mapping[ligand_original_id],
                        ':TYPE': 'L1_L2'
                    }
                    l1_l2_relationships.append(relationship)
        
        # Generate L2-L3 relationships (Ligand to Repaired Ligand) - from repaired ligand data
        logger.info("生成L2-L3关系...")
        l2_l3_relationships = []
        if not repaired_df.empty:
            # Use the same limited DataFrame for relationship generation
            repaired_df_for_rel = repaired_df.head(self.record_limit) if self.record_limit > 0 else repaired_df
            logger.info(f"处理 {len(repaired_df_for_rel)} 条修复配体记录...")
            for row in tqdm(repaired_df_for_rel.itertuples(), total=len(repaired_df_for_rel), desc="生成L2-L3关系"):
                ligand_original_id = getattr(row, 'source_did', '')  # source_did points to original ligand
                repaired_original_id = getattr(row, 'ligand_new_did', '')
                
                # Only create relationship if both nodes exist after deduplication
                if ligand_original_id in ligand_id_mapping and repaired_original_id in repaired_ligand_id_mapping:
                    relationship = {
                        'ligand_did:START_ID(Ligand)': ligand_id_mapping[ligand_original_id],
                        'repaired_ligand_did:END_ID(RepairedLigand)': repaired_ligand_id_mapping[repaired_original_id],
                        ':TYPE': 'L2_L3'
                    }
                    l2_l3_relationships.append(relationship)
        
        # Deduplicate L1-L2 and L2-L3 relationships
        l1_l2_relationships = self._deduplicate_relationships(l1_l2_relationships)
        l2_l3_relationships = self._deduplicate_relationships(l2_l3_relationships)
        
        # Save L1-L2 and L2-L3 relationships
        self.stats['total_l1_l2_relationships'] = len(l1_l2_relationships)
        self.stats['total_l2_l3_relationships'] = len(l2_l3_relationships)
        self._save_data_to_csv(l1_l2_relationships, self.l1_l2_relationships_output)
        self._save_data_to_csv(l2_l3_relationships, self.l2_l3_relationships_output)
        
        # Generate L1-L3 relationships by combining L1-L2 and L2-L3
        logger.info("生成L1-L3关系...")
        l1_l3_relationships = []
        
        # Create a mapping from L2 to L3 for efficient lookup
        l2_to_l3_mapping = {}
        for rel in l2_l3_relationships:
            l2_id = rel['ligand_did:START_ID(Ligand)']
            l3_id = rel['repaired_ligand_did:END_ID(RepairedLigand)']
            l2_to_l3_mapping[l2_id] = l3_id
        
        # Combine L1-L2 with L2-L3 to create L1-L3
        logger.info(f"处理 {len(l1_l2_relationships)} 条L1-L2关系...")
        for rel in tqdm(l1_l2_relationships, desc="生成L1-L3关系"):
            l1_id = rel['complex_did:START_ID(Complex)']
            l2_id = rel['ligand_did:END_ID(Ligand)']
            
            # If this L2 has a corresponding L3, create L1-L3 relationship
            if l2_id in l2_to_l3_mapping:
                l3_id = l2_to_l3_mapping[l2_id]
                relationship = {
                    'complex_did:START_ID(Complex)': l1_id,
                    'repaired_ligand_did:END_ID(RepairedLigand)': l3_id,
                    ':TYPE': 'L1_L3'
                }
                l1_l3_relationships.append(relationship)
        
        # Deduplicate L1-L3 relationships
        l1_l3_relationships = self._deduplicate_relationships(l1_l3_relationships)
        
        self.stats['total_l1_l3_relationships'] = len(l1_l3_relationships)
        self._save_data_to_csv(l1_l3_relationships, self.l1_l3_relationships_output)
        
        # Generate L3-L4 relationships (Repaired Ligand to Fragment) - from fragment data
        logger.info("生成L3-L4关系...")
        l3_l4_relationships = []
        if not fragment_df.empty:
            # Use the same limited DataFrame for relationship generation
            fragment_df_for_rel = fragment_df.head(self.record_limit) if self.record_limit > 0 else fragment_df
            logger.info(f"处理 {len(fragment_df_for_rel)} 条片段记录...")
            for row in tqdm(fragment_df_for_rel.itertuples(), total=len(fragment_df_for_rel), desc="生成L3-L4关系"):
                repaired_original_id = getattr(row, 'source_did', '')  # source_did points to repaired ligand
                fragment_original_id = getattr(row, 'fragment_DID', '')
                
                # Only create relationship if both nodes exist after deduplication
                if repaired_original_id in repaired_ligand_id_mapping and fragment_original_id in fragment_id_mapping:
                    relationship = {
                        'repaired_ligand_did:START_ID(RepairedLigand)': repaired_ligand_id_mapping[repaired_original_id],
                        'fragment_did:END_ID(Fragment)': fragment_id_mapping[fragment_original_id],
                        ':TYPE': 'L3_L4'
                    }
                    l3_l4_relationships.append(relationship)
        
        # Deduplicate L3-L4 relationships
        l3_l4_relationships = self._deduplicate_relationships(l3_l4_relationships)
        
        self.stats['total_l3_l4_relationships'] = len(l3_l4_relationships)
        self._save_data_to_csv(l3_l4_relationships, self.l3_l4_relationships_output)
        
        # Generate L4-L5 relationships (Fragment to IRL) - from fragment data
        logger.info("生成L4-L5关系...")
        l4_l5_relationships = []
        if not fragment_df.empty:
            # Use the same limited DataFrame for relationship generation
            fragment_df_for_rel = fragment_df.head(self.record_limit) if self.record_limit > 0 else fragment_df
            logger.info(f"处理 {len(fragment_df_for_rel)} 条片段记录...")
            for row in tqdm(fragment_df_for_rel.itertuples(), total=len(fragment_df_for_rel), desc="生成L4-L5关系"):
                fragment_original_id = getattr(row, 'fragment_DID', '')
                irl_did = getattr(row, 'fragment_IRL_did', '')
                
                if irl_did and irl_did != 'nan':
                    # Only create relationship if both nodes exist after deduplication
                    if fragment_original_id in fragment_id_mapping and irl_did in irl_id_mapping:
                        relationship = {
                            'fragment_did:START_ID(Fragment)': fragment_id_mapping[fragment_original_id],
                            'irl_did:END_ID(IRL)': irl_id_mapping[irl_did],
                            ':TYPE': 'L4_L5'
                        }
                        l4_l5_relationships.append(relationship)
        
        # Deduplicate L4-L5 relationships
        l4_l5_relationships = self._deduplicate_relationships(l4_l5_relationships)
        
        self.stats['total_l4_l5_relationships'] = len(l4_l5_relationships)
        self._save_data_to_csv(l4_l5_relationships, self.l4_l5_relationships_output)
        
        total_time = time.time() - start_time
        self.stats['processing_time'] = total_time
        
        self._save_stats()
        
        logger.info("=" * 50)
        logger.info("Neo4j CSV数据生成完成！最终统计:")
        logger.info(f"记录限制: {self.record_limit if self.record_limit > 0 else '无限制'}")
        logger.info(f"金属数: {self.stats['total_metals']}")
        logger.info(f"复合物数: {self.stats['total_complexes']}")
        logger.info(f"配体数: {self.stats['total_ligands']}")
        logger.info(f"修复配体数: {self.stats['total_repaired_ligands']}")
        logger.info(f"片段数: {self.stats['total_fragments']}")
        logger.info(f"IRL数: {self.stats['total_irl']}")
        logger.info(f"M-L1关系数: {self.stats['total_m_l1_relationships']}")
        logger.info(f"L1-L2关系数: {self.stats['total_l1_l2_relationships']}")
        logger.info(f"L2-L3关系数: {self.stats['total_l2_l3_relationships']}")
        logger.info(f"L1-L3关系数: {self.stats['total_l1_l3_relationships']}")
        logger.info(f"L3-L4关系数: {self.stats['total_l3_l4_relationships']}")
        logger.info(f"L4-L5关系数: {self.stats['total_l4_l5_relationships']}")
        logger.info(f"总处理时间: {total_time:.2f}秒")
        logger.info(f"输出目录: {self.output_dir}")
        logger.info("=" * 50)
        
        return self.stats

def main():
    try:
        import argparse
        import sys

        if len(sys.argv) > 1 and not any(a.startswith("-") for a in sys.argv[1:]):
            record_limit = DEFAULT_RECORD_LIMIT
            try:
                record_limit = int(sys.argv[1])
                print(f"使用指定的记录限制: {record_limit}")
            except ValueError:
                print("无效的记录限制，使用默认值")
            processor = Step8Processor(record_limit=record_limit)
            stats = processor.process_all_data()
        else:
            parser = argparse.ArgumentParser(description="ChemDB Step8")
            parser.add_argument("--input-dir", type=str, default=INPUT_DIR)
            parser.add_argument("--output-dir", type=str, default=OUTPUT_DIR)
            parser.add_argument("--metal-list", type=str, default=METAL_LIST_PATH)
            parser.add_argument("--limit", type=int, default=DEFAULT_RECORD_LIMIT)
            args = parser.parse_args()
            processor = Step8Processor(
                record_limit=args.limit,
                input_dir=args.input_dir,
                output_dir=args.output_dir,
                metal_list_path=args.metal_list,
            )
            stats = processor.process_all_data()
        
        print(f"\n处理完成！结果保存在: {processor.output_dir}")
        print(f"金属数: {stats['total_metals']}")
        print(f"复合物数: {stats['total_complexes']}")
        print(f"配体数: {stats['total_ligands']}")
        print(f"修复配体数: {stats['total_repaired_ligands']}")
        print(f"片段数: {stats['total_fragments']}")
        print(f"IRL数: {stats['total_irl']}")
        print(f"M-L1关系数: {stats['total_m_l1_relationships']}")
        print(f"L1-L2关系数: {stats['total_l1_l2_relationships']}")
        print(f"L2-L3关系数: {stats['total_l2_l3_relationships']}")
        print(f"L1-L3关系数: {stats['total_l1_l3_relationships']}")
        print(f"L3-L4关系数: {stats['total_l3_l4_relationships']}")
        print(f"L4-L5关系数: {stats['total_l4_l5_relationships']}")
        print(f"处理速度: {sum([stats['total_complexes'], stats['total_ligands'], stats['total_repaired_ligands'], stats['total_fragments'], stats['total_irl']])/stats['processing_time']:.1f} 记录/秒")
        
    except Exception as e:
        logger.error(f"处理失败: {e}")
        import sys
        sys.exit(1)

if __name__ == "__main__":
    main()
