#!/usr/bin/env python3
"""
ChemDB Step11 - L3 Embeddings and M-L3 Relationships
Outputs two CSV files:
1. L3 embeddings: did, smiles, embedding (from repaired_ligand_data.csv)
2. M-L3 relationships: symbol, did, smiles (constructed from M-L1 and L1-L3 relationships)

Uses single-process GPU batched processing for optimal performance.
"""

# Default configuration variables
INPUT_DIR = "./tmp"
OUTPUT_DIR = "./tmp"
DEFAULT_RECORD_LIMIT = -1  # -1 for unlimited
DEFAULT_BATCH_SIZE = 256  # Batch size for GPU processing

import json
import logging
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
import torch
from tqdm import tqdm

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# Import embedding functions from molclr_api (direct call, not via API)
from molclr_api import get_model, extract_embeddings, numpy_to_base64

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def get_most_free_gpu() -> str:
    """Select the GPU with most free memory using nvidia-smi."""
    if not torch.cuda.is_available():
        logger.info("CUDA不可用，使用CPU")
        return 'cpu'
    
    num_gpus = torch.cuda.device_count()
    if num_gpus == 0:
        return 'cpu'
    
    # Use nvidia-smi to get actual free memory
    import subprocess
    try:
        result = subprocess.run(
            ['nvidia-smi', '--query-gpu=index,name,memory.free,memory.total', '--format=csv,noheader,nounits'],
            capture_output=True, text=True, check=True
        )
        
        max_free_memory = 0
        best_gpu = 0
        
        for line in result.stdout.strip().split('\n'):
            parts = [p.strip() for p in line.split(',')]
            gpu_idx = int(parts[0])
            gpu_name = parts[1]
            free_mb = int(parts[2])
            total_mb = int(parts[3])
            
            logger.info(f"GPU {gpu_idx} ({gpu_name}): 空闲 {free_mb} MB / 总共 {total_mb} MB")
            
            if free_mb > max_free_memory:
                max_free_memory = free_mb
                best_gpu = gpu_idx
        
        device = f'cuda:{best_gpu}'
        logger.info(f"选择GPU {best_gpu}，空闲内存最大: {max_free_memory} MB")
        return device
        
    except Exception as e:
        logger.warning(f"nvidia-smi查询失败: {e}，使用默认GPU 0")
        return 'cuda:0'


class Step11Processor:
    def __init__(self, batch_size: int = DEFAULT_BATCH_SIZE, record_limit: int = DEFAULT_RECORD_LIMIT, 
                 model_type: str = "gin", device: str = None):
        self.input_dir = Path(INPUT_DIR)
        self.output_dir = Path(OUTPUT_DIR)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.batch_size = batch_size
        self.record_limit = record_limit
        self.model_type = model_type
        
        # Auto-select best GPU if not specified
        if device is None:
            self.device = get_most_free_gpu()
        else:
            self.device = device
        
        # Input files
        self.repaired_ligand_file = self.input_dir / "repaired_ligand_data.csv"
        self.m_l1_file = self.input_dir / "neo4j_m_l1_relationships.csv"
        self.l1_l3_file = self.input_dir / "neo4j_l1_l3_relationships.csv"
        self.metals_file = self.input_dir / "neo4j_metals.csv"
        
        # Output files
        self.output_embeddings_file = self.output_dir / "l3_embeddings.csv"
        self.output_m_l3_file = self.output_dir / "m_l3_relationships.csv"
        self.stats_output = self.output_dir / "step11_stats.json"
        
        # Statistics
        self.stats = {
            'total_ligands': 0,
            'successful_embeddings': 0,
            'failed_embeddings': 0,
            'm_l3_relationships': 0,
            'processing_time': 0,
            'batch_size': self.batch_size,
            'device': self.device,
            'model_type': self.model_type,
            'record_limit': self.record_limit
        }
        
        logger.info(f"Step11处理器初始化完成，输出目录: {self.output_dir}")
        logger.info(f"设备: {self.device}")
        logger.info(f"批处理大小: {self.batch_size}")
        logger.info(f"模型类型: {self.model_type}")
        logger.info(f"记录限制: {self.record_limit if self.record_limit > 0 else '无限制'}")
    
    def _load_repaired_ligands(self) -> List[Tuple[str, str]]:
        """Load repaired ligand data as list of (did, smiles) tuples."""
        logger.info(f"加载配体数据: {self.repaired_ligand_file}")
        df = pd.read_csv(self.repaired_ligand_file, dtype=str)
        
        # Keep only unique ligand_new_did with first occurrence
        df = df.drop_duplicates(subset=['ligand_new_did'], keep='first')
        
        # Select columns
        df = df[['ligand_new_did', 'new_smiles']].copy()
        df.columns = ['did', 'smiles']
        
        # Remove rows with missing values
        df = df.dropna(subset=['did', 'smiles'])
        df = df[df['smiles'].str.strip() != '']
        
        if self.record_limit > 0:
            df = df.head(self.record_limit)
        
        # Convert to list of tuples
        ligand_data = list(zip(df['did'].tolist(), df['smiles'].tolist()))
        
        logger.info(f"加载了 {len(ligand_data)} 个唯一配体")
        return ligand_data
    
    def _save_data_to_csv(self, data: List[Dict], output_file: Path):
        """Save data to CSV file"""
        if not data:
            logger.info(f"没有数据需要保存到: {output_file}")
            return
            
        df = pd.DataFrame(data)
        df.to_csv(output_file, index=False, encoding='utf-8')
        logger.info(f"数据已保存到: {output_file} ({len(data)} 条记录)")
    
    def _save_stats(self):
        """Save processing statistics"""
        stats_data = {
            'timestamp': datetime.now().isoformat(),
            'processing_time_seconds': self.stats['processing_time'],
            'statistics': self.stats,
            'output_files': {
                'embeddings': str(self.output_embeddings_file),
                'm_l3_relationships': str(self.output_m_l3_file)
            }
        }
        
        with open(self.stats_output, 'w', encoding='utf-8') as f:
            json.dump(stats_data, f, indent=2, ensure_ascii=False)
        
        logger.info(f"统计信息已保存到: {self.stats_output}")
    
    def _build_m_l3_relationships(self, ligand_data: List[Tuple[str, str]]) -> List[Dict]:
        """Build M-L3 relationships by joining M-L1 and L1-L3."""
        logger.info("构建M-L3关系")
        
        # Create ligand lookup for smiles
        ligand_smiles_map = {did: smiles for did, smiles in ligand_data}
        
        # Load M-L1 relationships
        logger.info(f"加载M-L1关系: {self.m_l1_file}")
        m_l1_df = pd.read_csv(self.m_l1_file, dtype=str)
        # Columns: metal_symbol:START_ID(Metal), complex_did:END_ID(Complex), :TYPE
        m_l1_df.columns = ['metal_id', 'complex_id', 'type']
        
        # Load L1-L3 relationships
        logger.info(f"加载L1-L3关系: {self.l1_l3_file}")
        l1_l3_df = pd.read_csv(self.l1_l3_file, dtype=str)
        # Columns: complex_did:START_ID(Complex), repaired_ligand_did:END_ID(RepairedLigand), :TYPE
        l1_l3_df.columns = ['complex_id', 'ligand_id', 'type']
        
        # Load metals for symbol lookup
        logger.info(f"加载金属列表: {self.metals_file}")
        metals_df = pd.read_csv(self.metals_file, dtype=str)
        # Columns: id:ID(Metal), symbol, :LABEL
        metals_df.columns = ['metal_id', 'symbol', 'label']
        metal_symbol_map = dict(zip(metals_df['metal_id'], metals_df['symbol']))
        
        # Join M-L1 with L1-L3 on complex_id
        logger.info("合并M-L1和L1-L3关系")
        m_l3_df = m_l1_df.merge(l1_l3_df, on='complex_id', how='inner', suffixes=('_m', '_l'))
        
        # Extract metal symbol from metal_id (M_Ca -> Ca)
        m_l3_df['symbol'] = m_l3_df['metal_id'].map(metal_symbol_map)
        
        # Extract did from ligand_id (L3_D177227588194109 -> D177227588194109)
        m_l3_df['did'] = m_l3_df['ligand_id'].str.replace('L3_', '', regex=False)
        
        # Map smiles from ligand data
        m_l3_df['smiles'] = m_l3_df['did'].map(ligand_smiles_map)
        
        # Select final columns and drop rows with missing smiles
        result_df = m_l3_df[['symbol', 'did', 'smiles']].copy()
        result_df = result_df.dropna(subset=['symbol', 'did', 'smiles'])
        
        # Remove duplicates (same metal-ligand pair)
        result_df = result_df.drop_duplicates(subset=['symbol', 'did'], keep='first')
        
        # Convert to list of dicts
        result_data = result_df.to_dict('records')
        
        logger.info(f"构建了 {len(result_data)} 个M-L3关系")
        return result_data
    
    def process_all(self) -> Dict:
        """Main processing function - uses GPU batched processing"""
        start_time = time.time()
        logger.info("开始L3嵌入向量和M-L3关系处理")
        
        # Step 1: Load repaired ligands
        ligand_data = self._load_repaired_ligands()
        
        if not ligand_data:
            logger.warning("没有找到有效的配体数据")
            return self.stats
        
        self.stats['total_ligands'] = len(ligand_data)
        
        # Step 2: Load model on selected device
        logger.info(f"加载模型到设备: {self.device}")
        model = get_model(self.model_type, device=self.device)
        logger.info("模型加载完成")
        
        # Step 3: Process in batches using GPU
        all_embedding_data = []
        total_failed = 0
        
        # Prepare data
        dids = [did for did, _ in ligand_data]
        smiles_list = [smiles for _, smiles in ligand_data]
        
        num_batches = (len(smiles_list) + self.batch_size - 1) // self.batch_size
        
        with tqdm(total=len(ligand_data), desc=f"GPU批处理嵌入向量 ({self.device})") as pbar:
            for batch_idx in range(num_batches):
                batch_start = batch_idx * self.batch_size
                batch_end = min(batch_start + self.batch_size, len(smiles_list))
                
                batch_dids = dids[batch_start:batch_end]
                batch_smiles = smiles_list[batch_start:batch_end]
                
                try:
                    # Process batch
                    embeddings, valid_indices, invalid_smiles = extract_embeddings(
                        model, batch_smiles, device=self.device
                    )
                    
                    # Map embeddings back to DIDs
                    for i, valid_idx in enumerate(valid_indices):
                        did = batch_dids[valid_idx]
                        smiles = batch_smiles[valid_idx]
                        embedding_base64 = numpy_to_base64(embeddings[i])
                        
                        all_embedding_data.append({
                            'did': did,
                            'smiles': smiles,
                            'embedding': embedding_base64
                        })
                    
                    total_failed += len(invalid_smiles)
                    
                except Exception as e:
                    logger.error(f"批次 {batch_idx} 处理失败: {e}")
                    total_failed += len(batch_smiles)
                
                pbar.update(len(batch_smiles))
                pbar.set_postfix({
                    'success': len(all_embedding_data),
                    'failed': total_failed
                })
        
        # Update stats
        self.stats['successful_embeddings'] = len(all_embedding_data)
        self.stats['failed_embeddings'] = total_failed
        
        # Step 4: Build M-L3 relationships
        m_l3_data = self._build_m_l3_relationships(ligand_data)
        self.stats['m_l3_relationships'] = len(m_l3_data)
        
        # Step 5: Save results
        logger.info("开始保存处理结果")
        
        # Sort results for consistent output order
        if all_embedding_data:
            all_embedding_data.sort(key=lambda x: x.get('did', ''))
            logger.info("嵌入向量数据已按DID排序")
        
        if m_l3_data:
            m_l3_data.sort(key=lambda x: (x.get('symbol', ''), x.get('did', '')))
            logger.info("M-L3关系数据已按symbol和DID排序")
        
        self._save_data_to_csv(all_embedding_data, self.output_embeddings_file)
        self._save_data_to_csv(m_l3_data, self.output_m_l3_file)
        
        total_time = time.time() - start_time
        self.stats['processing_time'] = total_time
        
        self._save_stats()
        
        logger.info("=" * 50)
        logger.info("L3嵌入向量和M-L3关系处理完成！最终统计:")
        logger.info(f"设备: {self.device}")
        logger.info(f"批处理大小: {self.batch_size}")
        logger.info(f"模型类型: {self.stats['model_type']}")
        logger.info(f"记录限制: {self.record_limit if self.record_limit > 0 else '无限制'}")
        logger.info(f"总配体数: {self.stats['total_ligands']}")
        logger.info(f"成功嵌入数: {self.stats['successful_embeddings']}")
        logger.info(f"失败数: {self.stats['failed_embeddings']}")
        logger.info(f"M-L3关系数: {self.stats['m_l3_relationships']}")
        logger.info(f"总处理时间: {total_time:.2f}秒")
        logger.info(f"平均处理速度: {self.stats['total_ligands']/total_time:.1f} 配体/秒")
        logger.info(f"输出目录: {self.output_dir}")
        logger.info("=" * 50)
        
        return self.stats


def main():
    """Main function
    
    Usage:
        python step11.py [batch_size] [record_limit] [model_type] [device]
        
    Examples:
        python step11.py                    # auto GPU, batch=128, unlimited
        python step11.py 256                # batch=256
        python step11.py 128 1000           # batch=128, limit=1000
        python step11.py 128 -1 gin cuda:0  # specific GPU
    """
    import sys
    
    try:
        batch_size = DEFAULT_BATCH_SIZE
        record_limit = DEFAULT_RECORD_LIMIT
        model_type = "gin"
        device = None  # Auto-select
        
        if len(sys.argv) > 1:
            try:
                batch_size = int(sys.argv[1])
                print(f"批处理大小: {batch_size}")
            except ValueError:
                print("无效的批处理大小，使用默认值")
        
        if len(sys.argv) > 2:
            try:
                record_limit = int(sys.argv[2])
                print(f"记录限制: {record_limit if record_limit > 0 else '无限制'}")
            except ValueError:
                print("无效的记录限制，使用默认值")
        
        if len(sys.argv) > 3:
            model_type = sys.argv[3]
            if model_type not in ['gin', 'gcn']:
                print(f"无效的模型类型 '{model_type}'，使用 'gin'")
                model_type = 'gin'
            else:
                print(f"模型类型: {model_type}")
        
        if len(sys.argv) > 4:
            device = sys.argv[4]
            print(f"指定设备: {device}")
        
        processor = Step11Processor(
            batch_size=batch_size, 
            record_limit=record_limit, 
            model_type=model_type,
            device=device
        )
        stats = processor.process_all()
        
        print(f"\n处理完成！结果保存在: {processor.output_dir}")
        print(f"设备: {stats['device']}")
        print(f"总配体数: {stats['total_ligands']}")
        print(f"成功嵌入数: {stats['successful_embeddings']}")
        print(f"失败数: {stats['failed_embeddings']}")
        print(f"M-L3关系数: {stats['m_l3_relationships']}")
        if stats['processing_time'] > 0:
            print(f"处理速度: {stats['total_ligands']/stats['processing_time']:.1f} 配体/秒")
        
    except Exception as e:
        logger.error(f"处理失败: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
