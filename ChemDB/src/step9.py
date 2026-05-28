#!/usr/bin/env python3
"""
ChemDB Step9 Processor - SMILES指纹生成
读取多个CSV文件中的SMILES字段，去重后生成Morgan指纹和DID索引
输出到./tmp/did_index.csv
使用流式处理避免内存瓶颈
"""

# Default configuration variables
INPUT_DIR = "./tmp"
DATA_DIR = "./data"
OUTPUT_DIR = "./tmp"
DEFAULT_RECORD_LIMIT = -1  # -1 for unlimited
DEFAULT_WORKERS = None  # Will use CPU count if None

import os
import json
import time
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit import RDLogger
import multiprocessing as mp
from multiprocessing import Process, Queue, cpu_count
from queue import Empty

# Disable RDKit logging
RDLogger.DisableLog('rdApp.*')
RDLogger.DisableLog('rdMol.*')
RDLogger.DisableLog('rdDecomposition.*')

from utils import calculate_did
from fingerprint_utils import generate_morgan_fingerprint, process_single_smiles

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def get_optimal_worker_count() -> int:
    """Get optimal worker count based on CPU cores"""
    cpu_cores = cpu_count()
    optimal_workers = min(cpu_cores, 200)
    logger.info(f"CPU核心数: {cpu_cores}, 使用CPU核心数作为工作进程数")
    logger.info(f"计算的最优工作进程数: {optimal_workers}")
    return optimal_workers

def process_single_smiles_wrapper(smiles_data: Tuple) -> Tuple[Optional[str], Optional[str], int]:
    """Wrapper for process_single_smiles to handle the tuple format"""
    try:
        smiles = smiles_data[0]
        return process_single_smiles(smiles)
    except Exception as e:
        logger.error(f"处理SMILES失败: {smiles_data}, 错误: {e}")
        return None, None, 0

def worker_process(task_queue: Queue, result_queue: Queue, worker_id: int, batch_size: int = 100):
    """Worker process for parallel SMILES processing"""
    processed_count = 0
    failed_count = 0
    batch_results = []
    
    while True:
        try:
            task = task_queue.get(timeout=1)
            
            if task is None:
                if batch_results:
                    for result in batch_results:
                        result_queue.put(result)
                break
            
            did, fingerprint, bit_count = process_single_smiles_wrapper(task)
            
            if did and fingerprint is not None:
                batch_results.append((did, fingerprint, bit_count))
            else:
                failed_count += 1
            
            processed_count += 1
            
            if len(batch_results) >= batch_size:
                for result in batch_results:
                    result_queue.put(result)
                batch_results = []
                
        except Empty:
            if batch_results:
                for result in batch_results:
                    result_queue.put(result)
                batch_results = []
            continue
        except Exception as e:
            logger.error(f"工作进程 {worker_id} 处理任务失败: {e}")
            failed_count += 1
            continue
    
    result_queue.put(('STATS', processed_count, failed_count))

class Step9Processor:
    def __init__(self, num_workers: int = None, record_limit: int = DEFAULT_RECORD_LIMIT):
        self.input_dir = Path(INPUT_DIR)
        self.data_dir = Path(DATA_DIR)
        self.output_dir = Path(OUTPUT_DIR)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        if num_workers is None:
            self.num_workers = get_optimal_worker_count()
        else:
            self.num_workers = num_workers
        
        self.record_limit = record_limit
        self.batch_size = max(10, 500 // self.num_workers)
        
        self.stats = {
            'total_smiles': 0,
            'unique_smiles': 0,
            'successful_fingerprints': 0,
            'failed_fingerprints': 0,
            'processing_time': 0,
            'num_workers': self.num_workers,
            'batch_size': self.batch_size,
            'record_limit': self.record_limit
        }
        
        self.output_file = self.output_dir / "did_index.csv"
        self.stats_output = self.output_dir / "step9_stats.json"
        
        logger.info(f"Step9处理器初始化完成，输入目录: {self.input_dir}")
        logger.info(f"输出目录: {self.output_dir}")
        logger.info(f"并行处理工作进程数: {self.num_workers}")
        logger.info(f"批处理大小: {self.batch_size}")
        logger.info(f"记录限制: {self.record_limit if self.record_limit > 0 else '无限制'}")
    
    def _get_total_smiles_count(self) -> int:
        """Get total count of SMILES from all files for progress tracking"""
        total_count = 0
        input_files = [
            self.input_dir / "complex_data.csv",
            self.input_dir / "ligand_data.csv", 
            self.input_dir / "repaired_ligand_data.csv",
            self.input_dir / "fragments.csv",
            self.data_dir / "IRL_filtered.csv",
        ]
        
        for file_path in input_files:
            if file_path.exists():
                try:
                    with open(file_path, 'r') as f:
                        next(f)  # Skip header
                        count = sum(1 for line in f)
                    total_count += count
                except:
                    continue
        
        return total_count

    def _producer_process(self, task_queue: Queue, num_workers: int):
        """Producer process that streams SMILES data to workers"""
        input_files = [
            self.input_dir / "complex_data.csv",
            self.input_dir / "ligand_data.csv", 
            self.input_dir / "repaired_ligand_data.csv",
            self.input_dir / "fragments.csv"
        ]
        
        logger.info("生产者进程开始流式加载数据...")
        
        chunk_size = 10000
        total_processed = 0
        all_smiles_seen = set()
        
        try:
            for file_path in input_files:
                if not file_path.exists():
                    logger.warning(f"文件不存在，跳过: {file_path}")
                    continue
                
                logger.info(f"从 {file_path.name} 流式提取SMILES...")
                file_smiles_count = 0
                
                for chunk in pd.read_csv(file_path, chunksize=chunk_size):
                    chunk_data = []
                    
                    for row in chunk.itertuples(index=False):
                        if self.record_limit > 0 and total_processed >= self.record_limit:
                            break
                        
                        # Check all columns for SMILES-like content
                        for col_name in chunk.columns:
                            if 'smiles' in col_name.lower():
                                value = getattr(row, col_name)
                                if pd.notna(value):
                                    smiles_str = str(value).strip()
                                    if smiles_str and smiles_str != 'nan' and smiles_str not in all_smiles_seen:
                                        all_smiles_seen.add(smiles_str)
                                        chunk_data.append(smiles_str)
                                        total_processed += 1
                                        file_smiles_count += 1
                        
                        if self.record_limit > 0 and total_processed >= self.record_limit:
                            break
                    
                    # Stream chunk data to workers
                    for smiles in chunk_data:
                        while task_queue.qsize() >= self.num_workers * 500:
                            time.sleep(0.01)
                        
                        task_queue.put((smiles,))
                    
                    if self.record_limit > 0 and total_processed >= self.record_limit:
                        break
                
                logger.info(f"文件 {file_path.name} 处理完成，新增: {file_smiles_count} 个唯一SMILES，当前总计: {len(all_smiles_seen)} 个")
            
            logger.info(f"生产者完成，总共处理了 {len(all_smiles_seen)} 个唯一SMILES")
            
        except Exception as e:
            logger.error(f"生产者进程错误: {e}")
        
        finally:
            for _ in range(num_workers):
                task_queue.put(None)
    
    def _init_csv_file(self, output_file: Path):
        """Initialize CSV file with headers"""
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write('did,fingerprint,ones\n')
        logger.info(f"CSV文件已初始化: {output_file}")
    
    def _append_to_csv(self, record: Dict, output_file: Path):
        """Append a single record to CSV file"""
        with open(output_file, 'a', encoding='utf-8') as f:
            f.write(f"{record['did']},{record['fingerprint']},{record['ones']}\n")

    def _save_stats(self):
        """Save processing statistics"""
        stats_data = {
            'timestamp': datetime.now().isoformat(),
            'processing_time_seconds': self.stats['processing_time'],
            'statistics': self.stats,
            'output_files': {
                'did_index': str(self.output_file)
            }
        }
        
        with open(self.stats_output, 'w', encoding='utf-8') as f:
            json.dump(stats_data, f, indent=2, ensure_ascii=False)
        
        logger.info(f"统计信息已保存到: {self.stats_output}")
    
    def process_all_smiles(self) -> Dict:
        """Main processing function"""
        start_time = time.time()
        logger.info("开始SMILES指纹生成处理")
        
        # Get total count for progress tracking
        total_count = self._get_total_smiles_count()
        logger.info(f"估计总SMILES数量: {total_count}")
        
        # Initialize CSV file
        self._init_csv_file(self.output_file)
        
        # Set up parallel processing with smaller queue sizes to prevent memory issues
        task_queue = Queue(maxsize=self.num_workers * 1000)
        result_queue = Queue(maxsize=self.num_workers * 1000)
        
        # Start worker processes
        logger.info(f"启动 {self.num_workers} 个工作进程...")
        processes = []
        for i in range(self.num_workers):
            p = Process(target=worker_process, args=(task_queue, result_queue, i, self.batch_size))
            p.start()
            processes.append(p)
        
        # Start producer process
        logger.info("启动生产者进程...")
        producer_process = Process(target=self._producer_process, args=(task_queue, self.num_workers))
        producer_process.start()
        
        # Stream results directly to CSV with deduplication
        seen_dids = set()  # In-memory set for deduplication
        total_failed = 0
        completed_tasks = 0
        
        with tqdm(desc="并行生成指纹") as pbar:
            while True:
                try:
                    result = result_queue.get(timeout=1)
                    
                    if result[0] == 'STATS':
                        _, processed, failed = result
                        total_failed += failed
                    else:
                        did, fingerprint, bit_count = result
                        
                        if did and fingerprint is not None:
                            # Check for duplicates
                            if did not in seen_dids:
                                seen_dids.add(did)
                                result_record = {
                                    'did': did,
                                    'fingerprint': fingerprint,
                                    'ones': bit_count
                                }
                                # Stream directly to CSV
                                self._append_to_csv(result_record, self.output_file)
                                self.stats['successful_fingerprints'] += 1
                        else:
                            total_failed += 1
                        
                        completed_tasks += 1
                        pbar.update(1)
                        
                        pbar.set_postfix({
                            '成功': self.stats['successful_fingerprints'],
                            '失败': total_failed,
                            '总计': completed_tasks,
                            '唯一': len(seen_dids)
                        })
                        
                except Empty:
                    # Check if all processes are done
                    if all(not p.is_alive() for p in processes) and not producer_process.is_alive():
                        break
                    # If producer is done but workers are still alive, continue waiting
                    if not producer_process.is_alive() and task_queue.empty():
                        # Give workers a bit more time to finish
                        time.sleep(0.1)
                    continue
        
        # Wait for all processes to finish
        producer_process.join()
        for p in processes:
            p.join()
        
        self.stats['total_smiles'] = completed_tasks
        self.stats['unique_smiles'] = len(seen_dids)
        self.stats['failed_fingerprints'] = total_failed
        
        total_time = time.time() - start_time
        self.stats['processing_time'] = total_time
        
        self._save_stats()
        
        logger.info("=" * 50)
        logger.info("SMILES指纹生成完成！最终统计:")
        logger.info(f"工作进程数: {self.stats['num_workers']}")
        logger.info(f"批处理大小: {self.stats['batch_size']}")
        logger.info(f"记录限制: {self.record_limit if self.record_limit > 0 else '无限制'}")
        logger.info(f"总SMILES数: {self.stats['total_smiles']}")
        logger.info(f"唯一SMILES数: {self.stats['unique_smiles']}")
        logger.info(f"成功生成指纹数: {self.stats['successful_fingerprints']}")
        logger.info(f"失败数: {self.stats['failed_fingerprints']}")
        logger.info(f"总处理时间: {total_time:.2f}秒")
        logger.info(f"平均处理速度: {self.stats['total_smiles']/total_time:.1f} SMILES/秒")
        logger.info(f"输出文件: {self.output_file}")
        logger.info("=" * 50)
        
        return self.stats

def main():
    """Main function"""
    try:
        import sys
        num_workers = None
        record_limit = DEFAULT_RECORD_LIMIT
        
        if len(sys.argv) > 1:
            try:
                num_workers = int(sys.argv[1])
                print(f"使用指定的工作进程数: {num_workers}")
            except ValueError:
                print("无效的工作进程数，使用默认值")
        
        if len(sys.argv) > 2:
            try:
                record_limit = int(sys.argv[2])
                print(f"使用指定的记录限制: {record_limit}")
            except ValueError:
                print("无效的记录限制，使用默认值")
        
        processor = Step9Processor(num_workers=num_workers, record_limit=record_limit)
        stats = processor.process_all_smiles()
        
        print(f"\n处理完成！结果保存在: {processor.output_file}")
        print(f"处理了 {stats['total_smiles']} 个SMILES")
        print(f"生成了 {stats['successful_fingerprints']} 个指纹")
        print(f"失败 {stats['failed_fingerprints']} 个")
        print(f"处理速度: {stats['total_smiles']/stats['processing_time']:.1f} SMILES/秒")
        
    except Exception as e:
        logger.error(f"处理失败: {e}")
        import sys
        sys.exit(1)

if __name__ == "__main__":
    main()
