#!/usr/bin/env python3
"""
ChemDB Step2 Processor - 配体结构修复
基于step1的输出，修复有问题的SMILES结构，输出到./tmp目录
支持基于行的并行处理以提高性能，充分利用所有CPU核心
"""

# Default configuration variables
INPUT_DATA_DIR = "./data/pubchem"
METAL_LIST_PATH = "./data/metal_list.txt"
P_ELEMENTS_LIST_PATH = "./data/p_elements_list.txt"
DEFAULT_RECORD_LIMIT = -1  # -1 for unlimited
OUTPUT_DIR = "./tmp"
INPUT_DIR = "./tmp"
DEFAULT_WORKERS = None  # Will use CPU count if None

import os
import json
import time
import yaml
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
import hashlib
import multiprocessing as mp
from multiprocessing import Process, Queue, Manager, cpu_count
from queue import Empty
import threading

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
RDLogger.DisableLog('rdMol.*')
RDLogger.DisableLog('rdDecomposition.*')

from utils import calculate_did, load_elements_list

from tools.fix_smiles import MolecularFormulaRepairer

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def get_optimal_worker_count() -> int:
    """CHEMDB_STEP2_WORKERS (integration) > min(cpu, 200) for full runs."""
    env = os.environ.get("CHEMDB_STEP2_WORKERS")
    if env:
        n = max(1, int(env))
        logger.info(f"使用 CHEMDB_STEP2_WORKERS={n}")
        return n
    cpu_cores = cpu_count()
    optimal_workers = min(cpu_cores, 200)
    logger.info(f"CPU核心数: {cpu_cores}, 使用CPU核心数作为工作进程数")
    logger.info(f"计算的最优工作进程数: {optimal_workers}")
    return optimal_workers

def process_single_ligand(ligand_data: Tuple) -> Tuple[Optional[str], str, str, int]:
    try:
        did, ligand_smiles = ligand_data
        
        if not ligand_smiles or ligand_smiles == 'nan' or ligand_smiles == '':
            return "INVALID", did, "", 0
        
        repairer = MolecularFormulaRepairer(ligand_smiles, debug=False, draw=False)
        new_smiles, need_kill, readable = repairer.repair_smiles()
        
        if need_kill or not readable:
            return "INVALID", did, "", 0
        
        if new_smiles == ligand_smiles:
            return did, did, new_smiles, 2
        else:
            new_did = calculate_did(new_smiles)
            if new_did:
                return new_did, did, new_smiles, 1
            else:
                return "INVALID", did, "", 0
                
    except Exception as e:
        logger.error(f"处理配体失败: {ligand_data}, 错误: {e}")
        return "INVALID", ligand_data[0], "", 0

def worker_process(task_queue: Queue, result_queue: Queue, worker_id: int, batch_size: int = 100):
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
            
            new_did, source_did, new_smiles, is_repaired = process_single_ligand(task)
            
            batch_results.append((new_did, source_did, new_smiles, task[1], is_repaired))
            
            processed_count += 1
            if is_repaired == 0:
                failed_count += 1
            
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

class Step2Processor:
    def __init__(
        self,
        num_workers: int = None,
        record_limit: int = DEFAULT_RECORD_LIMIT,
        input_dir: str = None,
        output_dir: str = None,
    ):
        self.input_dir = Path(input_dir or INPUT_DIR)
        self.output_dir = Path(output_dir or OUTPUT_DIR)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        if num_workers is None:
            self.num_workers = get_optimal_worker_count()
        else:
            self.num_workers = num_workers
        
        self.record_limit = record_limit
        self.batch_size = max(50, 1000 // self.num_workers)
        
        self.stats = {
            'total_ligands': 0,
            'repaired_ligands': 0,
            'unchanged_ligands': 0,
            'invalid_ligands': 0,
            'processing_time': 0,
            'num_workers': self.num_workers,
            'batch_size': self.batch_size,
            'record_limit': self.record_limit
        }
        
        self.repaired_output = self.output_dir / "repaired_ligand_data.csv"
        self.stats_output = self.output_dir / "repair_stats.json"
        
        logger.info(f"Step2处理器初始化完成，输出目录: {self.output_dir}")
        logger.info(f"并行处理工作进程数: {self.num_workers}")
        logger.info(f"批处理大小: {self.batch_size}")
        logger.info(f"记录限制: {self.record_limit if self.record_limit > 0 else '无限制'}")
    
    def _load_ligand_data(self) -> List[Tuple]:
        ligand_file = self.input_dir / "ligand_data.csv"
        if not ligand_file.exists():
            raise FileNotFoundError(f"Step1输出文件不存在: {ligand_file}")
        
        logger.info(f"从 {ligand_file} 加载配体数据...")
        
        chunk_size = 10000
        ligand_data = []
        
        for chunk in pd.read_csv(ligand_file, chunksize=chunk_size):
            for _, row in chunk.iterrows():
                if self.record_limit > 0 and len(ligand_data) >= self.record_limit:
                    break
                    
                did = row.get('did', '')
                ligand_smiles = row.get('ligand_smiles', '')
                if did and ligand_smiles:
                    ligand_data.append((did, ligand_smiles))
                    
                if self.record_limit > 0 and len(ligand_data) >= self.record_limit:
                    break
        
        logger.info(f"加载了 {len(ligand_data)} 个配体")
        return ligand_data
    
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
                'repaired_ligand_data': str(self.repaired_output)
            }
        }
        
        with open(self.stats_output, 'w', encoding='utf-8') as f:
            json.dump(stats_data, f, indent=2, ensure_ascii=False)
        
        logger.info(f"统计信息已保存到: {self.stats_output}")
    
    def process_all_ligands(self) -> Dict:
        start_time = time.time()
        logger.info("开始配体结构修复处理")
        
        task_queue = Queue(maxsize=self.num_workers * 5000)
        result_queue = Queue(maxsize=self.num_workers * 5000)
        
        logger.info(f"启动 {self.num_workers} 个工作进程...")
        processes = []
        for i in range(self.num_workers):
            p = Process(target=worker_process, args=(task_queue, result_queue, i, self.batch_size))
            p.start()
            processes.append(p)
        
        logger.info("开始流式加载和处理数据...")
        
        producer_process = Process(target=self._producer_process, args=(task_queue, self.num_workers))
        producer_process.start()
        
        all_repaired_data = []
        valid_repaired_data = []
        total_failed = 0
        completed_tasks = 0
        
        total_count = self._get_total_ligand_count()
        
        with tqdm(total=total_count, desc="并行修复配体") as pbar:
            while completed_tasks < total_count:
                try:
                    result = result_queue.get(timeout=1)
                    
                    if result[0] == 'STATS':
                        _, processed, failed = result
                        total_failed += failed
                    else:
                        new_did, source_did, new_smiles, old_smiles, is_repaired = result
                        
                        repaired_record = {
                            'ligand_new_did': new_did,
                            'source_did': source_did,
                            'new_smiles': new_smiles,
                            'old_smiles': old_smiles,
                            'is_repaired': is_repaired
                        }
                        all_repaired_data.append(repaired_record)
                        
                        if is_repaired > 0:  # 1 for repaired, 2 for unchanged
                            valid_repaired_data.append(repaired_record)
                        
                        if is_repaired == 0:
                            total_failed += 1
                        
                        completed_tasks += 1
                        pbar.update(1)
                        
                        pbar.set_postfix({
                            '成功': completed_tasks - total_failed,
                            '失败': total_failed,
                            '进度': f"{completed_tasks}/{total_count}"
                        })
                        
                except Empty:
                    if all(not p.is_alive() for p in processes) and not producer_process.is_alive():
                        break
                    continue
        
        producer_process.join()
        for p in processes:
            p.join()
        
        self.stats['total_ligands'] = len(all_repaired_data)
        self.stats['repaired_ligands'] = len([r for r in all_repaired_data if r['is_repaired'] == 1])
        self.stats['unchanged_ligands'] = len([r for r in all_repaired_data if r['is_repaired'] == 2])
        self.stats['invalid_ligands'] = len([r for r in all_repaired_data if r['is_repaired'] == 0])
        
        logger.info("开始保存修复结果")
        self._save_data_to_csv(valid_repaired_data, self.repaired_output)
        
        total_time = time.time() - start_time
        self.stats['processing_time'] = total_time
        
        self._save_stats()
        
        logger.info("=" * 50)
        logger.info("配体结构修复完成！最终统计:")
        logger.info(f"工作进程数: {self.stats['num_workers']}")
        logger.info(f"批处理大小: {self.stats['batch_size']}")
        logger.info(f"记录限制: {self.record_limit if self.record_limit > 0 else '无限制'}")
        logger.info(f"总配体数: {self.stats['total_ligands']}")
        logger.info(f"修复配体数: {self.stats['repaired_ligands']}")
        logger.info(f"未变化配体数: {self.stats['unchanged_ligands']}")
        logger.info(f"无效配体数: {self.stats['invalid_ligands']}")
        logger.info(f"有效数据保存数: {len(valid_repaired_data)}")  # New log line
        logger.info(f"总处理时间: {total_time:.2f}秒")
        logger.info(f"平均处理速度: {self.stats['total_ligands']/total_time:.1f} 配体/秒")
        logger.info(f"输出目录: {self.output_dir}")
        logger.info("=" * 50)
        
        return self.stats

    def _get_total_ligand_count(self) -> int:
        ligand_file = self.input_dir / "ligand_data.csv"
        if not ligand_file.exists():
            return 0
        
        try:
            with open(ligand_file, 'r') as f:
                next(f)
                count = sum(1 for line in f)
            return count
        except:
            return 0

    def _producer_process(self, task_queue: Queue, num_workers: int):
        ligand_file = self.input_dir / "ligand_data.csv"
        if not ligand_file.exists():
            return
        
        logger.info("生产者进程开始加载数据...")
        
        chunk_size = 10000
        total_processed = 0
        
        try:
            for chunk in pd.read_csv(ligand_file, chunksize=chunk_size):
                chunk_data = []
                
                for _, row in chunk.iterrows():
                    if self.record_limit > 0 and total_processed >= self.record_limit:
                        break
                        
                    did = row.get('ligand_did', '')
                    ligand_smiles = row.get('ligand_smiles', '')
                    if did and ligand_smiles:
                        chunk_data.append((did, ligand_smiles))
                        total_processed += 1
                        
                    if self.record_limit > 0 and total_processed >= self.record_limit:
                        break
                
                for item in chunk_data:
                    while task_queue.qsize() >= 4000:
                        time.sleep(0.01)
                    
                    task_queue.put(item)
                
                if self.record_limit > 0 and total_processed >= self.record_limit:
                    break
            
            logger.info(f"生产者完成，总共处理了 {total_processed} 个配体")
            
        except Exception as e:
            logger.error(f"生产者进程错误: {e}")
        
        finally:
            for _ in range(num_workers):
                task_queue.put(None)

def main():
    try:
        import argparse
        import sys

        if len(sys.argv) > 1 and not any(a.startswith("-") for a in sys.argv[1:]):
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
            processor = Step2Processor(num_workers=num_workers, record_limit=record_limit)
            stats = processor.process_all_ligands()
        else:
            parser = argparse.ArgumentParser(description="ChemDB Step2")
            parser.add_argument("--input-dir", type=str, default=INPUT_DIR)
            parser.add_argument("--output-dir", type=str, default=OUTPUT_DIR)
            parser.add_argument("--workers", type=int, default=None)
            parser.add_argument("--limit", type=int, default=DEFAULT_RECORD_LIMIT)
            args = parser.parse_args()
            processor = Step2Processor(
                num_workers=args.workers,
                record_limit=args.limit,
                input_dir=args.input_dir,
                output_dir=args.output_dir,
            )
            stats = processor.process_all_ligands()
        
        print(f"\n处理完成！结果保存在: {processor.output_dir}")
        print(f"处理了 {stats['total_ligands']} 个配体")
        print(f"修复了 {stats['repaired_ligands']} 个配体")
        print(f"未变化 {stats['unchanged_ligands']} 个配体")
        print(f"无效 {stats['invalid_ligands']} 个配体")
        print(f"处理速度: {stats['total_ligands']/stats['processing_time']:.1f} 配体/秒")
        
    except Exception as e:
        logger.error(f"处理失败: {e}")
        import sys
        sys.exit(1)

if __name__ == "__main__":
    main() 