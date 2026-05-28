#!/usr/bin/env python3
"""
ChemDB Step1 Processor - 独立版本
处理原始PubChem CSV数据，输出到./tmp目录
支持基于行的并行处理以提高性能
"""

# Default configuration variables
INPUT_DATA_DIR = "./data/pubchem"
METAL_LIST_PATH = "./data/metal_list.txt"
P_ELEMENTS_LIST_PATH = "./data/p_elements_list.txt"
DEFAULT_RECORD_LIMIT = -1  # -1 for unlimited
OUTPUT_DIR = "./tmp"
DEFAULT_WORKERS = None  # Resolved via _resolve_num_workers()
DEFAULT_MAX_WORKERS = 200  # Full-run default cap; override via --workers or CHEMDB_STEP1_WORKERS
SAVE_PUBCHEM_DATA = True  # Set to False to disable pubchem_data.csv output
SORT_COMPLEX_BY_SMILES_LENGTH = False  # Set to True to sort complex by SMILES length (longer first)


def _sort_complex_by_smiles_length() -> bool:
    """Env CHEMDB_STEP1_SORT_BY_LENGTH=1 enables sort (used by integration tests)."""
    if os.environ.get("CHEMDB_STEP1_SORT_BY_LENGTH", "").lower() in ("1", "true", "yes"):
        return True
    return SORT_COMPLEX_BY_SMILES_LENGTH


def _row_timeout_sec() -> int:
    raw = os.environ.get("CHEMDB_STEP1_ROW_TIMEOUT_SEC", "0")
    try:
        return max(0, int(raw))
    except ValueError:
        return 0

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
from multiprocessing import Process, Queue, Manager
from queue import Empty
from tools.CC_split import ChemicalComplex
            
import threading

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
RDLogger.DisableLog('rdMol.*')
RDLogger.DisableLog('rdDecomposition.*')

from utils import calculate_did, is_coordination_complex, load_elements_list

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def _process_single_row_impl(row_data: Tuple, metal_list: List[str]) -> Tuple[Optional[Dict], bool]:
    try:
        cid, smiles, source_file, row_index = row_data
        
        if not smiles or smiles == 'nan' or smiles == '':
            return None, False
        
        if is_coordination_complex(smiles, metal_list):
            try:
                # 创建ChemicalComplex实例进行分析
                complex_info = {
                    'cid': cid,
                    'smiles': smiles,
                    'corcomp': True
                }
                complex_instance = ChemicalComplex(complex_info, metal_list)
                
                # 获取完整的结构信息
                result = complex_instance.merge_neighbor_info()
                
                if result and result.get('complex_info', {}).get('DID'):
                    # Return the complete complex_info structure as provided by merge_neighbor_info()
                    return result, False
                    
            except Exception as e:
                return None, True
                    
        
        return None, False
        
    except Exception as e:
        return None, True


def process_single_row(row_data: Tuple, metal_list: List[str]) -> Tuple[Optional[Dict], bool]:
    timeout = _row_timeout_sec()
    if timeout <= 0:
        return _process_single_row_impl(row_data, metal_list)

    cid = row_data[0] if row_data else "?"
    out: Dict[str, object] = {"result": (None, True)}

    def _target() -> None:
        out["result"] = _process_single_row_impl(row_data, metal_list)

    thread = threading.Thread(target=_target, daemon=True)
    thread.start()
    thread.join(timeout)
    if thread.is_alive():
        logger.warning("行处理超时 (%ss)，跳过 CID=%s", timeout, cid)
        return None, True
    return out["result"]  # type: ignore[return-value]


def deduplicate_csv_data(all_rows: List[Tuple]) -> Tuple[List[Tuple], List[Tuple]]:
    """
    Deduplicate CSV data by CID immediately after loading.
    Returns (all_pubchem_rows, unique_complex_rows) where:
    - all_pubchem_rows: all rows for pubchem data output
    - unique_complex_rows: deduplicated rows for complex/ligand processing
    """
    seen_cids = set()
    unique_complex_rows = []
    
    for row_data in all_rows:
        cid = row_data[0]  # CID is the first element in the tuple
        if cid not in seen_cids:
            seen_cids.add(cid)
            unique_complex_rows.append(row_data)
    
    logger.info(f"CSV数据去重结果: {len(all_rows)} -> {len(unique_complex_rows)} 唯一CID")
    
    return all_rows, unique_complex_rows

def worker_process(task_queue: Queue, result_queue: Queue, metal_list: List[str], worker_id: int):
    processed_count = 0
    failed_count = 0
    task = None

    while True:
        try:
            task = task_queue.get(timeout=1)

            if task is None:
                break

            complex_record, is_failed = process_single_row(task, metal_list)
            result_queue.put((complex_record, is_failed))

            processed_count += 1
            if is_failed:
                failed_count += 1
            task = None

        except Empty:
            continue
        except Exception as e:
            logger.error(f"工作进程 {worker_id} 处理任务失败: {e}")
            failed_count += 1
            if task is not None:
                result_queue.put((None, True))
                processed_count += 1
            task = None

    result_queue.put(('STATS', processed_count, failed_count))


def _resolve_num_workers(num_workers: Optional[int]) -> int:
    """CLI --workers > CHEMDB_STEP1_WORKERS (integration) > min(cpu, 200) for full runs."""
    if num_workers is not None:
        return max(1, int(num_workers))
    env = os.environ.get("CHEMDB_STEP1_WORKERS")
    if env:
        return max(1, int(env))
    return max(1, min(mp.cpu_count(), DEFAULT_MAX_WORKERS))


def load_all_csv_data(csv_files: List[Path], record_limit: int = DEFAULT_RECORD_LIMIT) -> List[Tuple]:
    all_rows = []

    for csv_file in tqdm(csv_files, desc="加载CSV文件"):
        try:
            if record_limit > 0 and len(all_rows) >= record_limit:
                break

            df = pd.read_csv(csv_file, dtype=str, low_memory=False, encoding="utf-8-sig")
            df.columns = df.columns.str.strip()
            if "cid" not in df.columns or "isosmiles" not in df.columns:
                raise ValueError(f"missing cid/isosmiles columns: {list(df.columns[:5])}")
            df = df[["cid", "isosmiles"]]

            for row in df.itertuples(index=True):
                if record_limit > 0 and len(all_rows) >= record_limit:
                    break

                cid = row.cid.strip() if row.cid else f'unknown_{row.Index}'
                smiles = row.isosmiles

                if smiles and smiles != 'nan' and smiles.strip():
                    all_rows.append((cid, smiles, csv_file.name, row.Index))

        except Exception as e:
            logger.error(f"加载文件 {csv_file} 失败: {e}")

    return all_rows

class Step1Processor:
    def __init__(
        self,
        num_workers: int = None,
        record_limit: int = DEFAULT_RECORD_LIMIT,
        pubchem_dir: str = None,
        metal_list_path: str = None,
        p_elements_list_path: str = None,
        tmp_dir: str = None,
    ):
        self.output_dir = Path(tmp_dir or OUTPUT_DIR)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.pubchem_dir = Path(pubchem_dir or INPUT_DATA_DIR)
        self.metal_list_path = metal_list_path or METAL_LIST_PATH
        self.p_elements_list_path = p_elements_list_path or P_ELEMENTS_LIST_PATH
        
        self.num_workers = _resolve_num_workers(num_workers)
        
        self.record_limit = record_limit
        self.metal_list = load_elements_list(self.metal_list_path)
        
        self.stats = {
            'total_compounds': 0,
            'deduplicated_compounds': 0,
            'coordination_compounds': 0,
            'ligands': 0,
            'failed_compounds': 0,
            'processing_time': 0,
            'num_workers': self.num_workers,
            'record_limit': self.record_limit
        }
        
        self.pubchem_output = self.output_dir / "pubchem_data.csv"
        self.complex_output = self.output_dir / "complex_data.csv"
        self.ligand_output = self.output_dir / "ligand_data.csv"
        self.stats_output = self.output_dir / "processing_stats.json"
        
        logger.info(f"Step1处理器初始化完成，输出目录: {self.output_dir}")
        logger.info(f"并行处理工作进程数: {self.num_workers}")
        logger.info(f"记录限制: {self.record_limit if self.record_limit > 0 else '无限制'}")
    
    def _save_data_to_csv(self, data: List[Dict], output_file: Path):
        if not data:
            logger.info(f"没有数据需要保存到: {output_file}")
            return
            
        df = pd.DataFrame(data)
        df.to_csv(output_file, index=False, encoding='utf-8')
        logger.info(f"数据已保存到: {output_file} ({len(data)} 条记录)")
    
    def _save_stats(self):
        output_files = {
            'complex_data': str(self.complex_output),
            'ligand_data': str(self.ligand_output)
        }
        
        if SAVE_PUBCHEM_DATA:
            output_files['pubchem_data'] = str(self.pubchem_output)
        
        stats_data = {
            'timestamp': datetime.now().isoformat(),
            'processing_time_seconds': self.stats['processing_time'],
            'statistics': self.stats,
            'output_files': output_files
        }
        
        with open(self.stats_output, 'w', encoding='utf-8') as f:
            json.dump(stats_data, f, indent=2, ensure_ascii=False)
        
        logger.info(f"统计信息已保存到: {self.stats_output}")
    
    def process_all_files(self) -> Dict:
        start_time = time.time()
        logger.info("开始基于行的并行处理")
        
        input_dir = self.pubchem_dir
        if not input_dir.exists():
            raise FileNotFoundError(f"输入目录不存在: {input_dir}")
        
        csv_files = list(input_dir.glob("*.csv"))
        logger.info(f"找到 {len(csv_files)} 个CSV文件")
        
        logger.info("开始加载所有CSV文件数据...")
        all_rows = load_all_csv_data(csv_files, self.record_limit)
        logger.info(f"加载了 {len(all_rows)} 行数据")
        
        if not all_rows:
            logger.warning("没有找到有效的数据行")
            return self.stats
        
        # Deduplicate CSV data immediately after loading
        logger.info("开始CSV数据去重...")
        all_pubchem_rows, unique_complex_rows = deduplicate_csv_data(all_rows)
        
        # Create pubchem records from all rows (for pubchem_data.csv if enabled)
        all_pubchem_data = []
        for row_data in all_pubchem_rows:
            cid, smiles, source_file, row_index = row_data
            pubchem_record = {
                'cid': cid,
                'smiles': smiles,
                'source_file': source_file,
                'row_index': row_index
            }
            all_pubchem_data.append(pubchem_record)
        
        logger.info(f"PubChem记录数: {len(all_pubchem_data)}")
        logger.info(f"去重后待处理记录数: {len(unique_complex_rows)}")
        
        if _sort_complex_by_smiles_length():
            # Sort by SMILES length (longer first) for better parallel processing performance
            logger.info("按SMILES长度排序（较长优先）...")
            unique_complex_rows.sort(key=lambda x: len(x[1]), reverse=True)  # x[1] is the SMILES string
            logger.info("排序完成")
        
        # Update stats with deduplication info
        self.stats['total_compounds'] = len(all_pubchem_data)
        self.stats['deduplicated_compounds'] = len(unique_complex_rows)
        
        task_queue = Queue()
        result_queue = Queue()
        
        # Only process unique complex rows for complex/ligand analysis
        for row_data in unique_complex_rows:
            task_queue.put(row_data)
        
        for _ in range(self.num_workers):
            task_queue.put(None)
        
        logger.info(f"启动 {self.num_workers} 个工作进程...")
        processes = []
        for i in range(self.num_workers):
            p = Process(target=worker_process, args=(task_queue, result_queue, self.metal_list, i))
            p.start()
            processes.append(p)
        
        all_complex_data = []
        all_ligand_data = []
        total_failed = 0
        completed_tasks = 0
        last_progress_at = time.time()
        stall_log_interval = 120.0

        with tqdm(total=len(unique_complex_rows), desc="并行处理数据行") as pbar:
            while completed_tasks < len(unique_complex_rows):
                try:
                    result = result_queue.get(timeout=1)
                    
                    if result[0] == 'STATS':
                        _, processed, failed = result
                        total_failed += failed
                    else:
                        complex_info, is_failed = result
                        
                        if complex_info:
                            # Extract complex data for CSV output
                            complex_did = complex_info['complex_info'].get('DID')
                            complex_smiles = complex_info['complex_info'].get('complex_smiles', '')
                            inactive = complex_info['complex_info'].get('inactive', 0)
                            
                            # Generate metal_info from central_metal_info
                            metal_info_set = set()
                            for metal in complex_info.get('central_metal_info', []):
                                symbol = metal.get('central_metal')
                                valence = metal.get('valence')
                                if symbol:
                                    val_str = f"'{valence}'" if valence is not None else 'NULL'
                                    metal_info_set.add(f"{{'{symbol}': {val_str}}}")
                            metal_info_json = '[' + ','.join(metal_info_set) + ']'
                            
                            # Create complex record for CSV
                            complex_record = {
                                'did': complex_did,
                                'cid': complex_info.get('cid', ''),
                                'complex_smiles': complex_smiles,
                                'metal_info': metal_info_json,
                                'ligand_count': len(complex_info.get('neighbor_info', [])),
                                'from_pubchem': True,
                                'inactive': inactive
                            }
                            all_complex_data.append(complex_record)
                            
                            # Extract ligand data from neighbor_info
                            for i, ligand_info in enumerate(complex_info.get('neighbor_info', [])):
                                ligand_record = {
                                    'ligand_did': ligand_info.get('ligand_did', ''),
                                    'ligand_smiles': ligand_info.get('ligand_smiles', ''),
                                    # 'ligand_inactive': 0,
                                    'source_complex_did': complex_did,
                                    'ligand_index': i,
                                    # 'coordinating_atom_index': ligand_info.get('coordinating_atom_index'),
                                    # 'bond_type': str(ligand_info.get('bond_type', 'UNKNOWN')),
                                    # 'central_atom': ligand_info.get('central_atom', ''),
                                }
                                all_ligand_data.append(ligand_record)
                        
                        if is_failed:
                            total_failed += 1
                        
                        completed_tasks += 1
                        pbar.update(1)
                        last_progress_at = time.time()

                except Empty:
                    alive = [p for p in processes if p.is_alive()]
                    if not alive:
                        missing = len(unique_complex_rows) - completed_tasks
                        logger.error(
                            "并行主进程结束但结果不完整: %s/%s 已完成，%s 条无 worker 回传",
                            completed_tasks,
                            len(unique_complex_rows),
                            missing,
                        )
                        total_failed += missing
                        break
                    now = time.time()
                    if now - last_progress_at >= stall_log_interval:
                        logger.warning(
                            "并行处理 stall: %s/%s 已完成，%s 个 worker 仍存活",
                            completed_tasks,
                            len(unique_complex_rows),
                            len(alive),
                        )
                        last_progress_at = now
                    continue
        
        for p in processes:
            p.join()
        
        # Update stats with final counts
        self.stats['coordination_compounds'] = len(all_complex_data)
        self.stats['ligands'] = len(all_ligand_data)
        self.stats['failed_compounds'] = total_failed
        
        logger.info("开始保存处理结果")
        
        # Sort results for consistent output order
        if all_pubchem_data and SAVE_PUBCHEM_DATA:
            all_pubchem_data.sort(key=lambda x: (x.get('cid', ''), x.get('row_index', 0)))
            logger.info("PubChem数据已按CID和row_index排序")
        
        if all_complex_data:
            all_complex_data.sort(key=lambda x: (x.get('did', ''), x.get('cid', '')))
            logger.info("复合物数据已按DID和CID排序")
        
        if all_ligand_data:
            all_ligand_data.sort(key=lambda x: x.get('ligand_did', ''))
            logger.info("配体数据已按ligand_did排序")
        
        # Conditionally save pubchem data based on SAVE_PUBCHEM_DATA constant
        if SAVE_PUBCHEM_DATA:
            self._save_data_to_csv(all_pubchem_data, self.pubchem_output)
        else:
            logger.info("跳过PubChem数据保存 (SAVE_PUBCHEM_DATA = False)")
        
        self._save_data_to_csv(all_complex_data, self.complex_output)
        self._save_data_to_csv(all_ligand_data, self.ligand_output)
        
        total_time = time.time() - start_time
        self.stats['processing_time'] = total_time
        
        self._save_stats()
        
        logger.info("=" * 50)
        logger.info("基于行的并行处理完成！最终统计:")
        logger.info(f"工作进程数: {self.stats['num_workers']}")
        logger.info(f"记录限制: {self.record_limit if self.record_limit > 0 else '无限制'}")
        logger.info(f"总化合物数: {self.stats['total_compounds']}")
        logger.info(f"去重后待处理化合物数: {self.stats['deduplicated_compounds']}")
        logger.info(f"配位化合物数: {self.stats['coordination_compounds']}")
        logger.info(f"配体数: {self.stats['ligands']}")
        logger.info(f"失败化合物数: {self.stats['failed_compounds']}")
        logger.info(f"总处理时间: {total_time:.2f}秒")
        logger.info(f"平均处理速度: {self.stats['deduplicated_compounds']/total_time:.1f} 化合物/秒")
        logger.info(f"输出目录: {self.output_dir}")
        logger.info(f"PubChem数据保存: {'启用' if SAVE_PUBCHEM_DATA else '禁用'}")
        logger.info("=" * 50)
        
        return self.stats

def main():
    try:
        import argparse
        import sys

        # Legacy: python step1.py [workers] [record_limit]
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
            processor = Step1Processor(num_workers=num_workers, record_limit=record_limit)
            stats = processor.process_all_files()
        else:
            parser = argparse.ArgumentParser(description="ChemDB Step1")
            parser.add_argument("--pubchem-dir", type=str, default=INPUT_DATA_DIR)
            parser.add_argument("--metal-list", type=str, default=METAL_LIST_PATH)
            parser.add_argument("--p-elements-list", type=str, default=P_ELEMENTS_LIST_PATH)
            parser.add_argument("--tmp-dir", type=str, default=OUTPUT_DIR)
            parser.add_argument("--workers", type=int, default=None)
            parser.add_argument("--limit", type=int, default=DEFAULT_RECORD_LIMIT)
            args = parser.parse_args()
            processor = Step1Processor(
                num_workers=args.workers,
                record_limit=args.limit,
                pubchem_dir=args.pubchem_dir,
                metal_list_path=args.metal_list,
                p_elements_list_path=args.p_elements_list,
                tmp_dir=args.tmp_dir,
            )
            stats = processor.process_all_files()
        
        print(f"\n处理完成！结果保存在: {processor.output_dir}")
        print(f"总化合物数: {stats['total_compounds']}")
        print(f"去重后待处理化合物数: {stats['deduplicated_compounds']}")
        print(f"发现 {stats['coordination_compounds']} 个配位化合物")
        print(f"提取了 {stats['ligands']} 个配体")
        print(f"处理速度: {stats['deduplicated_compounds']/stats['processing_time']:.1f} 化合物/秒")
        print(f"PubChem数据保存: {'启用' if SAVE_PUBCHEM_DATA else '禁用'}")
        
    except Exception as e:
        logger.error(f"处理失败: {e}")
        import sys
        sys.exit(1)

if __name__ == "__main__":
    main() 