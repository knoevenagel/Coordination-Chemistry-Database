#!/usr/bin/env python3
"""
Step6+7处理器：Ligand标记和L4片段生成合并版
"""

# Default configuration variables
INPUT_DATA_DIR = "./data/pubchem"
METAL_LIST_PATH = "./data/metal_list.txt"
P_ELEMENTS_LIST_PATH = "./data/p_elements_list.txt"
GA_PATH = "./tmp/GA_with_id.csv"
IRL_PATH = "./tmp/IRL_filtered.csv"
DEFAULT_RECORD_LIMIT = -1  # -1 for unlimited
OUTPUT_DIR = "./tmp"
INPUT_DIR = "./tmp"
DEFAULT_WORKERS = None  # Will use CPU count if None

import os
import sys
import yaml
import json
import logging
import traceback
import multiprocessing as mp
from datetime import datetime
from pathlib import Path
import time
from multiprocessing import Process, Queue, Manager, cpu_count
from queue import Empty
from typing import Dict, List, Tuple, Optional, Set
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
#from tools.L4_create import MoleculeSplitter
from tools.L4_create_independent import IndependentL4Extractor as MoleculeSplitter

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
RDLogger.DisableLog('rdMol.*')
RDLogger.DisableLog('rdDecomposition.*')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def get_optimal_worker_count() -> int:
    cpu_cores = cpu_count()
    optimal_workers = min(cpu_cores, 200)
    logger.info(f"CPU核心数: {cpu_cores}, 使用CPU核心数作为工作进程数")
    logger.info(f"计算的最优工作进程数: {optimal_workers}")
    return optimal_workers

class MoleculeMarker:
    def __init__(self, ga_data: Dict[str, str], irl_data: Dict[str, str]):
        self.ga_data = ga_data
        self.irl_data = irl_data
    
    def mark_smiles(self, smiles: str) -> Optional[Chem.Mol]:
        try:
            input_mol = Chem.MolFromSmiles(smiles)
            if input_mol is None:
                return None
            
            for atom in input_mol.GetAtoms():
                atom.SetProp("atomNote", str(atom.GetIdx()))
            
            self.mark_IRL(input_mol)
            self.mark_GA(input_mol)
            
            return input_mol
        except Exception as e:
            logger.error(f"标记分子失败: {smiles}, 错误: {e}")
            return None
    
    def mark_IRL(self, molecule: Chem.Mol):
        for irl_smiles, irl_id in self.irl_data.items():
            try:
                irl_mol = Chem.MolFromSmiles(irl_smiles)
                if irl_mol:
                    matches = molecule.GetSubstructMatches(irl_mol)
                    for match in matches:
                        self._add_label_to_matched_atoms(molecule, irl_id.upper(), "IRL", match)
            except Exception as e:
                logger.warning(f"标记IRL失败: {irl_smiles}, 错误: {e}")
    
    def mark_GA(self, molecule: Chem.Mol):
        for ga_smiles, ga_id in self.ga_data.items():
            try:
                ga_mol = Chem.MolFromSmiles(ga_smiles)
                if ga_mol:
                    matches = molecule.GetSubstructMatches(ga_mol)
                    for match in matches:
                        self._add_label_to_matched_atoms(molecule, ga_id, "GA", match)
            except Exception as e:
                logger.warning(f"标记GA失败: {ga_smiles}, 错误: {e}")
    
    def _add_label_to_matched_atoms(self, molecule: Chem.Mol, atom_id: str, label_type: str, match: Tuple[int, ...]):
        for atom_idx in match:
            atom = molecule.GetAtomWithIdx(atom_idx)
            current = atom.GetPropsAsDict().get(f"{label_type}_ID", "")
            existing = current.split(";") if current else []
            if atom_id not in existing:
                existing.append(atom_id)
            atom.SetProp(f"{label_type}_ID", ";".join(existing))
    
    def get_marked_data(self, molecule: Chem.Mol) -> Dict:
        marked_data = {}
        for atom in molecule.GetAtoms():
            atom_idx = atom.GetIdx()
            atom_props = atom.GetPropsAsDict()
            
            marked_data[atom_idx] = {
                'element': atom.GetSymbol(),
                'IRL_ID': atom_props.get('IRL_ID', ''),
                'GA_ID': atom_props.get('GA_ID', ''),
                'atomNote': atom_props.get('atomNote', '')
            }
        
        return marked_data

def process_single_ligand(ligand_data: Tuple, marker: MoleculeMarker, irl_list: List[Dict], ga_data_dict: Dict) -> Tuple[bool, str, str, List[Dict]]:
    try:
        did, smiles = ligand_data
        
        marked_mol = marker.mark_smiles(smiles)
        if marked_mol is None:
            error_msg = f"分子解析失败: DID={did}, SMILES={smiles}"
            logger.error(error_msg)
            return False, did, error_msg, []
        
        marked_data = marker.get_marked_data(marked_mol)
        
        formatted_marked_data = {}
        for atom_idx, atom_info in marked_data.items():
            atom_key = f'atom_{atom_idx}'
            formatted_marked_data[atom_key] = {
                'IRL_ids': atom_info['IRL_ID'].split(';') if atom_info['IRL_ID'] else [],
                'ga_ids': atom_info['GA_ID'].split(';') if atom_info['GA_ID'] else [],
                'atom_symbol': atom_info['element']
            }
        
        try:
            try:
                splitter = MoleculeSplitter(
                    smiles=smiles,
                    marked_ligand_data=formatted_marked_data,
                    irl_list=irl_list,
                    ga_list=list(ga_data_dict.keys()),
                    source_did=did
                )
            except Exception as e:
                error_msg = f"创建MoleculeSplitter实例失败: {e}"
                logger.error(f"配体 {did} 创建MoleculeSplitter失败: {e}")
                return False, did, error_msg, []
            
            try:
                fragments_data = []
                for fragment in splitter.fragments:
                    fragment_info = {
                        'source_did': did,
                        'fragment_smiles': fragment['fragment_smiles'],
                        'fragment_DID': fragment['fragment_DID'],
                        'fragment_IRL_did': fragment['fragment_IRL_did'],
                        'fragment_IRL_smiles': fragment['fragment_IRL_smiles'],
                    }
                    fragments_data.append(fragment_info)
                
                return True, did, f"处理成功，生成 {len(fragments_data)} 个片段", fragments_data
                
            except Exception as e:
                error_msg = f"片段生成失败: {e}"
                # logger.error(f"配体 {did} 片段生成失败: {e}")
                # logger.error(f"错误详情: {traceback.format_exc()}")
                return False, did, error_msg, []
            
        except Exception as e:
            error_msg = f"L4片段生成失败: {e}"
            logger.error(f"配体 {did} L4片段生成失败: {e}")
            return False, did, error_msg, []
        
    except Exception as e:
        error_msg = f"处理异常: {e}"
        logger.error(f"配体 {ligand_data[0] if ligand_data else '未知'} 处理异常: {e}")
        logger.error(f"错误详情: {traceback.format_exc()}")
        return False, str(ligand_data[0]) if ligand_data else "未知", error_msg, []

def worker_process(task_queue: Queue, result_queue: Queue, worker_id: int, 
                  marker: MoleculeMarker, irl_list: List[Dict], ga_data_dict: Dict, 
                  batch_size: int = 100):
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
            
            success, did, message, fragments = process_single_ligand(task, marker, irl_list, ga_data_dict)
            
            batch_results.append((success, did, message, fragments))
            
            processed_count += 1
            if not success:
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
            logger.error(f"Worker {worker_id} 错误: {e}")
            failed_count += 1
            continue

class Step6_7Processor:
    def __init__(
        self,
        num_workers: int = None,
        record_limit: int = DEFAULT_RECORD_LIMIT,
        batch_size: int = None,
        max_queue_size: int = None,
        input_dir: str = None,
        output_dir: str = None,
    ):
        self.input_dir = Path(input_dir or INPUT_DIR)
        self.output_dir = Path(output_dir or OUTPUT_DIR)
        self.output_dir.mkdir(exist_ok=True)
        self.ga_path = self.input_dir / "GA_with_id.csv"
        self.irl_path = self.input_dir / "IRL_filtered.csv"
        
        self.stats_output = self.output_dir / "step6_7_stats.json"
        self.fragments_output = self.output_dir / "fragments.csv"
        
        if num_workers is None:
            self.num_workers = get_optimal_worker_count()
        else:
            self.num_workers = min(num_workers, 64)
        
        self.record_limit = record_limit
        
        if batch_size is None:
            self.batch_size = 50
        else:
            self.batch_size = batch_size
        
        if max_queue_size is None:
            self.max_queue_size = self.num_workers * 5000
        else:
            self.max_queue_size = max_queue_size
        
        self.stats = {
            'start_time': time.time(),
            'total_ligands': 0,
            'processed_ligands': 0,
            'failed_ligands': 0,
            'processing_time': 0,
            'num_workers': self.num_workers,
            'batch_size': self.batch_size,
            'max_queue_size': self.max_queue_size,
            'record_limit': self.record_limit
        }
        
        self.ga_data = self._load_ga_data()
        self.irl_list = self._load_irl_data()
        
        self.marker = MoleculeMarker(self.ga_data, {irl['complex_smiles']: irl['DID'] for irl in self.irl_list})
        
        logger.info(f"Step6+7合并处理器初始化完成，输出目录: {self.output_dir}")
        logger.info(f"并行处理工作进程数: {self.num_workers}")
        logger.info(f"批处理大小: {self.batch_size}")
        logger.info(f"最大队列大小: {self.max_queue_size}")
        logger.info(f"记录限制: {self.record_limit if self.record_limit > 0 else '无限制'}")
        logger.info(f"加载了 {len(self.ga_data)} 个GA模式和 {len(self.irl_list)} 个IRL模式")
    
    def _load_ga_data(self) -> Dict[str, str]:
        ga_path = self.ga_path
        if not ga_path.exists():
            logger.warning(f"GA数据文件不存在: {ga_path}")
            return {}
        
        try:
            df = pd.read_csv(ga_path)
            logger.info(f"GA数据列名: {list(df.columns)}")
            
            if 'GA_SMILES' in df.columns and 'GA_ID' in df.columns:
                ga_data = {row['GA_SMILES']: row['GA_ID'] for _, row in df.iterrows()}
            elif 'SMILES' in df.columns and 'ID' in df.columns:
                ga_data = {row['SMILES']: row['ID'] for _, row in df.iterrows()}
            elif 'smiles' in df.columns and 'id' in df.columns:
                ga_data = {row['smiles']: row['id'] for _, row in df.iterrows()}
            else:
                logger.warning(f"未找到预期的列名，使用前两列作为SMILES和ID")
                df = df.iloc[:, :2]
                df.columns = ['GA_SMILES', 'GA_ID']
                ga_data = {row['GA_SMILES']: row['GA_ID'] for _, row in df.iterrows()}
            
            logger.info(f"成功加载 {len(ga_data)} 个GA模式")
            return ga_data
        except Exception as e:
            logger.error(f"加载GA数据失败: {e}")
            logger.error(f"错误详情: {traceback.format_exc()}")
            return {}
    
    def _load_irl_data(self) -> List[Dict]:
        irl_path = self.irl_path
        if not irl_path.exists():
            logger.warning(f"IRL数据文件不存在: {irl_path}")
            return []
        
        try:
            df = pd.read_csv(irl_path)
            logger.info(f"IRL数据列名: {list(df.columns)}")
            
            if 'complex_smiles' in df.columns and 'DID' in df.columns:
                irl_list = df.to_dict('records')
            elif 'SMILES' in df.columns and 'DID' in df.columns:
                df = df.rename(columns={'SMILES': 'complex_smiles'})
                irl_list = df.to_dict('records')
            elif 'smiles' in df.columns and 'did' in df.columns:
                df = df.rename(columns={'smiles': 'complex_smiles', 'did': 'DID'})
                irl_list = df.to_dict('records')
            else:
                logger.warning(f"未找到预期的列名，使用前两列作为SMILES和DID")
                df = df.iloc[:, :2]
                df.columns = ['complex_smiles', 'DID']
                irl_list = df.to_dict('records')
            
            logger.info(f"成功加载 {len(irl_list)} 个IRL模式")
            return irl_list
        except Exception as e:
            logger.error(f"加载IRL数据失败: {e}")
            logger.error(f"错误详情: {traceback.format_exc()}")
            return []
    
    def _get_total_ligand_count(self) -> int:
        try:
            ligand_path = self.input_dir / "repaired_ligand_data.csv"
            if not ligand_path.exists():
                logger.error(f"配体数据文件不存在: {ligand_path}")
                return 0

            df = pd.read_csv(ligand_path, low_memory=False)
            logger.info(f"配体数据列名: {list(df.columns)}")

            if 'ligand_new_did' in df.columns and 'new_smiles' in df.columns:
                pass
            elif 'DID' in df.columns and 'complex_smiles' in df.columns:
                df = df.rename(columns={'DID': 'ligand_new_did', 'complex_smiles': 'new_smiles'})
            elif 'did' in df.columns and 'smiles' in df.columns:
                df = df.rename(columns={'did': 'ligand_new_did', 'smiles': 'new_smiles'})
            else:
                logger.warning(f"未找到预期的列名，使用前两列作为DID和SMILES")
                df = df.iloc[:, :2]
                df.columns = ['ligand_new_did', 'new_smiles']

            total_count = len(df)
            logger.info(f"成功加载 {total_count} 个配体记录")
            return total_count

        except Exception as e:
            logger.error(f"获取配体数量失败: {e}")
            logger.error(f"错误详情: {traceback.format_exc()}")
            return 0
    
    def _producer_process(self, task_queue: Queue, num_workers: int):
        try:
            ligand_path = self.input_dir / "repaired_ligand_data.csv"
            if not ligand_path.exists():
                logger.error(f"配体数据文件不存在: {ligand_path}")
                return

            df = pd.read_csv(ligand_path, low_memory=False)
            logger.info(f"配体数据列名: {list(df.columns)}")

            if 'ligand_new_did' in df.columns and 'new_smiles' in df.columns:
                pass
            elif 'DID' in df.columns and 'complex_smiles' in df.columns:
                df = df.rename(columns={'DID': 'ligand_new_did', 'complex_smiles': 'new_smiles'})
            elif 'did' in df.columns and 'smiles' in df.columns:
                df = df.rename(columns={'did': 'ligand_new_did', 'smiles': 'new_smiles'})
            else:
                df = df.iloc[:, :2]
                df.columns = ['ligand_new_did', 'new_smiles']

            total_processed = 0
            batch_size = 10000
            offset = 0

            while True:
                df_batch = df.iloc[offset:offset + batch_size]
                if df_batch.empty:
                    break

                for _, row in df_batch.iterrows():
                    if self.record_limit > 0 and total_processed >= self.record_limit:
                        break
                        
                    while task_queue.qsize() >= 4000:
                        time.sleep(0.01)
                    
                    task_queue.put((row['ligand_new_did'], row['new_smiles']))
                    total_processed += 1
                    
                    if self.record_limit > 0 and total_processed >= self.record_limit:
                        break
                
                if self.record_limit > 0 and total_processed >= self.record_limit:
                    break
                    
                offset += batch_size

            logger.info(f"生产者完成，总共处理了 {total_processed} 个配体")

        except Exception as e:
            logger.error(f"生产者进程错误: {e}")
            logger.error(f"错误详情: {traceback.format_exc()}")
        
        finally:
            for _ in range(num_workers):
                task_queue.put(None)
    
    def _save_stats(self):
        self.stats['processing_time'] = time.time() - self.stats['start_time']
        
        stats_data = {
            'timestamp': datetime.now().isoformat(),
            'processing_time_seconds': self.stats['processing_time'],
            'statistics': self.stats,
            'output_files': {
                'stats': str(self.stats_output)
            }
        }
        
        with open(self.stats_output, 'w', encoding='utf-8') as f:
            json.dump(stats_data, f, indent=2, ensure_ascii=False)
        
        logger.info(f"统计信息已保存到: {self.stats_output}")
        logger.info(f"处理总结: 总计 {self.stats['total_ligands']} 个配体，成功 {self.stats['processed_ligands']} 个，失败 {self.stats['failed_ligands']} 个")
    
    def process_all_ligands(self) -> Dict:
        start_time = time.time()
        logger.info("开始处理所有配体...")
        
        total_count = self._get_total_ligand_count()
        if total_count == 0:
            logger.warning("没有找到配体数据，跳过处理")
            return self.stats
        
        self.stats['total_ligands'] = total_count
        logger.info(f"找到 {total_count} 个配体记录")
        
        task_queue = Queue(maxsize=self.max_queue_size)
        result_queue = Queue(maxsize=self.max_queue_size)
        
        logger.info(f"启动 {self.num_workers} 个工作进程...")
        processes = []
        for i in range(self.num_workers):
            p = Process(target=worker_process, args=(
                task_queue, result_queue, i, self.marker, 
                self.irl_list, self.ga_data, self.batch_size
            ))
            p.start()
            processes.append(p)
        
        logger.info("开始流式加载和处理数据...")
        
        producer_process = Process(target=self._producer_process, args=(task_queue, self.num_workers))
        producer_process.start()
        
        processed_count = 0
        failed_count = 0
        all_fragments = []
        
        with tqdm(total=total_count, desc="处理配体") as pbar:
            while processed_count < total_count:
                try:
                    result = result_queue.get(timeout=1)
                    if result is None:
                        continue
                    
                    success, did, message, fragments = result
                    processed_count += 1
                    pbar.update(1)
                    
                    if success:
                        self.stats['processed_ligands'] += 1
                        all_fragments.extend(fragments)
                    else:
                        self.stats['failed_ligands'] += 1
                        failed_count += 1
                    
                    pbar.set_postfix({
                        '成功': self.stats['processed_ligands'],
                        '失败': self.stats['failed_ligands'],
                        '片段数': len(all_fragments)
                    })
                    
                except Empty:
                    if all(not p.is_alive() for p in processes) and not producer_process.is_alive():
                        break
                    continue
        
        producer_process.join()
        for p in processes:
            p.join()
        
        if all_fragments:
            logger.info(f"保存 {len(all_fragments)} 个片段到 {self.fragments_output}")
            df = pd.DataFrame(all_fragments)
            df.to_csv(self.fragments_output, index=False, encoding='utf-8')
            logger.info(f"片段数据已保存到: {self.fragments_output}")
        else:
            logger.warning("没有生成任何片段")
        
        self._save_stats()
        
        processing_time = time.time() - start_time
        logger.info(f"配体标记和L4片段生成完成，耗时: {processing_time:.2f} 秒")
        logger.info(f"总共生成 {len(all_fragments)} 个片段")
        
        return self.stats

def main():
    try:
        import argparse
        import sys

        if len(sys.argv) > 1 and not any(a.startswith("-") for a in sys.argv[1:]):
            num_workers = None
            record_limit = DEFAULT_RECORD_LIMIT
            batch_size = None
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
            if len(sys.argv) > 3:
                try:
                    batch_size = int(sys.argv[3])
                    print(f"使用指定的批处理大小: {batch_size}")
                except ValueError:
                    print("无效的批处理大小，使用默认值")
            processor = Step6_7Processor(
                num_workers=num_workers,
                record_limit=record_limit,
                batch_size=batch_size,
            )
        else:
            parser = argparse.ArgumentParser(description="ChemDB Step6+7")
            parser.add_argument("--input-dir", type=str, default=INPUT_DIR)
            parser.add_argument("--output-dir", type=str, default=OUTPUT_DIR)
            parser.add_argument("--workers", type=int, default=None)
            parser.add_argument("--limit", type=int, default=DEFAULT_RECORD_LIMIT)
            parser.add_argument("--batch-size", type=int, default=None)
            args = parser.parse_args()
            processor = Step6_7Processor(
                num_workers=args.workers,
                record_limit=args.limit,
                batch_size=args.batch_size,
                input_dir=args.input_dir,
                output_dir=args.output_dir,
            )

        logger.info("启动Step6+7合并处理器...")
        
        logger.info("开始处理配体数据...")
        stats = processor.process_all_ligands()
        
        print(f"\n{'='*50}")
        print(f"处理完成！")
        print(f"{'='*50}")
        print(f"总配体数: {stats['total_ligands']}")
        print(f"成功处理: {stats['processed_ligands']}")
        print(f"处理失败: {stats['failed_ligands']}")
        print(f"处理时间: {stats['processing_time']:.2f} 秒")
        print(f"成功率: {(stats['processed_ligands']/stats['total_ligands']*100):.1f}%" if stats['total_ligands'] > 0 else "成功率: N/A")
        print(f"{'='*50}")
        
        return 0
        
    except KeyboardInterrupt:
        logger.info("用户中断处理")
        print("\n处理被用户中断")
        return 1
        
    except Exception as e:
        logger.error(f"处理过程中发生错误: {e}")
        logger.error(f"错误详情: {traceback.format_exc()}")
        print(f"\n处理过程中发生错误: {e}")
        return 1

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code) 