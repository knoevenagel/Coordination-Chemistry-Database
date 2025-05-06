import os
import time
import json
import warnings
import pandas as pd
from tqdm import tqdm
from rdkit import RDLogger
from database import ComplexDatabase
from split_complex import ChemicalComplex
from insert_new_attribute import add_new_attribute_to_table
from repair_ligand import RepairLigandObject
from utils import generate_element_symbol_with_did, insert_data_to_table, generate_info_dict_list, load_elements_list, generate_ligand_did_list_from_database
from concurrent.futures import ProcessPoolExecutor, as_completed
from mysql.connector.pooling import MySQLConnectionPool

warnings.filterwarnings("ignore", category=UserWarning, module='rdkit')
RDLogger.DisableLog('rdApp.warning')

# 全局变量：子进程中的连接池实例
conn_pool = None

def init_process_pool(db_config):
    """
    在子进程中初始化连接池，并对 ComplexDatabase 类进行 monkey-patch 使其使用连接池。
    """
    global conn_pool
    conn_pool = MySQLConnectionPool(
        pool_name="mypool",
        pool_size=3,  # 可根据需求调整连接池大小
        **db_config
    )

    # 对 ComplexDatabase 进行 monkey-patch，使其使用连接池中的连接
    original_init = ComplexDatabase.__init__
    def new_init(self, config, table_name):
        # 先调用原有的初始化逻辑，这样ComplexDatabase原有的属性（如db_host）会被正常初始化
        original_init(self, config, table_name)
        # 然后再使用连接池获取连接，覆盖原有连接
        self.conn = conn_pool.get_connection()

    ComplexDatabase.__init__ = new_init


def process_one_dict(dict_info, metal_list, table_name, attribute_name, db_config):
    """
    处理单个 dict_info 对象。
    包括创建 ChemicalComplex 实例、合并邻居信息、插入数据到数据库。
    """
    complex_db = ComplexDatabase(db_config, table_name)
    complex_instance = ChemicalComplex(dict_info, metal_list)
    origin_dict = complex_instance.merge_neighbor_info()
    if origin_dict is None:
        smiles = dict_info['smiles']
        print(f"Skipping SMILES '{smiles}' due to missing data.")
        return None

    # 插入复合物信息和 replicate 分子信息
    complex_db.insert_complex_info(origin_dict)
    complex_db.insert_rep_mol(origin_dict)
    return None


def main(metal_file_path, p_elements_path, cid_file_name, table_name, attribute_name, db_config):
    start_time = time.time()

    # 创建数据库实例（此处不会真正创建连接，因为子进程中才初始化连接池）
    complex_db = ComplexDatabase(db_config, table_name)
    # 插入新属性到表中
    add_new_attribute_to_table(table_name, attribute_name, db=complex_db)

    # 加载金属和 p 元素列表
    metal_list = load_elements_list(metal_file_path)
    p_elements = load_elements_list(p_elements_path)

    # 生成需要处理的 dict_list
    dict_list = generate_info_dict_list(cid_file_name, metal_list, p_elements)

    # 插入元素信息到表中（一次性操作）
    for atomic_number in range(1, 119):
        DID, _ = generate_element_symbol_with_did(atomic_number)
        insert_data_to_table(atomic_number, DID, table_name=table_name, db=complex_db)

    # 动态调整进程数：至少8个，不超过 CPU 核数和任务数
    # cpu_count = os.cpu_count() or 8
    # if len(dict_list) > 0:
    #     max_workers = min(cpu_count, len(dict_list))
    #     max_workers = max(max_workers, 4)  # 至少使用8个进程
    # else:
    #     max_workers = 1
    max_workers = 4

    futures = []
    with ProcessPoolExecutor(
        max_workers=max_workers,
        initializer=init_process_pool,
        initargs=(db_config,)
    ) as executor:
        for dict_info in dict_list:
            futures.append(
                executor.submit(process_one_dict, dict_info, metal_list, table_name, attribute_name, db_config)
            )

        # 使用 tqdm 显示进度，desc中显示当前处理的文件名
        for _ in tqdm(as_completed(futures), total=len(futures), desc=f"Processing {os.path.basename(cid_file_name)}"):
            pass

    end_time = time.time()
    processing_time = end_time - start_time
    line_count = len(dict_list)
    return line_count, processing_time


if __name__ == "__main__":
    with open('./config/config.json', 'r') as f:
        config = json.load(f)

    mysql_config = config['database']
    mysql_config_path = "./config/config.json"
    metal_path = './data/metal_list.txt'
    table_name = "Database80"
    attribute_name = "level"
    p_block_elements_path = './data/p_elements_list.txt'

    # 批处理指定文件夹下所有 CSV 文件
    cid_file_path = "/Users/weilinzou/Code/DatabaseProject/project_data/pubchem_raw_data_processed"
    file_name_list = [f for f in os.listdir(cid_file_path) if f.endswith('.csv')]


    summary_dir = "./log_file"
    summary_file = os.path.join(summary_dir, "csv_processing_summary.csv")

    # 如果目录不存在则创建
    if not os.path.exists(summary_dir):
        os.makedirs(summary_dir)

    with open(summary_file, 'w', encoding='utf-8') as sf:
        sf.write("csv_name,line_count,processing_time\n")

    for file_name in file_name_list:
        cid_file_name = os.path.join(cid_file_path, file_name)
        line_count, processing_time = main(metal_path, p_block_elements_path, cid_file_name, table_name, attribute_name, mysql_config)
        # 将统计信息写入 summary 文件
        with open(summary_file, 'a', encoding='utf-8') as sf:
            sf.write(f"{file_name},{line_count},{processing_time}\n")

    # 可选：处理结束后统计数据
    result = generate_ligand_did_list_from_database(mysql_config, table_name)
    print(len(result))