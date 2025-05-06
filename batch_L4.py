
import io
import csv
import sys
import random
import os
import json
import time
import mysql.connector
import matplotlib.pyplot as plt
from tqdm import tqdm
from rdkit import Chem
from PIL import Image
from rdkit import Chem
from collections import deque
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from collections import deque
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

sys.path.append('/Users/weilinzou/Code/DatabaseProject/ChemDB/Core')

import Core.L4_create as genL4

# 增加字段大小限制
csv.field_size_limit(sys.maxsize)  # 增加字段大小限制

def read_csv_to_dict(csv_file, id_column, smiles_column):
    """
    读取CSV文件，只提取所需的SMILES和ID列数据，并将其存储为字典。
    csv_file: 输入的CSV文件路径
    id_column: 包含ID的列名（如IRL的DID或GA的ID）
    smiles_column: 包含SMILES的列名
    return: 字典格式的数据，{SMILES: ID}
    """
    smiles_id_dict = {}

    with open(csv_file, mode='r', newline='', encoding='utf-8') as infile:
        csv_reader = csv.reader(infile)  # 使用csv.reader读取文件
        
        # 获取列名，确定SMILES和ID列的位置
        header = next(csv_reader)  # 读取第一行（列名）
        id_idx = header.index(id_column)  # 获取ID列索引
        smiles_idx = header.index(smiles_column)  # 获取SMILES列索引

        # 遍历每一行，读取需要的列数据
        for row in csv_reader:
            if len(row) > max(id_idx, smiles_idx):  # 确保行数据完整
                smiles = row[smiles_idx]  # 获取SMILES列的数据
                id_value = row[id_idx]  # 获取ID列的数据
                smiles_id_dict[smiles] = id_value  # 将SMILES作为键，ID作为值存储到字典中

    return smiles_id_dict

def read_csv_to_irl_list(csv_file, id_column, smiles_column):
    """
    读取CSV文件，将SMILES和ID列数据存储为字典列表，表头不包含在数据中。

    csv_file: 输入的CSV文件路径
    id_column: 包含ID的列名（如IRL的DID）
    smiles_column: 包含SMILES的列名
    return: 列表格式的数据，每个元素是一个字典，字典包含SMILES和ID的映射。
    """
    irl_list = []
    
    with open(csv_file, mode='r', newline='', encoding='utf-8') as infile:
        csv_reader = csv.reader(infile)  # 使用csv.reader读取文件
        
        # 获取列名，确定SMILES和ID列的位置
        header = next(csv_reader)  # 读取第一行（列名）
        id_idx = header.index(id_column)  # 获取ID列索引
        smiles_idx = header.index(smiles_column)  # 获取SMILES列索引

        # 遍历每一行，读取需要的列数据并将其存储为字典
        for row in csv_reader:
            if len(row) > max(id_idx, smiles_idx):  # 确保行数据完整
                irl_data = {
                    id_column: row[id_idx],  # 获取ID列的数据
                    smiles_column: row[smiles_idx]  # 获取SMILES列的数据
                }
                irl_list.append(irl_data)  # 将每行数据作为字典添加到列表中

    return irl_list

def batch_process_smiles(config_path, marked_db_table, output_path, irl_list, ga_list,fragment_table):
    """
    批量处理分子数据并保存分割后的图像。

    :param config_path: 数据库配置文件路径
    :param db_table: 数据库表名
    :param output_path: 图像保存路径
    :param irl_list: IRL 列表
    :param ga_list: GA 列表
    """
    import os
    import json
    import mysql.connector

    # 加载数据库配置
    db_config = load_db_config(config_path)

    # 创建输出目录
    os.makedirs(output_path, exist_ok=True)

    # 连接到数据库
    connection = mysql.connector.connect(
        host=db_config['host'],
        user=db_config['user'],
        password=db_config['password'],
        database=db_config['database']
    )
    cursor = connection.cursor()

    try:
        # 从数据库查询数据
        query = f"SELECT DID, molecule_smiles, molecule_data FROM {marked_db_table}"
        cursor.execute(query)
        rows = cursor.fetchall()

        # 使用 tqdm 添加进度条
        for index, row in tqdm(enumerate(rows), total=len(rows), desc="Processing molecules"):
            try:
                DID, smiles, marked_data = row
                # 如果 `marked_data` 是 JSON 字符串，则解析为字典
                if isinstance(marked_data, str):
                    marked_data = json.loads(marked_data)

                # 初始化 MoleculeSplitter
                splitter = genL4.MoleculeSplitter(smiles, marked_data, irl_list, ga_list)
                fragment_smiles = splitter.split_and_get_smiles()
                source_DID = DID
                genL4.insert_fragments_to_mysql(source_DID, fragment_smiles, config_path, fragment_table)

            except Exception as e:
                print(f"Error processing row {index}: {e}")

    finally:
        # 关闭数据库连接
        cursor.close()
        connection.close()


# 示例使用
if __name__ == "__main__":
        # 读取GA的CSV文件
    time1 = time.time()
    ga_csv_path = '/Users/weilinzou/Code/DatabaseProject/scidata/GA_IRL_data/SemiGA_with_id.csv'
    ga_id_column = 'GA_ID'  # GA的ID列名是GA_ID
    ga_smiles_column = 'GA_SMILES'  # GA的SMILES列名是Aromatic Rings SMILES
    ga_data_dict = read_csv_to_dict(ga_csv_path, ga_id_column, ga_smiles_column)
    IRL_csv_file = "/Users/weilinzou/Code/DatabaseProject/scidata/GA_IRL_data/IRL_sorted_PdZn_1234.csv"
    sorted_irl_list = read_csv_to_irl_list(IRL_csv_file, 'DID', 'complex_smiles')


    config_path = "./config/config.json"
    db_table = "ligand03GAC"
    fragment_db_table = 'ligandtotal'
    output_path = "/Users/weilinzou/Code/DatabaseProject/project_data/fragment_figures"

    batch_process_smiles(config_path, db_table, output_path, sorted_irl_list, ga_data_dict,fragment_db_table)
    
    time2 = time.time()
    print(time2-time1)
