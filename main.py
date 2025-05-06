import os
import yaml
from tqdm import tqdm
from datetime import datetime
from Utils.load_pubchem_data import ProcessPubchemRawData
from Utils.load_database_config import  load_mysql_config,load_neo4j_config
from Core.CC_split import ChemicalComplex,load_elements_list
from Core.neo4j_insert import Neo4jWriter
from Core.mysql_insert import MySQLLigandComplexDB

def load_yaml_config(yaml_config_path):
    with open(yaml_config_path, "r") as f:
        cfg = yaml.safe_load(f)

    return cfg

def generate_pubchem_table_name(prefix):
    today = datetime.now().strftime("%Y%m%d")
    table_name = f"{prefix}_{today}"
    return table_name

def process_one_complex(complex_info,metal_list):

    complex_instance = ChemicalComplex(complex_info,metal_list) 
    merged_dict = complex_instance.merge_neighbor_info()

    return merged_dict
    # print(merged_dict)

def process_one_batch(batch_df, yaml_config, show_progress=False):
    """
    处理单个 batch 的 coordination complex 数据。
    返回一个包含结构化信息的字典列表。
    """
    metal_list_path = yaml_config['paths']['metal_list_path']
    metal_list = load_elements_list(metal_list_path)

    info_dicts = generate_info_dict_list(batch_df)

    results = []
    iterator = tqdm(info_dicts) if show_progress else info_dicts

    for info_dict in iterator:
        try:
            result = process_one_complex(info_dict, metal_list)
            if not isinstance(result, dict):
                raise ValueError(f"Expected dict return, got {type(result)}")
            results.append(result)
        except Exception as e:
            print(f"[ERROR] Failed processing CID {info_dict['cid']}: {e}")
            continue

    return results

def process_all_batches_and_insert(df, yaml_config, batch_size=100):
    """
    主控函数：将整个 DataFrame 分批处理并插入数据库。

    Args:
        df (pd.DataFrame): 包含所有 coordination complex 数据。
        yaml_config (dict): 配置文件。
        batch_size (int): 每个批次的大小。
    """
    total = len(df)
    print(f"Total entries to process: {total}")
    


    for i in range(0, total, batch_size):
        print(f"\nProcessing batch {i // batch_size + 1} ({i} to {min(i+batch_size, total)})...")
        batch_df = df.iloc[i:i+batch_size]

        batch_results = process_one_batch(batch_df, yaml_config, show_progress=True)

        # 在这里对处理后的 batch_results 做插入数据库操作
        insert_batch_results_to_database(batch_results,yaml_config)

def process_all_and_insert_sequentially(df, yaml_config):
    """
    顺序处理整个 DataFrame，每次只处理一条记录并插入数据库。
    """
    metal_list_path = yaml_config['paths']['metal_list_path']
    metal_list = load_elements_list(metal_list_path)
    mysql_processor = MySQLLigandComplexDB(yaml_config)
    # 使用你的类，通过 yaml_config 初始化 Neo4jWriter
    neo4j_writer = Neo4jWriter(yaml_config)
    info_dicts = generate_info_dict_list(df)

    for info_dict in tqdm(info_dicts):
        try:
            result = process_one_complex(info_dict, metal_list)
            if isinstance(result, dict):
                mysql_processor.insert_simple_complex_info(result)
                neo4j_writer.insert_complex_and_ligands(result)
        except Exception as e:
            print(f"[ERROR] Failed processing CID {info_dict['cid']}: {e}")


def insert_batch_results_to_database(batch_results,yaml_config):
    """
    将一个 batch 的处理结果插入数据库。
    """
    mysql_processor = MySQLLigandComplexDB(yaml_config)  

    for result in batch_results:
        try:
            # 自行替换为你的具体插入逻辑
            # print(result)
            mysql_processor.insert_simple_complex_info(result)
        except Exception as e:
            print(f"[DB ERROR] Failed to insert result with CID {result.get('cid')}: {e}")

def generate_info_dict_list(df):
    """
    将 DataFrame 中的每一行转换为标准的 info_dict 结构。

    Args:
        df (pd.DataFrame): 包含 'cid' 和 'isosmiles' 列的一个 batch。

    Returns:
        List[Dict]: 每行转换后的字典结构。
    """
    info_list = []

    for _, row in df.iterrows():
        cid = str(row['cid'])
        smiles = row['isosmiles']

        info_dict = {
            'cid': cid,
            'smiles': smiles,
            'corcomp': True  # 可自定义判断
        }
        info_list.append(info_dict)

    return info_list

def main(yaml_config,generate_new_pubchemdata_table = False):

    # 1. load parameters
    mysql_config_path = yaml_config['database']['mysql_config_path']
    mysql_config = load_mysql_config(mysql_config_path)
    pubchem_data_path = yaml_config['paths']['pubchem_data_path']
    pubchem_table_sample_path = yaml_config['paths']['pubchem_table_sample']
    pubchem_table_name = yaml_config['database']['pubchem_table_name']
    pubchem_table_name_with_date = generate_pubchem_table_name(pubchem_table_name)

    # 2. initialize pubchemdata class

    PubChemProcessor = ProcessPubchemRawData(mysql_config)
    # generate tabel and insert pubchem raw data to mysql table
    if generate_new_pubchemdata_table == False:
        valid_pubchem_tabel_name = pubchem_table_name
    else:
        valid_pubchem_tabel_name = pubchem_table_name_with_date
        PubChemProcessor.create_pubchem_data_table(pubchem_table_sample_path,valid_pubchem_tabel_name)
        PubChemProcessor.process_folder(pubchem_data_path,valid_pubchem_tabel_name)

    # 3. get cid and smiles DataFrame
    df = PubChemProcessor.fetch_data_parallel_from_mysql(valid_pubchem_tabel_name)

    # process_all_batches_and_insert(df,yaml_config)
    process_all_and_insert_sequentially(df,yaml_config)

    
    
if __name__ == "__main__":
    yaml_config_path = '/Users/weilinzou/Code/DatabaseProject/ChemDB/Config/pipline_config.yaml'
    yaml_config = load_yaml_config(yaml_config_path)
    main(yaml_config)
    # test_info_dict = {'cid': '111111', 'smiles': 'C1=CC=C2C(=C1)C=CC=C2O.C1=CC=C2C(=C1)C=CC=C2O.[Mn]','corcomp':True }
    # # process_one_complex(test_info_dict,yaml_config)
    


    
    