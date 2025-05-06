import mysql.connector
from rdkit.Chem import GetPeriodicTable
import csv
import json
from rdkit import Chem
def load_elements_list(file_path)  :
    with open(file_path) as f:  
        lines = f.readlines()
    elements_list = list(map(lambda x: x.rstrip(),lines)) 
    return elements_list

def generate_element_symbol_with_did(atomic_number):
    """
    根据给定的原子序数生成元素符号和 DID。

    参数：
    atomic_number (int)：要生成元素符号和 DID 的原子序数。

    返回：
    DID (str)：生成的 DID。
    element_symbol (str)：生成的元素符号。
    """
    # 从原子序数获取元素符号
    pt = GetPeriodicTable()
    element_symbol = pt.GetElementSymbol(atomic_number)
    
    # 生成 DID
    DID = 'D' + str(atomic_number) + 'E'
    
    return DID, element_symbol

def add_did_info(input_dict, csv_filename):
    new_dict = input_dict.copy()  # 复制输入的字典以避免修改原始数据
    complex_smiles = input_dict['complex_info']['complex_smiles']
    try:
        # 尝试从 CSV 文件中获取对应的 DID
        with open(csv_filename, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if row['isosmiles'] == complex_smiles:
                    # 如果在 CSV 中找到了 DID，则直接使用该 DID
                    did = row.get('DID')
                    if did:
                        new_dict['complex_info']['DID'] = did
                        return new_dict  # 返回更新后的字典
                    else:
                        cid = row.get('cid')
                        if cid:
                            # 如果 CSV 中没有 DID，但有 CID，则创建一个类似于 "D" 开头的 DID
                            did = 'D' + cid
                            new_dict['complex_info']['DID'] = did
                            return new_dict  # 返回更新后的字典
                        else:
                            break  # 如果 CSV 中既没有 DID 也没有 CID，则继续下面的计算 DID 步骤
    
    except Exception as e:
        print(f"Error while searching for DID in CSV: {e}")

    # 如果以上步骤都没有找到 DID，则返回原始字典
    return new_dict
def insert_data_to_table(atomic_number, DID, db, table_name):
    try:
        # 连接到数据库
        conn = mysql.connector.connect(host=db.db_host, user=db.db_user, password=db.db_password, database=db.db_name)
        cursor = conn.cursor()

        # 准备 SQL 查询
        query = f"""
        INSERT INTO {table_name} (DID, complex_smiles, metal_info, FromPubChem, rep_mol, CorComp, ligand, source_mol, label_list, level)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """

        # 生成元素符号
        pt = GetPeriodicTable()
        element_symbol = pt.GetElementSymbol(atomic_number)

        # 准备数据
        metal_info = f"[{{'{element_symbol}': NULL}}]"
        complex_smiles = element_symbol
        FromPubChem = 1
        CorComp = 0
        ligand = 0
        source_mol = None
        label_list = None
        level = 0
        rep_mol = None

        # 执行查询
        data = (DID, complex_smiles, metal_info, FromPubChem, rep_mol, CorComp, ligand, source_mol, label_list, level)
        cursor.execute(query, data)

        # 提交更改并关闭连接
        conn.commit()
        conn.close()
        print(f"Data for element {atomic_number} inserted successfully.")

    except mysql.connector.Error as error:
        print("Error:", error)


def check_corcomp(input_smiles,metal_elements, p_block_elements):
    has_metal = any(element in input_smiles for element in metal_elements)
    has_p_block = any(element in input_smiles for element in p_block_elements)
    if has_metal and has_p_block:
        return True
    else:
        return False

def check_gibberish(input_smiles):
    mol = Chem.MolFromSmiles(input_smiles)
    if mol is None:
        return True
    else:
        return False


def generate_info_dict_list(csv_filename,metal_elements, p_block_elements):
    """
    从CSV文件中生成包含每行数据的字典列表。

    Args:
        csv_filename (str): CSV文件的路径。
        metal_elements(str): 金属中心文件的路径,用来确认是否是corcomp
        p_block_elements(str): p区非金属元素的文件路径
    Returns:
        list: 包含每行数据字典的列表,包含cid, smiles和corcomp
    """
    cid_smiles_list = []
    try:
        with open(csv_filename, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                cid = row.get('cid')
                smiles = row.get('isosmiles')
                # print(cid)
                if cid and smiles:
                    corcomp = check_corcomp(smiles,metal_elements, p_block_elements)

                    # gibberish = check_gibberish(smiles)
                    # if gibberish == True:
                    #     print(f"Skipping SMILES '{smiles}' due to invalid SMILES.")
                    cid_smiles_list.append({'cid': cid, 'smiles': smiles,'corcomp':corcomp}) 
    except Exception as e:
        print(f"Error while generating CID-SMILES list: {e}")
    return cid_smiles_list

def generate_ligand_did_list_from_database(db_config,table_name):
    try:
        # 连接到数据库
        conn = mysql.connector.connect(host=db_config['host'], user=db_config['user'], password=db_config['password'], database=db_config['database'])
        cursor = conn.cursor()

        # 准备 SQL 查询
        query = f"""
        SELECT DID
        FROM {table_name}
        WHERE ligand = 1
        AND inactive = 0
        AND CorComp = 0
        """

        # 执行查询
        cursor.execute(query)
        result = cursor.fetchall()

        # 关闭连接
        conn.close()

        return result

    except mysql.connector.Error as error:
        print("Error:", error)
        return None

if __name__ == "__main__":
    
    csv_file_path = '/home/weilin/ScientificProject/scientific/PubChemData/processed_data/Ac.csv'
    metal_elements_path = './data/metal_list.txt'
    p_block_elements_path = './data/p_elements_list.txt'
    db_config_path = "/home/weilin/ScientificProject/scientific/ligand_process/config/config.json"
    table_name = "test0916Ac"
    # metal_elements = load_elements_list(metal_elements_path)
    # p_elements = load_elements_list(p_block_elements_path)
    # cid_smiles_dict = generate_info_dict_list(csv_file_path,metal_elements, p_elements)
    # print(cid_smiles_dict)