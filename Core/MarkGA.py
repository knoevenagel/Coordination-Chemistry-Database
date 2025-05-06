import pandas as pd
from rdkit import Chem
import json
import mysql.connector
from tqdm import tqdm  # 加入 tqdm 用于显示进度条

class DatabaseHandler:
    def __init__(self, config_file):
        self.config = self.load_config(config_file)
        self.conn = mysql.connector.connect(
            host=self.config['host'],
            user=self.config['user'],
            password=self.config['password'],
            database=self.config['database']
        )
        self.cursor = self.conn.cursor()

    def load_config(self, config_file):
        with open(config_file, 'r') as file:
            return json.load(file)['database']

    def create_table(self, table_name):
        create_table_query = f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            DID VARCHAR(255) PRIMARY KEY,
            molecule_smiles TEXT,
            molecule_data JSON
        )
        """
        self.cursor.execute(create_table_query)
        self.conn.commit()

    def insert_molecule_data(self, table_name, did, smiles, marked_data):
        insert_query = f"""
        INSERT INTO {table_name} (DID, molecule_smiles, molecule_data)
        VALUES (%s, %s, %s)
        ON DUPLICATE KEY UPDATE molecule_smiles = VALUES(molecule_smiles), molecule_data = VALUES(molecule_data)
        """
        self.cursor.execute(insert_query, (did, smiles, marked_data))
        self.conn.commit()

    def close(self):
        self.cursor.close()
        self.conn.close()


class MoleculeMarker:
    def __init__(self, ga_data, irl_data):
        self.ga_data = ga_data
        self.irl_data = irl_data

    def mark_smiles(self, smiles):
        input_mol = Chem.MolFromSmiles(smiles)
        if input_mol is None:
            return Chem.MolFromSmiles('')
        for atom in input_mol.GetAtoms():
            atom.SetProp("atomNote", str(atom.GetIdx()))
        self.mark_IRL(input_mol)
        self.mark_GA(input_mol)
        return input_mol

    def mark_IRL(self, molecule):
        for irl_smiles, irl_id in self.irl_data.items():
            irl_mol = Chem.MolFromSmiles(irl_smiles)
            if irl_mol:
                matching_positions = molecule.GetSubstructMatches(irl_mol)
                for match in matching_positions:
                    self._add_label_to_matched_atoms(molecule, irl_id.upper(), "IRL", irl_mol, match)
        return molecule

    def mark_GA(self, molecule):
        for ga_smiles, ga_id in self.ga_data.items():
            ga_mol = Chem.MolFromSmiles(ga_smiles)
            if ga_mol:
                matching_positions = molecule.GetSubstructMatches(ga_mol)
                for match in matching_positions:
                    self._add_label_to_matched_atoms(molecule, ga_id, "GA", ga_mol, match)
        return molecule

    def _add_label_to_matched_atoms(self, molecule, atom_id, label_type, substruct_mol, match):
        for atom_idx in match:
            atom = molecule.GetAtomWithIdx(atom_idx)
            current_labels = atom.GetPropsAsDict().get(f"{label_type}_ID", "").split(";")
            if atom_id not in current_labels:
                current_labels.append(atom_id)
            atom.SetProp(f"{label_type}_ID", ";".join(current_labels))

    def get_marked_data(self, molecule):
        marked_data = {}
        for atom in molecule.GetAtoms():
            atom_idx = atom.GetIdx()
            atom_symbol = atom.GetSymbol()
            irl_ids = [irl_id for irl_id in atom.GetProp("IRL_ID").split(";") if irl_id] if atom.HasProp("IRL_ID") else []
            ga_ids = [ga_id for ga_id in atom.GetProp("GA_ID").split(";") if ga_id] if atom.HasProp("GA_ID") else []
            marked_data[f"atom_{atom_idx}"] = {
                "ga_ids": ga_ids,
                "IRL_ids": irl_ids,
                "atom_symbol": atom_symbol
            }
        return marked_data


def read_csv_to_dict(csv_path, id_column, smiles_column):
    df = pd.read_csv(csv_path)
    data_dict = {row[smiles_column]: row[id_column] for _, row in df.iterrows()}
    return data_dict


def process_molecules(config_file, ga_csv_path, irl_csv_path, input_csv_path, database_table_name):
    # 读取GA和IRL数据
    ga_data = read_csv_to_dict(ga_csv_path, 'GA_ID', 'GA_SMILES')
    irl_data = read_csv_to_dict(irl_csv_path, 'DID', 'complex_smiles')

    # 初始化标记类和数据库处理类
    molecule_marker = MoleculeMarker(ga_data, irl_data)
    db_handler = DatabaseHandler(config_file)

    # 创建数据库表
    db_handler.create_table(database_table_name)

    # 读取分子数据
    df = pd.read_csv(input_csv_path)

    # 使用 tqdm 显示进度条
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Processing molecules", unit="row"):
        did = row['DID']
        complex_smiles = row['complex_smiles']

        mol = Chem.MolFromSmiles(complex_smiles)
        if mol is None:
            tqdm.write(f"SMILES {complex_smiles} 无法解析，跳过该行")
            continue

        # 标记分子
        marked_molecule = molecule_marker.mark_smiles(complex_smiles)

        # 获取标记数据
        marked_data = molecule_marker.get_marked_data(marked_molecule)

        # 将标记数据转换为JSON格式
        marked_data_json = json.dumps(marked_data)

        # 将标记后的数据插入数据库
        db_handler.insert_molecule_data(database_table_name, did, complex_smiles, marked_data_json)

    # 关闭数据库连接
    db_handler.close()

    print("所有分子已成功标记并插入到数据库中。")


# 参数化调用
if __name__ == "__main__":
    config_file = './config/config.json'
    ga_csv_path = './GA_IRL_data/SemiGA_with_id.csv'
    irl_csv_path = './GA_IRL_data/IRL_sorted_PdZn_1234.csv'
    input_csv_path = '/Users/weilinzou/Code/DatabaseProject/project_data/all_ligand/ligands03_with_GAC.csv'
    database_table_name = 'ligand03GAC'

    process_molecules(config_file, ga_csv_path, irl_csv_path, input_csv_path, database_table_name)