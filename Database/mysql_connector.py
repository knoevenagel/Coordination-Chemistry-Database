import json
import mysql.connector
from rdkit.Chem import GetPeriodicTable

class ComplexDatabase:
    def __init__(self, db_config, table_name):
        self.db_host = db_config['host']
        self.db_user = db_config['user']
        self.db_password = db_config['password']
        self.db_name = db_config['database']
        self.table_name = table_name
        self.create_table()

    def create_table(self):
        try:    
            conn = mysql.connector.connect(host=self.db_host, user=self.db_user, password=self.db_password, database=self.db_name)
            cursor = conn.cursor()
            cursor.execute(f'''
                CREATE TABLE IF NOT EXISTS {self.table_name} (
                    DID VARCHAR(255) PRIMARY KEY,
                    complex_smiles TEXT,
                    metal_info TEXT,
                    FromPubChem BOOLEAN,
                    rep_mol TEXT,
                    CorComp BOOLEAN,    
                    ligand INT,  -- 将ligand的类型改为整数
                    source_mol LONGTEXT,
                    label_list LONGTEXT,
                    inactive INT -- 添加新的整数列
                );

            ''')    
            conn.commit()
            conn.close()
            # print("Table created successfully.")
        except Exception as e:
            print(f"Error creating table: {e}")

    def insert_complex_info(self, complex_info):
        try:
            conn = mysql.connector.connect(host=self.db_host, user=self.db_user, password=self.db_password, database=self.db_name)
            cursor = conn.cursor()

            complex_did = complex_info['complex_info'].get('DID', None)
            if complex_did is None:
                raise ValueError("DID不存在")

            complex_smiles = complex_info['complex_info']['complex_smiles']
            #metal_info = json.dumps([{metal['central_metal']: metal.get('valence', None)} for metal in complex_info['central_metal_info']])
            inactive = complex_info['complex_info']['inactive']
            metal_info_elements = set()

            for metal in complex_info['central_metal_info']:
                element_symbol = metal['central_metal']
                valence = metal.get('valence', None)
                valence_str = f"'{valence}'" if valence is not None else 'NULL'
                metal_element = f"{{'{element_symbol}': {valence_str}}}"
                metal_info_elements.add(metal_element)

            metal_info_json = '[' + ','.join(metal_info_elements) + ']'

            FromPubChem = True
            CorComp = True
            ligand = 0
            source_mol = None

            rep_mol = list(set([neighbor['ligand_did'] for neighbor in complex_info['neighbor_info']]))

            # 获取所有不同的 metal_did
            metal_did_set = set()
            for metal_atom in complex_info['central_metal_info']:
                metal_did = metal_atom.get('metal_did')
                if metal_did is not None:
                    if metal_did not in metal_did_set:
                        metal_did_set.add(metal_did)
                        rep_mol.append(metal_did)
            #print(metal_did_set)
            # 更新数据库中的 source_mol
            for metal_did in metal_did_set:
                cursor.execute(f'''
                    SELECT source_mol
                    FROM {self.table_name}
                    WHERE DID = %s
                ''', (metal_did,))
                existing_source_mol = cursor.fetchone()[0]  # 获取现有的 source_mol 值

                if existing_source_mol:
                    source_mol_list = existing_source_mol.split(',')
                    source_mol_list.append(complex_did)  # 将新的 complex_did 添加到列表中
                    source_mol_str = ','.join(source_mol_list)  # 将列表转换为字符串
                else:
                    source_mol_str = complex_did  # 如果原先的 source_mol 为空，则直接使用新的 complex_did

                cursor.execute(f'''
                    UPDATE {self.table_name}
                    SET source_mol = %s
                    WHERE DID = %s
                ''', (source_mol_str, metal_did))

            # 正常部分，不包括给金属添加 source_mol    
            cursor.execute(f'SELECT COUNT(*) FROM {self.table_name} WHERE DID = %s', (complex_did,))
            count = cursor.fetchone()[0]

            if count > 0:
                # cursor.execute(f'SELECT * FROM {self.table_name} WHERE DID = %s', (complex_did,))
                # existing_entry = cursor.fetchone()
                # existing_label_list = json.loads(existing_entry[-1]) if existing_entry[-1] else []
                # new_label_list = []

                # for neighbor in complex_info['neighbor_info']:
                #     new_label = [neighbor['coordinating_atom_index'], str(neighbor['bond_type']).split('.')[-1], neighbor['central_atom'][0], complex_did]
                #     if new_label not in existing_label_list:
                #         new_label_list.append(new_label)

                # merged_label_list = existing_label_list + new_label_list

                cursor.execute(f'''
                    UPDATE {self.table_name}
                    SET complex_smiles = %s, metal_info = %s, FromPubChem = %s, rep_mol = %s, CorComp = %s, ligand = %s, source_mol = %s, inactive = %s
                    WHERE DID = %s
                ''', (complex_smiles, metal_info_json, FromPubChem, json.dumps(rep_mol), CorComp, ligand, source_mol, complex_did,inactive))
            else:
                cursor.execute(f'''
                    INSERT INTO {self.table_name} (DID, complex_smiles, metal_info, FromPubChem, rep_mol, CorComp, ligand, source_mol,inactive)
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
                ''', (complex_did, complex_smiles, metal_info_json, FromPubChem, json.dumps(rep_mol), CorComp, ligand, source_mol,inactive))

            conn.commit()
            conn.close()

        except Exception as e:
            print(f"Error inserting complex information: {e}")

    def insert_rep_mol(self, complex_info):
        try:
            conn = mysql.connector.connect(host=self.db_host, user=self.db_user, password=self.db_password, database=self.db_name)
            cursor = conn.cursor()

            complex_did = complex_info['complex_info'].get('DID', None)
            neighbor_info = complex_info.get('neighbor_info', [])
            complex_inactive = complex_info['complex_info']['inactive']
            if complex_inactive != 0:
                return
            
            if neighbor_info:
                for neighbor in neighbor_info:
                    ligand_did = neighbor['ligand_did']
                    ligand_smiles = neighbor['ligand_smiles']
                    label_list = self.prepare_label_list(neighbor, complex_did, complex_info)
                    FromPubChem = 0
                    CorComp = 0
                    ligand = 1
                    inactive = 0

                    cursor.execute(f'SELECT * FROM {self.table_name} WHERE DID = %s', (ligand_did,))
                    existing_entry = cursor.fetchone()

                    if existing_entry:
                        # 注意，下面两个列表的索引是根据我们刚才定义的table的结构来的，如果table的结构发生变化，这里的索引也需要相应调整
                        existing_label_list_str = existing_entry[8]
                        existing_label_list = json.loads(existing_label_list_str) if existing_label_list_str else []

                        if not any(label == label_list for label in existing_label_list):
                            existing_label_list.extend(label_list)
                        existing_source_mol = existing_entry[7]
                        if existing_source_mol:
                            source_mol_list = existing_source_mol.split(',')
                            source_mol_list.append(complex_did)
                            source_mol_str = ','.join(source_mol_list)
                        cursor.execute(f'''
                            UPDATE {self.table_name}
                            SET complex_smiles = %s,
                                source_mol = %s,
                                label_list = %s,
                                FromPubChem = %s,
                                CorComp = %s,
                                ligand = %s
                            WHERE DID = %s
                        ''', (ligand_smiles, source_mol_str, json.dumps(existing_label_list), FromPubChem, CorComp, ligand, ligand_did))
                    else:
                        cursor.execute(f'''
                            INSERT INTO {self.table_name} (DID, complex_smiles, source_mol, label_list, FromPubChem, CorComp, ligand, inactive)
                            VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
                        ''', (ligand_did, ligand_smiles, complex_did, json.dumps(label_list), FromPubChem, CorComp, ligand, inactive))
                
                conn.commit()
                conn.close()
            else:
                print(f"No neighbor information available for complex with DID: {complex_did}")

        
        except Exception as e:
            print(f"Error inserting rep_mol information: {e}")

    def prepare_label_list(self, neighbor, complex_did, complex_info):
        label_list = []
        coordinating_atom_index = neighbor.get('coordinating_atom_index')
        central_atom_did = neighbor['central_atom_did']

        # 如果 central_atom_did 是字符串，则直接放入集合中
        if isinstance(central_atom_did, str):
            central_atom_did_set = {central_atom_did}
        # 如果 central_atom_did 是列表，则使用 set() 函数将其转换为集合，自动去除重复项
        elif isinstance(central_atom_did, list):
            central_atom_did_set = set(central_atom_did)
        # 转换为列表，确保顺序不变
        central_atom_did_unique = list(central_atom_did_set)
        # 是否配位键判断
        if coordinating_atom_index is not None:
            bond_type = str(neighbor['bond_type']).split('.')[-1]
            central_atom = neighbor['central_atom'][0]

            #print(central_atom)
            label_list.append([str(coordinating_atom_index), bond_type, central_atom_did_unique, str(complex_did)])
        else:
            bond_type = 'IONIC'
            coordinating_atom_index = None
            central_atom = neighbor['central_atom']

            #print(central_atom)
            label_list.append([str(coordinating_atom_index), bond_type, central_atom_did_unique, str(complex_did)])
        
        return label_list


if __name__ == "__main__":

    with open('./config/config.json', 'r') as f:
        config = json.load(f)

    mysql_config = config['database']

    # 实例化ComplexDatabase类
    db = ComplexDatabase(mysql_config, "test123")

    # 假设你有一个名为complex_info的字典，它包含要插入数据库的复杂信息
    # 这里是一个示例调用
    complex_info = {'complex_info': {'DID': 'D167530506', 'complex_smiles': 'C1=CC=C2C(=C1)C=CC=C2O.C1=CC=C2C(=C1)C=CC=C2O.[Mn].[Sn+2]'}, 'central_metal_info': [{'atom_index': 22, 'central_metal': 'Mn', 'valence': 0, 'metal_did': 'D25E'}, {'atom_index': 23, 'central_metal': 'Sn', 'valence': 0, 'metal_did': 'D50E'}], 'neighbor_info': [{'coordinating_atom_index': None, 'ligand_smiles': 'C1=CC=C2C(=C1)C=CC=C2O', 'bond_type': 'IONIC', 'central_atom': ['Mn', 'Sn'], 'central_atom_did': ['D25E', 'D50E'], 'ligand_did': 'D8808444227'}]}
    complex_info2 = {'complex_info': {'DID': 'D10051606', 'complex_smiles': 'CC1N(C=CN1CC2=CC=CC(=N2)CN3CN(C=C3)CCCO)CCCO.[OH-].[Ag+]'}, 'central_metal_info': [{'atom_index': 28, 'central_metal': 'Ag', 'valence': 0, 'metal_did': 'D47E'}], 'neighbor_info': [{'coordinating_atom_index': None, 'ligand_smiles': 'CC1N(C=CN1CC2=CC=CC(=N2)CN3CN(C=C3)CCCO)CCCO', 'bond_type': 'IONIC', 'central_atom': ['Ag'], 'central_atom_did': ['D47E'], 'ligand_did': 'D1434359447'}, {'coordinating_atom_index': None, 'ligand_smiles': '[OH-]', 'bond_type': 'IONIC', 'central_atom': ['Ag'], 'central_atom_did': ['D47E'], 'ligand_did': 'D1662042841'}]}
    # complex_info3 = 
    # 插入复杂信息
    db.insert_complex_info(complex_info)

    # 插入重复分子信息（如果有需要）
    # db.insert_rep_mol(complex_info)
