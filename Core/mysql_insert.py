import mysql.connector
import json
import yaml

class MySQLLigandComplexDB:
    def __init__(self, yaml_config):

        mysql_config_path = yaml_config['database']['mysql_config_path']
        with open(mysql_config_path, 'r') as f:
            config = yaml.safe_load(f)

        db_config = config['database']


        self.db_host = db_config['host']
        self.db_user = db_config['user']
        self.db_password = db_config['password']
        self.db_name = db_config['database']

        self.table_name_CC = yaml_config['database']['coordinate_compound_table_name']
        self.table_name_ligand = yaml_config['database']['ligand_table_name']

        self.create_CC_table()
        self.create_ligand_table

    def create_CC_table(self):
        """
        创建配合物与配体共享的表结构。
        """
        try:
            conn = mysql.connector.connect(
                host=self.db_host,
                user=self.db_user,
                password=self.db_password,
                database=self.db_name
            )
            cursor = conn.cursor()

            cursor.execute(f'''
                CREATE TABLE IF NOT EXISTS {self.table_name_CC} (
                    DID VARCHAR(255) PRIMARY KEY,
                    complex_smiles TEXT,
                    metal_info TEXT,
                    FromPubChem BOOLEAN,
                    label_list LONGTEXT,
                    ligand INT,
                    inactive INT
                );
            ''')

            conn.commit()
            conn.close()
        except Exception as e:
            print(f"[ERROR] Failed to create table: {e}")

    def create_ligand_table(self):
        """
        创建配体表，仅包含 DID、ligand_smiles、ligand_inactive 三个字段。
        """
        try:
            conn = mysql.connector.connect(
                host=self.db_host,
                user=self.db_user,
                password=self.db_password,
                database=self.db_name
            )
            cursor = conn.cursor()

            cursor.execute(f'''
                CREATE TABLE IF NOT EXISTS {self.table_name_ligand} (
                    DID VARCHAR(255) PRIMARY KEY,
                    ligand_smiles TEXT,
                    ligand_inactive INT
                );
            ''')

            conn.commit()
            conn.close()
            print(f"[INFO] Ligand table '{self.table_name_ligand}' created successfully.")
        except Exception as e:
            print(f"[ERROR] Failed to create ligand table '{self.table_name_ligand}': {e}")

    def insert_simple_complex_info(self, complex_info):
        """
        插入配合物的简化信息。
        """
        try:
            conn = mysql.connector.connect(
                host=self.db_host,
                user=self.db_user,
                password=self.db_password,
                database=self.db_name
            )
            cursor = conn.cursor()

            did = complex_info['complex_info'].get('DID')
            smiles = complex_info['complex_info']['complex_smiles']
            inactive = complex_info['complex_info']['inactive']
            FromPubChem = True

            # 构建 metal_info
            metal_info_set = set()
            for metal in complex_info['central_metal_info']:
                symbol = metal['central_metal']
                valence = metal.get('valence')
                val_str = f"'{valence}'" if valence is not None else 'NULL'
                metal_info_set.add(f"{{'{symbol}': {val_str}}}")
            metal_info_json = '[' + ','.join(metal_info_set) + ']'

            # 插入或更新
            cursor.execute(f"SELECT COUNT(*) FROM {self.table_name_CC} WHERE DID = %s", (did,))
            if cursor.fetchone()[0]:
                cursor.execute(f'''
                    UPDATE {self.table_name_CC}
                    SET complex_smiles = %s, metal_info = %s, FromPubChem = %s, inactive = %s
                    WHERE DID = %s
                ''', (smiles, metal_info_json, FromPubChem, inactive, did))
            else:
                cursor.execute(f'''
                    INSERT INTO {self.table_name_CC} (DID, complex_smiles, metal_info, FromPubChem, inactive)
                    VALUES (%s, %s, %s, %s, %s)
                ''', (did, smiles, metal_info_json, FromPubChem, inactive))

            conn.commit()
            conn.close()
        except Exception as e:
            print(f"[ERROR] Failed to insert complex info: {e}")

    def insert_ligand_info(self, ligand_info):
        """
        插入配体信息，包括 label_list。
        """
        try:
            conn = mysql.connector.connect(
                host=self.db_host,
                user=self.db_user,
                password=self.db_password,
                database=self.db_name
            )
            cursor = conn.cursor()

            did = ligand_info['DID']
            smiles = ligand_info['smiles']
            label_list = json.dumps(ligand_info.get('label_list', []))
            ligand = 1
            inactive = ligand_info.get('inactive', 0)

            cursor.execute(f"SELECT COUNT(*) FROM {self.table_name_ligand} WHERE DID = %s", (did,))
            if cursor.fetchone()[0]:
                cursor.execute(f'''
                    UPDATE {self.table_name_ligand}
                    SET complex_smiles = %s, label_list = %s, ligand = %s, inactive = %s
                    WHERE DID = %s
                ''', (smiles, label_list, ligand, inactive, did))
            else:
                cursor.execute(f'''
                    INSERT INTO {self.table_name_ligand} (DID, complex_smiles, label_list, ligand, inactive)
                    VALUES (%s, %s, %s, %s, %s)
                ''', (did, smiles, label_list, ligand, inactive))

            conn.commit()
            conn.close()
        except Exception as e:
            print(f"[ERROR] Failed to insert ligand info: {e}")