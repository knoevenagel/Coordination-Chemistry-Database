import os
import json
import pandas as pd
import mysql.connector
from mysql.connector import pooling
from tqdm import tqdm
from Utils.load_database_config import load_mysql_config
from concurrent.futures import ProcessPoolExecutor, as_completed


class ProcessPubchemRawData:
    def __init__(self, db_config):
        """
        初始化类，传入 MySQL 数据库配置（字典），并建立连接池
        """
        self.db_config = db_config
        self.checked_tables = set()
        self.pool = pooling.MySQLConnectionPool(
            pool_name="mypool",
            pool_size=8,
            **db_config
        )

    def create_pubchem_data_table(self, csv_file_path, table_name):
        """
        从 CSV 文件读取表头字段，并创建 MySQL 表：
        - 'cid' 为 VARCHAR(255) 主键
        - 其余字段为 TEXT 类型
        """
        if table_name in self.checked_tables:
            return

        try:
            df = pd.read_csv(csv_file_path, nrows=0)
            field_names = list(df.columns)

            if "cid" not in field_names:
                raise ValueError("字段 'cid' 不存在，无法设置主键")

            conn = self.pool.get_connection()
            cursor = conn.cursor()

            column_defs = []
            for col in field_names:
                if col == "cid":
                    column_defs.append(f"{col} VARCHAR(255) PRIMARY KEY")
                else:
                    column_defs.append(f"{col} TEXT")
            columns = ', '.join(column_defs)

            sql = f"CREATE TABLE IF NOT EXISTS {table_name} ({columns})"
            cursor.execute(sql)
            conn.commit()
            self.checked_tables.add(table_name)
            print(f"表 `{table_name}` 创建成功")
        except Exception as e:
            print(f"建表失败: {e}")
        finally:
            if 'cursor' in locals():
                cursor.close()
            if 'conn' in locals():
                conn.close()

    def insert_to_mysql(self, data, table_name):
        """
        插入单条数据到 MySQL（使用连接池连接）
        """
        try:
            conn = self.pool.get_connection()
            cursor = conn.cursor()

            columns = ', '.join(data.keys())
            placeholders = ', '.join(['%s'] * len(data))
            sql = f"REPLACE INTO {table_name} ({columns}) VALUES ({placeholders})"
            values = tuple(data.values())

            cursor.execute(sql, values)
            conn.commit()
        except mysql.connector.Error as err:
            print(f"插入失败: {err}")
            with open("failed_rows.log", "a") as f:
                f.write(json.dumps(data, ensure_ascii=False) + "\n")
        finally:
            if 'cursor' in locals():
                cursor.close()
            if 'conn' in locals():
                conn.close()

    def insert_data_parallel_to_mysql(self, data_list, table_name):
        """
        串行插入数据并显示 tqdm 进度条（使用连接池连接）
        """
        if not data_list:
            return

        for row in tqdm(data_list, desc=f"插入数据到 {table_name}", unit="条"):
            try:
                conn = self.pool.get_connection()
                cursor = conn.cursor()

                columns = ', '.join(row.keys())
                placeholders = ', '.join(['%s'] * len(row))
                sql = f"REPLACE INTO {table_name} ({columns}) VALUES ({placeholders})"
                values = tuple(row.values())

                cursor.execute(sql, values)
                conn.commit()
            except mysql.connector.Error as err:
                print(f"插入失败: {err}")
                with open("failed_rows.log", "a") as f:
                    f.write(json.dumps(row, ensure_ascii=False) + "\n")
            finally:
                if 'cursor' in locals():
                    cursor.close()
                if 'conn' in locals():
                    conn.close()

    def extract_csv_data_parallel(self, file_path, batch_size=1000):
        """
        批量读取 CSV 数据为字典列表
        """
        df_iterator = pd.read_csv(file_path, chunksize=batch_size)
        all_data = []
        for chunk in df_iterator:
            rows = chunk.to_dict(orient="records")
            all_data.extend(rows)
        return all_data

    def process_folder(self, folder_path, table_name, batch_size=1000):
        """
        遍历文件夹下所有 CSV 文件，并插入 MySQL（显示 tqdm 进度条）
        """
        csv_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.csv')]
        for file_path in tqdm(csv_files, desc="处理CSV文件"):
            print(f"\n开始处理文件: {file_path}")
            data = self.extract_csv_data_parallel(file_path, batch_size)
            self.insert_data_parallel_to_mysql(data, table_name)
            print(f"处理完成: {file_path}")


    @staticmethod
    def fetch_cid_isosmiles_batch(db_config, table_name, offset, batch_size):
        """
        提取一批 cid 和 isosmiles，直接返回 DataFrame 行（dict 列表）
        """
        try:
            conn = mysql.connector.connect(**db_config)
            cursor = conn.cursor()
            cursor.execute(f"SELECT cid, isosmiles FROM {table_name} LIMIT {batch_size} OFFSET {offset}")
            rows = cursor.fetchall()
            cursor.close()
            conn.close()

            return [{"cid": str(cid), "isosmiles": isosmiles}
                    for cid, isosmiles in rows if cid and isosmiles]
        except Exception as e:
            print(f"offset={offset} 提取失败: {e}")
            return []

    def fetch_data_parallel_from_mysql(self, table_name, batch_size=1000, max_workers=4):
            """
            并发提取所有 cid 和 isosmiles，直接返回 pandas.DataFrame
            """
            try:
                conn = self.pool.get_connection()
                cursor = conn.cursor()
                cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
                total_rows = cursor.fetchone()[0]
                cursor.close()
                conn.close()

                offsets = list(range(0, total_rows, batch_size))
                task_args = [(self.db_config, table_name, offset, batch_size) for offset in offsets]

                all_data = []
                with ProcessPoolExecutor(max_workers=max_workers) as executor:
                    futures = [executor.submit(self.fetch_cid_isosmiles_batch, *args) for args in task_args]
                    for future in tqdm(as_completed(futures), total=len(futures), desc="提取 cid 和 isosmiles"):
                        all_data.extend(future.result())

                # 直接转 DataFrame
                df = pd.DataFrame(all_data)
                return df
            except Exception as e:
                print(f"提取失败: {e}")
                return pd.DataFrame(columns=["cid", "isosmiles"])


if __name__ == "__main__":
    mysql_config_path = '/Users/weilinzou/Code/DatabaseProject/ChemDB/Config/mysql_config.json'
    pubchem_data_path = '/Users/weilinzou/Code/DatabaseProject/project_data/PubchemTestData'
    pubchem_table_sample = '/Users/weilinzou/Code/DatabaseProject/project_data/pubchem_sample.csv'
    pubchem_table_name = 'PubChemRawData'

    mysql_config = load_mysql_config(mysql_connect_config_path=mysql_config_path)
    processor = ProcessPubchemRawData(mysql_config)

    # # 1. 建表
    # processor.create_pubchem_data_table(pubchem_table_sample, pubchem_table_name)

    # # 2. 插入数据（带进度条）
    # processor.process_folder(pubchem_data_path, pubchem_table_name)

    df = processor.fetch_data_parallel_from_mysql("PubChemRawData", batch_size=1000)
    print(df.head())
    print("提取数据条数:", df.shape[0])
    print(df.memory_usage(deep=True))  # 各列占用情况
    print("总内存占用：", df.memory_usage(deep=True).sum() / 1024**2, "MB")