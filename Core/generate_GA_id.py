import csv
import hashlib
import os

csv.field_size_limit(1048576)  # 设置字段大小限制为1MB

def generate_GA_id(smiles):
    """
    根据GA的SMILES生成一个唯一的身份标识符（G开头的六位十进制编号）
    """
    hash_object = hashlib.md5(smiles.encode())
    hash_hex = hash_object.hexdigest()
    first_six_hex = hash_hex[:12]
    hash_decimal = int(first_six_hex, 16)
    final_ga_id = f"G{hash_decimal % 1000000:06d}"
    return final_ga_id

def process_and_check_GA_file(input_csv_path, output_csv_path, error_csv_path):
    """
    处理输入的CSV文件，生成GA的ID，检查重复，并保存结果。
    如果有重复ID，将不重复的部分保存到指定输出文件，并将重复的SMILES保存到另一个文件。

    参数：
        input_csv_path: 输入CSV文件路径
        output_csv_path: 输出CSV文件路径（不重复部分）
        error_csv_path: 报错SMILES的输出CSV文件路径
    """
    output_rows = []  # 存储不重复的行
    error_rows = []  # 存储重复的SMILES行
    ga_id_set = set()  # 用于检测重复的GA ID

    # 处理输入文件，生成GA ID并检查重复
    with open(input_csv_path, mode='r', newline='', encoding='utf-8') as infile:
        csv_reader = csv.DictReader(infile)
        
        for row in csv_reader:
            smiles = row["GA_SMILES"]
            ga_id = generate_GA_id(smiles)

            # 检查GA ID是否重复
            if ga_id in ga_id_set:
                error_rows.append([smiles, ga_id])
            else:
                ga_id_set.add(ga_id)
                output_rows.append([smiles, ga_id])

    # 将不重复的行保存到输出文件
    with open(output_csv_path, mode='w', newline='', encoding='utf-8') as outfile:
        csv_writer = csv.writer(outfile)
        csv_writer.writerow(["GA_SMILES", "GA_ID"])
        csv_writer.writerows(output_rows)

    # 如果有重复，将重复的SMILES保存到错误文件
    if error_rows:
        with open(error_csv_path, mode='w', newline='', encoding='utf-8') as errorfile:
            csv_writer = csv.writer(errorfile)
            csv_writer.writerow(["GA_SMILES", "GA_ID"])
            csv_writer.writerows(error_rows)

        print(f"发现重复的GA ID，重复SMILES已保存到: {error_csv_path}")
    else:
        print("没有发现重复的GA ID")

    print(f"处理完成，不重复的GA ID已保存到: {output_csv_path}")

# 示例调用
if __name__ == "__main__":
    # 文件路径设置
    input_csv = './ligand_data/GA_PdZnDatabase.csv'
    output_csv = './GA_IRL_data/SemiGA_with_id.csv'
    error_csv = './GA_IRL_data/SemiGA_duplicates.csv'

    # 调用函数
    process_and_check_GA_file(input_csv, output_csv, error_csv)