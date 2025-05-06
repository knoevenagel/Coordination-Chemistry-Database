import pandas as pd
from rdkit import Chem
from tqdm import tqdm
# 这个代码需要修改下，让进度条的显示更加合理
class ExtractIRLObject:
    def __init__(self, csv_path: str) -> None:
        """
        初始化方法，从 CSV 文件加载分子数据
        :param csv_path: 包含 'DID', 'complex_smiles', 'GAC' 列的 CSV 文件路径
        """
        self.csv_path = csv_path
        self.molecule_list = self.load_csv()
        self.compute_heavy_elements()  # 计算非氢元素列表

    def load_csv(self):
        """从 CSV 文件中加载数据并返回包含必要字段的字典列表"""
        df = pd.read_csv(self.csv_path, usecols=['DID', 'complex_smiles', 'GAC'])
        molecule_list = df.to_dict(orient='records')
        # 按 GAC 从小到大排序
        molecule_list.sort(key=lambda x: x['GAC'])
        return molecule_list

    def compute_heavy_elements(self):
        """
        计算每个分子包含的非氢元素，并将其添加到分子字典中
        """
        for item in tqdm(self.molecule_list, desc="Calculating non-hydrogen elements"):
            smiles = item['complex_smiles']
            
            if not isinstance(smiles, str) or pd.isna(smiles):
                item['heavy_elements'] = []  # 如果 SMILES 无效，则添加空列表
                continue

            mol = Chem.MolFromSmiles(smiles)
            
            if mol is None:
                item['heavy_elements'] = []  # 无法解析的 SMILES
                continue

            # 提取非氢元素并去重
            heavy_elements = {atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetSymbol() != 'H'}
            item['heavy_elements'] = sorted(heavy_elements)  # 排序以统一格式
            
    def contains_required_elements(self, current_elements, target_elements):
        """
        检查 target_elements 是否包含 current_elements 中的所有元素
        :param current_elements: 当前分子包含的非氢元素列表
        :param target_elements: 目标分子包含的非氢元素列表
        :return: 如果 target_elements 包含所有 current_elements，则返回 True；否则返回 False
        """
        return set(current_elements).issubset(set(target_elements))

    def check_substructure(self, mol1, mol2):
        """
        检查 mol1 是否是 mol2 的子结构
        :param mol1: RDKit Mol 对象 (子结构)
        :param mol2: RDKit Mol 对象 (被检查对象)
        :return: 如果 mol1 是 mol2 的子结构，返回 True，否则返回 False
        """
        return mol2.HasSubstructMatch(mol1)
    
    def filter_by_GAC(self):
        """
        筛选出 GAC 大于等于 2 的分子
        :return: 仅包含 GAC >= 2 的分子列表
        """
        self.molecule_list = [item for item in self.molecule_list if item['GAC'] >= 2]

    def extract_TIRL(self):
        """
        根据提取规则，提取 TIRL 并检查是否被包含在其他分子结构中
        :return: 处理后的分子列表，更新 'flag', 'substructure', 'parent_structure'
        """
        # 先筛选 GAC >= 2 的分子
        self.filter_by_GAC()

        for i in tqdm(range(len(self.molecule_list)), desc="Processing TIRL extraction"):
            current_item = self.molecule_list[i]
            current_GAC = current_item['GAC']
            current_smiles = current_item['complex_smiles']
            current_mol = Chem.MolFromSmiles(current_smiles)
            
            if not current_mol:
                continue

            # 检查当前分子是否已经被标记为子结构（即 flag 为 True），如果是，则跳过
            if current_item.get('flag', False):
                continue  # 如果 flag 为 True，跳过当前分子

            for j in tqdm(range(i + 1, len(self.molecule_list)), desc=f"Checking substructures for item {i+1}/{len(self.molecule_list)}", leave=False):
                target_item = self.molecule_list[j]
                if target_item['GAC'] <= current_GAC:
                    continue

                target_smiles = target_item['complex_smiles']
                target_mol = Chem.MolFromSmiles(target_smiles)
                
                if not target_mol:
                    continue

                if self.check_substructure(current_mol, target_mol) and self.contains_required_elements(current_item['heavy_elements'], target_item['heavy_elements']):
                    target_item['flag'] = True
                    if 'substructure' not in target_item:
                        target_item['substructure'] = []
                    target_item['substructure'].append(current_item['DID'])

                    if 'parent_structure' not in current_item:
                        current_item['parent_structure'] = []
                    current_item['parent_structure'].append(target_item['DID'])

            # 在每次外层循环结束时，打印 'flag' 为 True 的分子数量
            self.print_flag_count()

        return self.molecule_list

    def print_flag_count(self):
        """
        检查并打印 'flag' 为 True 的分子数量
        """
        flag_true_count = sum(1 for item in self.molecule_list if item.get('flag', False))
        print(f"Number of molecules with flag set to True: {flag_true_count}")

if __name__ == "__main__":
    # 这里后续需要处理下，让IRL确定为在某个GA下的IRL，需要把代码联动一下
    ligand_csv_path = './ligand_data/ligand_GAC_PdZnDatabase.csv'  # 替换为实际 CSV 文件路径
    IRL_output_csv_path = './GA_IRL_data/IRL_output.csv'  # 替换为所需输出路径
    extractor = ExtractIRLObject(ligand_csv_path)
    #print(extractor.molecule_list[10000])

    result = extractor.extract_TIRL()

    # 打印结果或保存到新 CSV 文件
    result_df = pd.DataFrame(result)
    result_df.to_csv(IRL_output_csv_path, index=False)  # 替换为所需输出路径
