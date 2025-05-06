import csv
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Draw

class FindGA:
    def __init__(self, smiles, object_DID):
        self.smiles = smiles
        self.DID = object_DID
        self.mol = Chem.MolFromSmiles(smiles)
        if not self.mol:
            raise ValueError(f"无效的 SMILES 字符串: {smiles}")
        self.marked_atoms = []  # 添加属性来存储被标记的原子索引

    def get_ring_smiles(self, ring_size):
        """查找所有指定大小的环并返回其 SMILES 结构式"""
        rings = self.find_all_rings()
        specified_size_rings = [ring for ring in rings if len(ring) == ring_size]

        ring_smiles = []
        for ring in specified_size_rings:
            submol = self.get_ring_from_indices(ring)
            smiles = Chem.MolToSmiles(submol, canonical=True)
            ring_smiles.append(smiles)

        return ring_smiles
    # 三四元环的smiles也需要典范化过程
    def get_three_and_four_membered_rings(self):
        """同时查找三元环和四元环并返回其 SMILES 结构式"""
        three_membered_smiles = self.get_ring_smiles(3)
        four_membered_smiles = self.get_ring_smiles(4)
        return three_membered_smiles, four_membered_smiles

    def find_all_rings(self):
        # 找到所有环结构
        # 这里存在的问题是，可能找到的环很奇怪，并没有芳香性
        return Chem.GetSymmSSSR(self.mol)

    def find_aromatic_rings(self):
        """查找所有芳香环"""
        rings = self.find_all_rings()
        aromatic_rings = [ring for ring in rings if all(self.mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)]
        return aromatic_rings

    def get_ring_from_indices(self, ring):
        """根据环的原子索引提取环的完整子分子"""
        new_mol = Chem.RWMol()
        atom_mapping = {}
        
        # 添加原子
        for idx in ring:
            atom = self.mol.GetAtomWithIdx(idx)
            new_idx = new_mol.AddAtom(atom)
            atom_mapping[idx] = new_idx
        
        # 添加键
        for bond in self.mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if begin_idx in atom_mapping and end_idx in atom_mapping:
                new_mol.AddBond(atom_mapping[begin_idx], atom_mapping[end_idx], bond.GetBondType())
        
        new_mol.UpdatePropertyCache()
        return new_mol.GetMol()
    
    def get_aromatic_ring_smiles(self):
        """查找所有芳香环并返回其 SMILES 结构式"""
        rings = self.find_all_rings()  # 获取所有环
        aromatic_rings = [
            ring for ring in rings if all(self.mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        ]  # 筛选出芳香环

        aromatic_ring_smiles = []
        for ring in aromatic_rings:
            # 获取完整的环结构
            submol = self.get_ring_from_indices(ring)  # 调用新的函数
            if submol is None:
                continue  # 如果获取的子分子无效，跳过此环
        # 这里重复了smiles和mol的转变过程一次，是因为直接转成smiles无法典范化
            # 转换为标准化 SMILES
            smiles = Chem.MolToSmiles(submol, canonical=True)
            ring_mol = Chem.MolFromSmiles(smiles)
            if ring_mol is None:
                continue
            ring_smiles = Chem.MolToSmiles(ring_mol, canonical=True)
            if ring_smiles:  # 确保 SMILES 不为空
                aromatic_ring_smiles.append(ring_smiles)  # 添加到结果列表

        # 返回唯一的 SMILES 列表
        return list(set(aromatic_ring_smiles))


    def find_shared_atoms(self, aromatic_rings):
        """查找所有共享原子的稠环"""
        shared_atoms = set()
        for i in range(len(aromatic_rings)):
            for j in range(i + 1, len(aromatic_rings)):
                if set(aromatic_rings[i]).intersection(set(aromatic_rings[j])):
                    shared_atoms.update(aromatic_rings[i])
                    shared_atoms.update(aromatic_rings[j])
        return shared_atoms

    def mark_GA_in_mol(self):

        marked_atoms = set()

        # 标记三元环
        for ring in self.find_three_membered_rings():
            marked_atoms.update(ring)

        # 标记四元环
        for ring in self.find_four_membered_rings():
            marked_atoms.update(ring)

        # 查找所有芳香环
        aromatic_rings = self.find_aromatic_rings()

        # 查找稠环的共享原子
        shared_atoms = self.find_shared_atoms(aromatic_rings)

        # 标记非稠环芳香环
        for ring in aromatic_rings:
            if not set(ring).intersection(shared_atoms):
                marked_atoms.update(ring)
            #print(ring.MolToSmiles())

        self.marked_atoms = list(marked_atoms)  # 保存被标记的原子索引

        # 设置 GA 属性
        for atom in self.mol.GetAtoms():
            atom.SetProp("GA", "1" if atom.GetIdx() in marked_atoms else "0")

        return self.mol  # 返回标记后的分子对象


    def visualize(self):
        # 可视化分子并高亮标记的原子
        highlighted_atoms = self.marked_atoms
        colors = {idx: (1.0, 0.0, 0.0) for idx in highlighted_atoms}  # 红色标记

        # 创建原子编号的标签字典
        atom_labels = {idx: str(idx) for idx in highlighted_atoms}  # 将原子索引转换为字符串

        img = Draw.MolToImage(
            self.mol, 
            size=(300, 300), 
            highlightAtoms=highlighted_atoms, 
            highlightAtomColors=colors,
            highlightAtomLabels=atom_labels  # 显示原子编号
        )
        return img

def get_smiles_DID_list(csv_file_path):
    smiles_with_did = []
    with open(csv_file_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            did = row["DID"]  # 假设你的 CSV 中有一列是 DID
            smiles = row["complex_smiles"]
            smiles_with_did.append((did, smiles))  # 将 DID 和 SMILES 作为元组添加到列表中
    return smiles_with_did


def find_GA_list(smiles_and_did_list, output_csv, non_kekulizable_csv=False):
    GA_dict = {}  # 使用字典来存储 SMILES 和对应的 DID 列表
    non_kekulizable_smiles = []  # 存储不能凯库勒化的 SMILES

    # 使用 tqdm 包装循环以显示进度条
    for DID_smiles_pair in tqdm(smiles_and_did_list, desc="Processing SMILES"):
        mol_DID = DID_smiles_pair[0]
        mol_smiles = DID_smiles_pair[1]
        ga_finder = FindGA(mol_smiles, mol_DID)  # 创建 GA 查找器
        
        # 获取芳香环的 SMILES
        aromatic_ring_smiles = ga_finder.get_aromatic_ring_smiles()
        
        # 处理芳香环
        for ring_smiles in aromatic_ring_smiles:
            ring_mol = Chem.MolFromSmiles(ring_smiles)
            if ring_mol is None:
                non_kekulizable_smiles.append(ring_smiles)  # 记录无效的 SMILES
                continue
            
            # 检查重原子数目是否超过7
            heavy_atom_count = ring_mol.GetNumHeavyAtoms()
            if heavy_atom_count > 7:
                continue  # 跳过重原子数目大于7的环

            try:
                Chem.Kekulize(ring_mol, clearAromaticFlags=True)
                if ring_smiles not in GA_dict:
                    GA_dict[ring_smiles] = set()  # 使用集合去除重复的 DID
                GA_dict[ring_smiles].add(mol_DID)  # 添加 DID 到集合中

            except Exception:
                non_kekulizable_smiles.append(ring_smiles)  # 无法凯库勒化

        # 获取三元环和四元环的 SMILES
        three_membered_smiles, four_membered_smiles = ga_finder.get_three_and_four_membered_rings()

        # 处理三元环
        for ring_smiles in three_membered_smiles:
            ring_mol = Chem.MolFromSmiles(ring_smiles)
            if ring_mol is None:
                non_kekulizable_smiles.append(ring_smiles)  # 记录无效的 SMILES
                continue
            
            heavy_atom_count = ring_mol.GetNumHeavyAtoms()
            if heavy_atom_count > 7:
                continue  # 跳过重原子数目大于7的环

            # 将三元环直接添加到 GA_dict
            if ring_smiles not in GA_dict:
                GA_dict[ring_smiles] = set()
            GA_dict[ring_smiles].add(mol_DID)

        # 处理四元环
        for ring_smiles in four_membered_smiles:
            ring_mol = Chem.MolFromSmiles(ring_smiles)
            if ring_mol is None:
                non_kekulizable_smiles.append(ring_smiles)  # 记录无效的 SMILES
                continue
            
            heavy_atom_count = ring_mol.GetNumHeavyAtoms()
            if heavy_atom_count > 7:
                continue  # 跳过重原子数目大于7的环

            # 将四元环直接添加到 GA_dict
            if ring_smiles not in GA_dict:
                GA_dict[ring_smiles] = set()
            GA_dict[ring_smiles].add(mol_DID)

    # 将能凯库勒化的结果转换为 DataFrame 并保存到 CSV
    GA_list = []
    for ring_smiles, dids in GA_dict.items():
        GA_list.append({"GA_SMILES": ring_smiles, "GA_source": ', '.join(dids)})
    
    if GA_list:  # 如果有可用的凯库勒化结果
        df = pd.DataFrame(GA_list)
        df.to_csv(output_csv, index=False)

    if non_kekulizable_csv:
        non_kekulizable_df = pd.DataFrame(non_kekulizable_smiles, columns=["Non-Kekulizable SMILES"])
        non_kekulizable_df.to_csv(non_kekulizable_csv, index=False)




# 示例使用
if __name__ == "__main__":
    test_smiles = "C12=NC=CC=C1C=CC=C2"
    # 示例使用
    input_csv_path = './ligand_data/ligand_PdZnDatabase.csv'  # 替换为你的输入文件路径
    output_csv_path = './ligand_data/GA_PdZnDatabase.csv'  # 替换为你的输出文件路径
    non_kekulizable_csv_path = "/Users/weilinzou/Code/DatabaseProject/project_data/error_file/non_kekulizable_smiles.csv"
    smiles_and_did_list = get_smiles_DID_list(input_csv_path)
    find_GA_list(smiles_and_did_list, output_csv_path,non_kekulizable_csv=non_kekulizable_csv_path)

    # input_csv = '/home/weilin/ScientificProject/scientific/ligand_process/aromatic_rings.csv'
    # kekulizable_csv = '/home/weilin/ScientificProject/scientific/ligand_process/aromatic_GA.csv'
    # non_kekulizable_csv = '/home/weilin/ScientificProject/scientific/ligand_process/non_kekulizable_smiles.csv'

    # mol_finder = FindGA(test_smiles)
    # print(mol_finder.analyze())
    # visual_img = mol_finder.visualize()
    # visual_img.show()
