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

class MoleculeSplitter:
    def __init__(self, smiles, marked_ligand_data, irl_list, ga_list=None):
        """
        初始化MoleculeSplitter类，接受SMILES字符串、标记信息以及IRL列表。

        :param smiles: 输入的SMILES字符串
        :param marked_ligand_data: 含有标记的原子数据字典
        :param irl_list: IRL标记列表
        :param ga_list: 可选的广义原子数据字典
        """
        self.smiles = smiles  # 存储分子的SMILES字符串
        self.ga_data = ga_list  # 存储广义原子数据
        self.marked_ligand_data = marked_ligand_data  # 存储标记的配体原子数据
        self.irl_list = irl_list  # 存储IRL标记列表
        self.mol = Chem.MolFromSmiles(smiles)  # 从SMILES字符串生成RDKit分子对象
        self.marked_molecule_data = self.initialize_atom_status()  # 初始化原子的状态信息
        self.IRL_in_ligand = self.get_all_IRL_DIDs()  # 获取配体中的所有IRL的DID
        # print(self.IRL_in_ligand)  # 打印标记的原子信息
        # print(self.IRL_in_ligand)  # 打印IRL DID信息

    def initialize_atom_status(self):
        """
        初始化分子中每个原子的标记数据，默认所有原子都激活。

        :return: 原子的标记状态字典
        """
        atom_status = {}
        # 遍历分子中的每个原子，初始化原子标记状态
        for atom_idx in range(self.mol.GetNumAtoms()):
            atom_status[f'atom_{atom_idx}'] = {
                'IRL_ids': [],  # 初始化IRL_ids为空
                'GA_ids': [],   # 初始化GA_ids为空
                'atom_symbol': self.mol.GetAtomWithIdx(atom_idx).GetSymbol(),  # 获取原子符号
                'active': True,   # 默认所有原子都激活
                'atom_idx': atom_idx  # 保存原子的索引
            }

        # 根据marked_ligand_data填充每个原子的IRL_ids和GA_ids
        for idx, atom in enumerate(self.mol.GetAtoms()):
            atom_idx = f'atom_{idx}'  # 生成原子索引
            if atom_idx in self.marked_ligand_data:
                atom_data = self.marked_ligand_data[atom_idx]
                atom_status[atom_idx]['IRL_ids'] = atom_data.get('IRL_ids', [])  # 填充IRL_ids
                atom_status[atom_idx]['GA_ids'] = atom_data.get('ga_ids', [])  # 填充GA_ids
                atom_status[atom_idx]['atom_symbol'] = atom_data.get('atom_symbol', '')  # 填充原子符号

        return atom_status

    def get_all_IRL_DIDs(self):
        """
        提取配体分子中所有的IRL DID。
        
        :return: 一个包含所有IRL DID的唯一列表
        """
        irl_dids = set()  # 使用集合来避免重复的DID
        marked_ligand_data = self.marked_ligand_data  # 获取标记的原子数据

        # 遍历标记的原子数据，提取IRL_ids中的DID
        for atom_data in marked_ligand_data.values():
            irl_ids = atom_data.get('IRL_ids', [])  # 获取当前原子的IRL_ids
            for irl_id in irl_ids:  # 将IRL_ids中的每个ID添加到irl_dids集合中
                irl_dids.add(irl_id)

        return list(irl_dids)  # 返回IRL DID的列表

    def bfs_find_substructure(self, start_atom_indices):
        """
        使用BFS算法寻找一组原子的所有一阶邻居，包括其所在的GA或普通原子。

        :param start_atom_indices: 起始原子的索引集合（元组），表示一个IRL
        :return: 一个包含所有相关原子索引的集合，表示完整的分子片段
        """
        from collections import deque

        fragment_atoms = set()  # 用于存储片段中的所有原子索引
        visited = set()         # 用于记录已访问的原子
        queue = deque(start_atom_indices)  # 初始化队列，包含起始原子

        while queue:
            current_idx = queue.popleft()  # 获取队列中的原子索引
            if current_idx in visited:     # 如果已访问，跳过
                continue
            visited.add(current_idx)       # 标记为已访问
            fragment_atoms.add(current_idx)  # 添加当前原子到片段集合

            # 获取当前原子的邻居
            if isinstance(current_idx, str):  # 检查索引是否是字符串类型
                atom_idx = int(current_idx.split('_')[1])  # 提取原子编号
            else:
                atom_idx = current_idx  # 如果索引已经是整数，直接使用它

            current_atom = self.mol.GetAtomWithIdx(atom_idx)
            for neighbor in current_atom.GetNeighbors():
                neighbor_idx = f'atom_{neighbor.GetIdx()}'

                # if neighbor_idx in visited:  # 如果邻居已经访问，跳过
                #     continue

                # 获取邻居原子的信息
                neighbor_data = self.marked_molecule_data.get(neighbor_idx, {})
                # print(neighbor_data)
                neighbor_ga_ids = neighbor_data.get('GA_ids', [])
                # print(neighbor_ga_ids)
                if neighbor_ga_ids:  # 如果邻居属于 GA
                    # 找到整个 GA 的所有原子
                    ga_id = neighbor_ga_ids[0]
                    connected_ga_atoms = self._find_connected_ga_atoms(neighbor_idx, ga_id)
                    for ga_atom in connected_ga_atoms:
                        if self.marked_molecule_data[ga_atom]['active']:
                            fragment_atoms.add(ga_atom)
                else:  # 如果邻居不属于 GA
                    if self.marked_molecule_data[neighbor_idx]['active']:
                        fragment_atoms.add(neighbor_idx)

        return fragment_atoms



    def _find_connected_ga_atoms(self, start_atom_idx, ga_id):
        """
        找到与指定原子连接的、属于同一GA实例的所有原子。

        :param start_atom_idx: 起始原子的索引
        :param ga_id: 当前GA的ID
        :return: 一个包含与该GA实例相关的所有原子的索引列表
        """
        connected_ga_atoms = set()  # 用于存储属于同一GA实例的原子
        queue = deque([start_atom_idx])  # 初始化队列
        visited = set()  # 记录已访问的原子

        while queue:
            current_idx = queue.popleft()
            if current_idx in visited:
                continue
            visited.add(current_idx)

            # 目前默认同样的GA只在一个分子里出现一次
            # 检查当前原子是否属于指定的GA实例
            current_data = self.marked_molecule_data.get(current_idx, {})
            if ga_id in current_data.get('GA_ids', []):
                connected_ga_atoms.add(current_idx)

                # 将该原子的邻接原子加入队列
                current_atom = self.mol.GetAtomWithIdx(int(current_idx.split('_')[1]))
                for neighbor in current_atom.GetNeighbors():
                    neighbor_idx = f'atom_{neighbor.GetIdx()}'
                    if neighbor_idx not in visited:
                        queue.append(neighbor_idx)

        return connected_ga_atoms

    
    def match_substructure_and_get_indices(self, substructure_smiles):
        """
        在目标分子中查找所有匹配的子结构，并返回匹配的原子索引列表。

        :param substructure_smiles: 要匹配的子结构的SMILES字符串
        :return: 匹配的原子索引列表（每个元素都是一个匹配的原子索引元组）
        """
        substructure_mol = Chem.MolFromSmiles(substructure_smiles)  # 生成子结构的分子对象
        if substructure_mol is None:
            print(f"无法解析子结构：{substructure_smiles}")
            return []

        # 执行子结构匹配，找出所有匹配的子结构
        matches = self.mol.GetSubstructMatches(substructure_mol)
        if matches:
            # print(f"匹配到的原子索引：{matches}")
            return matches  # 返回所有匹配的原子索引列表
        else:
            print(f"没有找到子结构：{substructure_smiles}")
            return []
        
    def split_molecule_based_on_IRL(self):
        """
        根据IRL标记分割分子，返回分割后的分子列表，每个分子包括原子索引和标记信息。
        """
        split_molecules = []  # 用于存储分割后的分子片段
        IRL_unprocessed_list = self.IRL_in_ligand

        # 遍历IRL列表，按照顺序处理每个IRL
        for IRL_info in self.irl_list:
            IRL_did = IRL_info['DID']  # 获取IRL的DID
            IRL_smiles = IRL_info['complex_smiles']  # 获取IRL的SMILES

            if IRL_did not in IRL_unprocessed_list:
                continue
            # print(f"处理 IRL: {IRL_did}")
            # print(f"IRL SMILES: {IRL_smiles}")

            # 找出所有包含该IRL的原子索引
            atom_index_related_to_IRL = [
                atom_idx
                for atom_idx, atom_data in self.marked_molecule_data.items()
                if IRL_did in atom_data['IRL_ids']
            ]

            # 如果 IRL 中有任何一个原子处于失活状态，则跳过该 IRL
            if any(not self.marked_molecule_data[atom_idx]['active'] for atom_idx in atom_index_related_to_IRL):
                continue

            # 获取分子中 IRL 的具体匹配原子索引
            IRL_index_list = self.match_substructure_and_get_indices(IRL_smiles)
            for IRL_index in IRL_index_list:
                # 确保 IRL 中的原子均为活跃状态
                if any(not self.marked_molecule_data[f'atom_{idx}']['active'] for idx in IRL_index):
                    continue
                formatted_IRL_index = [f'atom_{idx}' for idx in IRL_index]
                # 使用 BFS 找到该 IRL 的相关片段
                fragment_atoms = self.bfs_find_substructure(formatted_IRL_index)
                # 标记片段中的所有原子为失活
                for atom in fragment_atoms:
                    self.mark_inactive(atom)
                
                # 将分割结果添加到分子列表
                # 将 IRL 的元数据与片段原子索引一起保存
                split_molecules.append({
                    'IRL_did': IRL_did,  # 添加 IRL DID
                    'IRL_smiles': IRL_smiles,  # 添加 IRL SMILES
                    'fragment_atoms': fragment_atoms  # 片段的原子索引
                })

        return split_molecules
    
    def get_substructure_smiles_from_custom_indices(self, atom_indices):
        """
        根据自定义格式的原子索引列表（如 'atom_13'）提取分子子结构的完整拓扑结构，
        并返回其SMILES表示。

        :param atom_indices: 包含自定义格式原子索引的列表（如 'atom_13'）
        :return: 子结构的SMILES字符串
        """
        if not atom_indices:
            return ""  # 如果索引列表为空，返回空字符串

        # 提取原子索引的数值部分
        numeric_indices = [int(idx.split('_')[1]) for idx in atom_indices]

        # 创建一个新的分子对象
        new_mol = Chem.RWMol()
        atom_mapping = {}

        # 添加原子到新分子
        for idx in numeric_indices:
            atom = self.mol.GetAtomWithIdx(idx)
            new_idx = new_mol.AddAtom(atom)  # 创建新的原子
            atom_mapping[idx] = new_idx  # 建立原子索引映射

        # 添加键到新分子
        for bond in self.mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            # 如果键的两个原子都在子分子中，添加该键
            if begin_idx in atom_mapping and end_idx in atom_mapping:
                new_mol.AddBond(
                    atom_mapping[begin_idx], atom_mapping[end_idx], bond.GetBondType()
                )

        # 更新分子属性缓存并返回SMILES
        new_mol.UpdatePropertyCache()
        substructure = new_mol.GetMol()
        return Chem.MolToSmiles(substructure)

    def mark_inactive(self, atom_idx):
        """
        将给定的原子标记为不可激活，即已经被处理过，不再作为IRL的起点。
        """
        atom_data = self.marked_molecule_data.get(atom_idx)  # 获取原子的标记数据
        if atom_data:
            atom_data['active'] = False  # 标记该原子为失活状态

    def get_active_atoms(self):
        """
        获取所有处于激活状态的原子索引。

        :return: 激活原子的索引列表
        """
        return [atom_idx for atom_idx, atom_data in self.marked_molecule_data.items() if atom_data['active']]  # 返回激活原子的索引列表

    def get_inactive_atoms(self):
        """
        获取所有处于失活状态的原子索引。

        :return: 失活原子的索引列表
        """
        return [atom_idx for atom_idx, atom_data in self.marked_molecule_data.items() if not atom_data['active']]  # 返回失活原子的索引列表

    def split_and_get_smiles(self):
        """
        分割分子并获取每个片段的SMILES表示。

        :return: 包含每个片段SMILES的列表
        """
        # 调用分割方法，获得片段的原子索引列表
        # split_fragments is a list with element of atom indices
        split_fragments = self.split_molecule_based_on_IRL()
        # print(split_fragments)
        # 获取每个片段的SMILES表示
        fragment_smiles_list = []
        for fragment in split_fragments:
            fragment_atoms = fragment['fragment_atoms']
            fragment_smiles = self.get_substructure_smiles_from_custom_indices(fragment_atoms)

            fragment_smiles_list.append({'fragment_IRL_did': fragment['IRL_did'],
                                         'fragment_IRL_smiles': fragment['IRL_smiles'],
                                         'fragment_smiles': fragment_smiles,
                                         'fragment_atoms': fragment_atoms})
        # print(split_fragments)
        # print(fragment_smiles_list)
        # print(1)
        return fragment_smiles_list

    def visualize_split_molecules(self, split_molecules, save_path=None):
        """
        可视化分割后的分子，并保存生成的图像（如果提供保存路径）。
        """
        num_fragments = len(split_molecules)

        # 创建子图布局
        fig, axes = plt.subplots(nrows=3, ncols=max(num_fragments, 1), figsize=(max(num_fragments, 1) * 4, 14))

        # 第一行：高亮分子图
        mol = Chem.Mol(self.mol)
        atom_colors = {}
        highlight_atoms = []

        # 定义或扩展颜色
        colors = [
            (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
            (1.0, 1.0, 0.0), (1.0, 0.0, 1.0), (0.0, 1.0, 1.0)
        ]
        while len(colors) < num_fragments:
            colors.append((random.random(), random.random(), random.random()))  # 随机生成颜色

        # 分配颜色给片段
        for idx, fragment in enumerate(split_molecules):
            for atom_tag in fragment['fragment_atoms']:
                try:
                    atom_idx = int(atom_tag.split('_')[1])
                    atom_colors[atom_idx] = colors[idx % len(colors)]
                    highlight_atoms.append(atom_idx)
                except (IndexError, ValueError):
                    print(f"跳过无法解析的 atom_tag: {atom_tag}")

        rdDepictor.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DCairo(1000, 1000)
        rdMolDraw2D.PrepareAndDrawMolecule(
            drawer, mol, highlightAtoms=highlight_atoms, highlightAtomColors=atom_colors
        )
        drawer.FinishDrawing()
        png_data = drawer.GetDrawingText()
        image = Image.open(io.BytesIO(png_data))

        axes[0, 0].imshow(image)
        axes[0, 0].axis('off')
        axes[0, 0].set_title("Ligand Structure (With Highlight)", loc='center')

        # 第二行和第三行：片段 SMILES 和结构
        for idx, fragment in enumerate(split_molecules):
            # 第二行
            fragment_smiles = fragment['fragment_smiles']
            fragment_mol = Chem.MolFromSmiles(fragment_smiles)
            if fragment_mol:
                img = Draw.MolToImage(fragment_mol, size=(400, 400))
                axes[1, idx].imshow(img)
                axes[1, idx].axis('off')
                axes[1, idx].text(0.5, -0.2, f"SMILES: {fragment_smiles}", ha='center', fontsize=8, bbox=dict(facecolor='white', alpha=0.7))

            # 第三行
            fragment_IRL_smiles = fragment['fragment_IRL_smiles']
            fragment_IRL_mol = Chem.MolFromSmiles(fragment_IRL_smiles)
            if fragment_IRL_mol:
                img_IRL = Draw.MolToImage(fragment_IRL_mol, size=(200, 200))
                axes[2, idx].imshow(img_IRL)
                axes[2, idx].axis('off')
                axes[2, idx].text(0.5, -0.2, f"SMILES: {fragment_IRL_smiles}\nDID: {fragment['fragment_IRL_did']}", ha='center', fontsize=8, bbox=dict(facecolor='white', alpha=0.7))

        plt.subplots_adjust(hspace=0.3)
        plt.tight_layout(pad=3.0)

        # 保存图像
        if save_path:
            plt.savefig(save_path)
            print(f"Image saved to {save_path}")
        plt.close(fig)
        # plt.show()

def load_db_config(config_path):
    """
    加载数据库配置。
    """
    with open(config_path, 'r') as f:
        config = json.load(f)
    return config['database']


def insert_fragments_to_mysql(source_DID, data, config_path, table_name):
    """
    将片段数据插入到 MySQL 数据库中的新表。
    
    :param source_DID: 片段的来源 DID
    :param data: 列表，每个元素是一个字典，包含 fragment_smiles, fragment_IRL_did, fragment_IRL_smiles
    :param config_path: 数据库配置 JSON 文件的路径
    :param table_name: 目标表名称
    """
    # 读取数据库配置
    db_config = load_db_config(config_path)
    
    # 连接数据库
    conn = mysql.connector.connect(
        host=db_config['host'],
        user=db_config['user'],
        password=db_config['password'],
        database=db_config['database']
    )
    cursor = conn.cursor()
    
    # 创建表
    create_table_sql = f'''
    CREATE TABLE IF NOT EXISTS {table_name} (
        DID INT AUTO_INCREMENT PRIMARY KEY,
        fragment_smiles VARCHAR(255),
        source_DID VARCHAR(50),
        fragment_IRL_did VARCHAR(50),
        fragment_IRL_smiles VARCHAR(50)
    );
    '''
    cursor.execute(create_table_sql)
    
    # 插入数据
    insert_sql = f'''
    INSERT INTO {table_name} (source_DID, fragment_smiles, fragment_IRL_did, fragment_IRL_smiles)
    VALUES (%s, %s, %s, %s);
    '''
    
    for entry in data:
        cursor.execute(insert_sql, (source_DID, entry['fragment_smiles'], entry['fragment_IRL_did'], entry['fragment_IRL_smiles']))
    
    # 提交更改并关闭连接
    conn.commit()
    cursor.close()
    conn.close()
    
    # print(f"数据已成功插入表 {table_name}")
