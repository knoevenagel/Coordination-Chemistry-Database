import re
import os
import sys
import csv
import hashlib
from collections import deque
from rdkit import Chem
from rdkit.Chem import AllChem      
from rdkit.Chem import Draw
from rdkit.Chem import GetPeriodicTable
from Utils.utils import load_elements_list
from Utils.DID_calculate import calculate_canonical_did


class ChemicalComplex:
    def __init__(self, pubchem_info_dict ,metal_element_list, error_file="error.out"):
        self.cid = pubchem_info_dict['cid']
        # print(self.cid)
        self.smiles = pubchem_info_dict['smiles']
        #self.gibberish = pubchem_info_dict['gibberish']
        self.corcomp = pubchem_info_dict['corcomp']
        # mol = Chem.MolFromSmiles(smiles)
        # self.smiles = Chem.MolToSmiles(mol,canonical = True)
        self.inactive = 0
        #self.illegitimate = 0
        #self.metalcompound = 0
        #self.questionable = 0
        self.did = self.calculate_DID_from_fingerprint(self.smiles)

        # print(self.did)
        self.metal_element_list = metal_element_list
        self.error_file = error_file
        self.complex_info = {
            "complex_info": {
                "DID": self.did,  
                "complex_smiles": self.smiles,
                #"gibberish": self.gibberish,
                "corcomp": self.corcomp,
                "inactive": 0
            },
            "central_metal_info": [],
            "neighbor_info": []
        }
        if self.did is None:
            print('SB SMILES')
            return
        
        #print(self.questionable)
        #self.did = self.generate_did_from_cid(self.cid)  # 生成 DID


        self.metal_info_list = self.extract_metal_content()
        self.metal_did_list = [self.generate_metal_did(symbol) for symbol in self.metal_info_list]

        self.process_metal_structure()
        self.process_complex()
        # questionable的边界条件可以通过修改参数来改变
        self.questionable_detection()
        self.is_metal_only(self.smiles)
        self.complex_info["complex_info"]["inactive"] = self.inactive
        # print(aaa)
        #self.merge_neighbor_info()
        
    # @staticmethod
    # def generate_did_from_cid(cid):
    #     """
    #     根据给定的 CID 生成对应的 DID。

    #     Args:
    #         cid (str): 化合物的 CID。

    #     Returns:
    #         str: 生成的 DID。
    #     """
    #     return 'D' + cid

        
    def generate_metal_did(self, symbol):
        if not symbol:
            print("Element symbol is not provided.")
            return None
        
        # 创建 PeriodicTable 实例
        pt = Chem.GetPeriodicTable()

        # 使用 PeriodicTable 实例获取元素符号对应的原子序数
        atomic_number = pt.GetAtomicNumber(symbol)

        # 检查原子序数是否有效
        if atomic_number is None:
            print("Element symbol is invalid.")
            return None

        # 在原子序数前加上 'D' 来生成 metal_did
        metal_did = f"D{atomic_number}E"
        return metal_did

    
    def mark_smiles(self):
        '''
        Mark atoms of the input molecule.
        Input: SMILES
        Return: RDKit molecule object
        '''
        input_smiles = self.smiles
        input_mol = Chem.MolFromSmiles(input_smiles)
        if input_mol is None:
            with open(self.error_file, "a") as f:
                f.write(f"Error creating molecule from SMILES '{input_smiles}'\n")
            # 返回一个默认的空分子对象，以便代码继续执行
            return Chem.MolFromSmiles('')
        for atom in input_mol.GetAtoms():
            atom.SetProp("atomNote", str(atom.GetIdx()))
        return input_mol


    def extract_metal_content(self):
        '''
        Return: match metal list such as ['Sn', 'Sn']
        Remark that this function is very important
        注意双金属的情况
        '''
        # 正则表达式为关键部分，最好注意这里
        #pattern = re.compile(r'\[(\d*?)([A-Za-z]+)(\+|\-)?\d*\]') # 修改正则表达式以匹配括号中内容
        pattern = re.compile(r'\[([^\]]+)\]')
        matches = pattern.findall(self.smiles)
        #print(matches)
        metal_matches = []
        # for match in matches:
        #     # 获取括号中的金属符号
        #     content_inside_brackets = match
        #     metal_symbol = ''.join([char for char in content_inside_brackets if char.isalpha()])
            
        #     # 检查提取出来的金属符号是否在金属元素列表中
        #     if metal_symbol in self.metal_element_list:
        #         metal_matches.append(metal_symbol)
        for match in matches:
            # 去除非大小写字母的字符，只保留字母部分
            clean_match = ''.join([char for char in match if char.isalpha()])
            
            # 按照大写字母进行拆分成列表
            metal_symbols = re.findall('[A-Z][a-z]*', clean_match)
            
            # 检查每个拆分出的金属符号是否在金属元素列表中
            for metal_symbol in metal_symbols:
                if metal_symbol in self.metal_element_list:
                    metal_matches.append(metal_symbol)
                    # 只要匹配到一个金属元素就跳出循环，因为我们只需检查是否存在
                    break
        #print(metal_matches)
        return metal_matches
    
    def questionable_detection(self,ligand_num_boundary = 11):
        """
        这里默认参数的boundary是可以修改的
        """
        cc_smiles = self.smiles
        ligand_num = cc_smiles.count('.')
        if ligand_num > ligand_num_boundary:
            #self.questionable = 1
            self.inactive = 3
        # print(ligand_num)

    def calculate_DID_from_fingerprint(self, smiles):
        result = calculate_canonical_did(smiles)
        #print(result)
        return result

    # def calculate_DID_from_fingerprint(self, smiles):
    #     try:
    #         # 使用 RDKit 验证 SMILES 字符串
    #         mol = Chem.MolFromSmiles(smiles)
    #         if mol is None:
    #             print(f"Invalid SMILES: '{smiles}'")
    #             self.illegitimate = 1  # 设置 self.illegitimate 为 1
    #             return None
            
    #         # 获取标准化的 SMILES
    #         canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
    #         # 生成分子指纹，参数还需要修改
    #         canonical_mol = Chem.MolFromSmiles(canonical_smiles)
    #         fingerprint = AllChem.GetMorganFingerprintAsBitVect(canonical_mol, radius=2, nBits=2048)
    #         fingerprint_bytes = fingerprint.ToBitString().encode()
    #         # 使用 hashlib 生成哈希值

    #         hash_value = hashlib.sha256(fingerprint_bytes).hexdigest()
    #         hash_decimal = int(hash_value, 16)
    #         truncated_hash_decimal = hash_decimal % (10 ** 10)
    #         truncated_hash_decimal_str = 'D' + str(truncated_hash_decimal)
            
    #         return truncated_hash_decimal_str
    #     except Exception as e:
    #         print(f"Error calculating DID for SMILES '{smiles}': {e}")
    #         self.illegitimate = 1  # 设置 self.illegitimate 为 1
    #         return None

    def generate_ligand_smiles(self, coordinate_atom_index, central_metal_index, original_mol, central_metals_set):
        """
        这段代码记得加注释
        """
        bonds = original_mol.GetBonds()
        connected_atoms = set()
        visited = set()
        queue = deque([coordinate_atom_index])
        while queue:
            current_atom_idx = queue.popleft()
            visited.add(current_atom_idx)
            if current_atom_idx == central_metal_index:
                continue
            connected_atoms.add(current_atom_idx)
            for bond in bonds:
                if current_atom_idx == bond.GetBeginAtomIdx():
                    neighbor_idx = bond.GetEndAtomIdx()
                elif current_atom_idx == bond.GetEndAtomIdx():
                    neighbor_idx = bond.GetBeginAtomIdx()
                else:
                    continue
                if neighbor_idx in central_metals_set:
                    continue
                if neighbor_idx not in visited:
                    queue.append(neighbor_idx)
        ligand_index_list = list(connected_atoms)
        new_mol = self.get_mol_from_indices(original_mol, ligand_index_list)
        coordinating_atom_new_index = ligand_index_list.index(coordinate_atom_index)
        ligand_smiles = Chem.MolToSmiles(new_mol, canonical=True)
        did = self.calculate_DID_from_fingerprint(ligand_smiles)
        return ligand_smiles, coordinating_atom_new_index, did

    def find_metals(self, smiles):
        """
        寻找SMILES结构中的金属元素
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print("Invalid SMILES input!")
            return
        metals_found = []
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol in self.metal_element_list:
                metals_found.append(symbol)
        return metals_found

    def get_mol_from_indices(self, original_mol, indices):
        """
        这段代码也要加注释
        """
        new_mol = Chem.RWMol()
        atom_mapping = {}
        for idx in indices:
            atom = original_mol.GetAtomWithIdx(idx)
            new_idx = new_mol.AddAtom(atom)
            atom_mapping[idx] = new_idx
        for bond in original_mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if begin_idx in indices and end_idx in indices:
                new_mol.AddBond(atom_mapping[begin_idx], atom_mapping[end_idx], bond.GetBondType())
        new_mol.UpdatePropertyCache()
        return new_mol.GetMol()

    def get_neighbors_info(self, central_metals):
        """
        代码记得加注释
        """
        marked_rdkit_mol = self.mark_smiles()
        bonds = marked_rdkit_mol.GetBonds()
        metal_atom_indices = set()
        for i in range(len(central_metals)):
            central_metal = central_metals[i]
            atom_index_tuples = marked_rdkit_mol.GetSubstructMatches(central_metal)
            metal_atom_indices.update([atom_tuple[0] for atom_tuple in atom_index_tuples])
        #print(metal_atom_indices)
        for atom_index in metal_atom_indices:
            central_metal_elements = marked_rdkit_mol.GetAtomWithIdx(atom_index).GetSymbol()
            central_metal_elements_did = self.generate_metal_did(central_metal_elements)
            #print(central_metal_elements)
            for bond in bonds:
                if bond.GetBeginAtomIdx() == atom_index:
                    coordinating_atom_index = bond.GetEndAtomIdx()
                elif bond.GetEndAtomIdx() == atom_index:
                    coordinating_atom_index = bond.GetBeginAtomIdx()
                else:
                    continue
                ligand_smiles, coordinating_atom_new_index, ligand_did = self.generate_ligand_smiles(
                    coordinating_atom_index, atom_index, marked_rdkit_mol, metal_atom_indices)
                bond_type = bond.GetBondType()
                self.complex_info["neighbor_info"].append({
                    "coordinating_atom_index": coordinating_atom_new_index,
                    "central_atom": central_metal_elements,
                    "central_atom_did": central_metal_elements_did,
                    "bond_type": bond_type,
                    "ligand_smiles": ligand_smiles,
                    "ligand_did": ligand_did
                })

    def is_metal_only(self, smiles):
        """
        检测是否是金属互化物，不需要返回值，只在全部原子都是金属时设置 self.inactive 为 2
        """
        try:
            # 解析 SMILES 字符串生成分子对象
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"Invalid SMILES: '{smiles}'")
                self.inactive = 1  # 设置 self.inactive 为 1
                return 
            
            # 获取分子中的所有原子
            atoms = mol.GetAtoms()

            # 检查所有原子是否都是金属元素
            for atom in atoms:
                symbol = atom.GetSymbol()
                if symbol not in self.metal_element_list:
                    return  # 如果有非金属元素，直接退出，不做任何设置

            # 如果所有原子都是金属元素
            self.inactive = 2

        except Exception as e:
            print(f"Error processing SMILES '{smiles}': {e}")

    
    def calculate_valence(self, bonds, metal_atom_index):
        """
        这个函数还需要修改，现在没啥用
        """
        total_bond_state = 0
        for bond in bonds:
            atom1_index = bond.GetBeginAtom().GetIdx()
            atom2_index = bond.GetEndAtom().GetIdx()
            if atom1_index == metal_atom_index or atom2_index == metal_atom_index:
                bond_order = bond.GetBondTypeAsDouble()
                total_bond_state += bond_order
        return total_bond_state


    def get_metal_info(self, central_metals):
        """
        这个函数需要写注释+优化
        """
        marked_rdkit_mol = self.mark_smiles()
        bonds = marked_rdkit_mol.GetBonds()
        metal_atom_indices = set()
        #print(metal_atom_indices)
        for central_metal in central_metals:
            atom_index_tuples = marked_rdkit_mol.GetSubstructMatches(central_metal)
            metal_atom_indices.update([atom_tuple[0] for atom_tuple in atom_index_tuples])
        #print(metal_atom_indices)
        for atom_index in metal_atom_indices:
            #print(atom_index)
            #central_metal_elements = self.find_metals(Chem.MolToSmiles(central_metals[0]))
            central_metal_elements = marked_rdkit_mol.GetAtomWithIdx(atom_index).GetSymbol()
            #atom_number = Chem.PeriodicTable.GetAtomicNumber(central_metal_elements)
            # print(type(central_metal_elements))
            # print(central_metal_elements)
            #print(atom_number)
            #print(type(central_metal_elements))
            metal_did = self.generate_metal_did(central_metal_elements)
            central_metal_info = {
                "atom_index": atom_index,
                "central_metal": central_metal_elements,
                "valence": 0,
                "metal_did": metal_did
            }
            total_bond_state = self.calculate_valence(bonds, atom_index)
            central_metal_info["valence"] = total_bond_state
            self.complex_info["central_metal_info"].append(central_metal_info)

    
    def process_metal_structure(self):
        metal_symbols = self.metal_info_list
        if metal_symbols:
            metal_substructures = [Chem.MolFromSmarts('[' + symbol + ']') for symbol in metal_symbols]
            self.get_metal_info(metal_substructures)    

    def process_covalent_structure(self):
        metal_symbols = self.metal_info_list
        if metal_symbols:
            metal_substructures = [Chem.MolFromSmarts('[' + symbol + ']') for symbol in metal_symbols]
            #print(metal_substructures)
            self.get_neighbors_info(metal_substructures)

    def process_ion_structure(self):
        metal_symbols = self.metal_info_list
        metal_symbols_did = self.metal_did_list
        #print(metal_symbols)
        if metal_symbols:
            substructure_pattern = re.compile(r'[^\.]+')
            substructures = substructure_pattern.findall(self.smiles)
            for substructure in substructures:
                if any(metal in substructure for metal in metal_symbols):
                    self.process_covalent_structure()
                else:
                    did = self.calculate_DID_from_fingerprint(substructure)
                    ligand_info = {'coordinating_atom_index': None, 'ligand_smiles': substructure, 'bond_type': 'IONIC',
                                'central_atom': metal_symbols, 'central_atom_did': metal_symbols_did, 'ligand_did': did}
                    self.complex_info['neighbor_info'].append(ligand_info)

    def merge_neighbor_info(self):
        merged_info = self.complex_info.copy()  # 复制输入字典以避免修改原始数据
        merged_neighbor_info_list = []
        seen = set()  # 用于跟踪已经见过的配体信息
        try:
            for ligand_info in self.complex_info["neighbor_info"]:
                # 处理 central_atom_did 可能是列表的情况
                central_atom_did = ligand_info.get("central_atom_did", "")
                if isinstance(central_atom_did, list):
                    central_atom_did = ",".join(map(str, central_atom_did))  # 将列表转换为逗号分隔的字符串
                tuple_key = (ligand_info.get("coordinating_atom_index", ""), central_atom_did, ligand_info.get("ligand_did", ""))
                
                # 如果配体信息相同，则不重复添加
                if tuple_key not in seen:
                    merged_neighbor_info_list.append(ligand_info)
                    seen.add(tuple_key)
                else:
                    # 如果配体信息已经存在，则合并字典
                    for existing_ligand_info in merged_neighbor_info_list:
                        if (existing_ligand_info.get("coordinating_atom_index") == ligand_info.get("coordinating_atom_index") and
                            existing_ligand_info.get("central_atom_did") == central_atom_did and
                            existing_ligand_info.get("ligand_did") == ligand_info.get("ligand_did")):
                            # 合并字典
                            existing_ligand_info.update(ligand_info)
                            break
            merged_info["neighbor_info"] = merged_neighbor_info_list
        except Exception as e:
            print(f"Error processing SMILES: {e}")
        return merged_info


    def is_multi_molecule(self):
        return '.' in self.smiles


    def process_complex(self):
        if self.is_multi_molecule():
            return self.process_ion_structure()
        else:
            return self.process_covalent_structure()


# 主函数的调用也需要更新
# 注意，这里有一个事情还没有解决，可能配体会被重复地插入到字典里
if __name__ == "__main__":
    metal_path = './data/metal_list.txt'
    metal_list = load_elements_list(metal_path)
    # print(metal_list)
    # print(len(metal_list))
    test_smiles = 'CC[Sn+2]CC.C(=O)[O-].C(=O)[O-].[Sn]'
    test_smiles2 = 'CCCC[Sn](CCCC)CCC(C)COC(=C)C(=O)OC'
    test_smiles3 = '[2H][C-]([2H])[2H].[2H]C([2H])([2H])C(C1=CC=NC=C1)(C2=NN(C=C2)CC(F)F)O.[Li+].CC(C1=CC=NC=C1)C2=NN(C=C2)CC(F)F.C1=CN=CC=C1C(C2=NN(C=C2)CC(F)F)O.C1=CN=CC=C1C(=O)C2=NN(C=C2)CC(F)F.C1=CN=CC=C1I.C1=CN(N=C1C=O)CC(F)F.C1=C(NN=C1)C=O.C(C(F)F)OS(=O)(=O)C(F)(F)F.O=[Mn]=O.[Zn]'
    test_smiles4 = 'C1=CC=C2C(=C1)C=CC=C2O.C1=CC=C2C(=C1)C=CC=C2O.[Mn]'
    cid_file_name = "/home/weilin/ScientificProject/scientific/PubChemData/processed_data/Mn.csv"  # 请替换为您的 CID 文件名
    # 创建 ChemicalComplex 实例
    print('ok')
    test_info_dict = {'cid': '111111', 'smiles': 'CC1=C2CC(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(CC1OC(=O)[C@@H]([C@@H](C)C=C(C)C)O)O)OC(=O)C)(CO4)O)O)C.[Ac]','corcomp':True}
    complex_instance = ChemicalComplex(test_info_dict, metal_list)
    #print(complex_instance.complex_info)
    merged_dict = complex_instance.merge_neighbor_info()
    print(merged_dict)
    print(complex_instance.inactive)
    # 使用 add_did_info 函数处理 DID 信息
    #complex_instance_with_did = complex_instance.add_did_info(cid_file_name)
    # print(complex_instance_with_did)

