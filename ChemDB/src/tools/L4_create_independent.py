#!/usr/bin/env python3
"""
独立L4提取器：每个IRL独立提取，互不干扰
不使用原子active/inactive状态，每个IRL都能完整提取
"""

import sys
from rdkit import Chem
from collections import deque

sys.path.insert(0, "/home/chenghua/ZhuoLi/ChemDB_restructured/ChemDB/src")
from tools.DID_calculate import calculate_canonical_did

class IndependentL4Extractor:
    """独立L4提取器：每个IRL独立处理，不共享原子状态
    
    兼容原始MoleculeSplitter接口，可直接替换使用
    """
    
    def __init__(self, smiles, marked_ligand_data, irl_list, ga_list=None, source_did=None):
        self.smiles = smiles
        self.marked_ligand_data = marked_ligand_data
        self.irl_list = irl_list
        self.ga_data = ga_list or []  # 兼容原始接口的ga_data属性名
        self.ga_list = ga_list or []
        self.source_did = source_did  # 添加source_did支持
        
        self.mol = Chem.MolFromSmiles(smiles)
        if self.mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        # 获取配体中的所有IRL（排除CC）
        self.irls_in_ligand = self._get_all_irl_dids()
        
        # 兼容属性：原始接口使用IRL_in_ligand
        self.IRL_in_ligand = self.irls_in_ligand
        
        # 缓存提取结果
        self._extraction_results = None
        self._fragments_cache = None
    
    def _get_all_irl_dids(self):
        """提取配体中所有IRL的DID（排除CC）"""
        CC_DID = 'D218966891838592'
        irl_dids = set()
        for atom_data in self.marked_ligand_data.values():
            for irl_id in atom_data.get('IRL_ids', []):
                if irl_id != CC_DID:  # 排除CC
                    irl_dids.add(irl_id)
        return list(irl_dids)
    
    def _find_connected_ga_atoms(self, start_atom_idx, ga_id):
        """找到GA组的所有原子"""
        connected = set()
        queue = deque([start_atom_idx])
        visited = set()
        
        while queue:
            current_idx = queue.popleft()
            if current_idx in visited:
                continue
            visited.add(current_idx)
            
            current_data = self.marked_ligand_data.get(current_idx, {})
            if ga_id in current_data.get('ga_ids', []):
                connected.add(current_idx)
                
                atom_idx = int(current_idx.split('_')[1])
                atom = self.mol.GetAtomWithIdx(atom_idx)
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = f'atom_{neighbor.GetIdx()}'
                    if neighbor_idx not in visited:
                        queue.append(neighbor_idx)
        
        return connected
    
    def _bfs_expand_from_irl(self, irl_atom_indices):
        """从IRL原子BFS扩展到一级邻居"""
        fragment_atoms = set()
        visited = set()
        queue = deque(irl_atom_indices)
        
        while queue:
            current_idx = queue.popleft()
            if current_idx in visited:
                continue
            visited.add(current_idx)
            fragment_atoms.add(current_idx)
            
            # 获取邻居
            atom_idx = int(current_idx.split('_')[1]) if isinstance(current_idx, str) else current_idx
            current_atom = self.mol.GetAtomWithIdx(atom_idx)
            
            for neighbor in current_atom.GetNeighbors():
                neighbor_idx = f'atom_{neighbor.GetIdx()}'
                neighbor_data = self.marked_ligand_data.get(neighbor_idx, {})
                neighbor_ga_ids = neighbor_data.get('ga_ids', [])
                
                if neighbor_ga_ids:
                    # 包含整个GA组
                    ga_id = neighbor_ga_ids[0]
                    ga_atoms = self._find_connected_ga_atoms(neighbor_idx, ga_id)
                    fragment_atoms.update(ga_atoms)
                else:
                    # 普通邻居原子
                    fragment_atoms.add(neighbor_idx)
        
        return fragment_atoms
    
    def _get_substructure_smiles(self, atom_indices):
        """从原子索引生成SMILES"""
        if not atom_indices:
            return ""
        
        numeric_indices = [int(idx.split('_')[1]) for idx in atom_indices]
        
        new_mol = Chem.RWMol()
        atom_mapping = {}
        
        for idx in numeric_indices:
            atom = self.mol.GetAtomWithIdx(idx)
            new_idx = new_mol.AddAtom(atom)
            atom_mapping[idx] = new_idx
        
        for bond in self.mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if begin_idx in atom_mapping and end_idx in atom_mapping:
                new_mol.AddBond(
                    atom_mapping[begin_idx], 
                    atom_mapping[end_idx], 
                    bond.GetBondType()
                )
        
        new_mol.UpdatePropertyCache()
        return Chem.MolToSmiles(new_mol.GetMol())
    
    def extract_l4_for_each_irl(self):
        """为每个IRL独立提取L4片段"""
        if self._extraction_results is not None:
            return self._extraction_results
        
        results = []
        
        for irl_info in self.irl_list:
            irl_did = irl_info['DID']
            irl_smiles = irl_info.get('complex_smiles', irl_info.get('SMILES', ''))
            
            # 只处理在配体中的IRL（CC已在_get_all_irl_dids中排除）
            if irl_did not in self.irls_in_ligand:
                continue
            
            # 找到标记为这个IRL的所有原子
            irl_atoms = [
                atom_idx for atom_idx, atom_data in self.marked_ligand_data.items()
                if irl_did in atom_data['IRL_ids']
            ]
            
            if not irl_atoms:
                continue
            
            # 尝试子结构匹配
            irl_mol = Chem.MolFromSmiles(irl_smiles)
            if irl_mol is None:
                results.append({
                    'irl_did': irl_did,
                    'irl_smiles': irl_smiles,
                    'status': 'SMILES_PARSE_FAILED',
                    'l4_fragments': []
                })
                continue
            
            matches = self.mol.GetSubstructMatches(irl_mol)
            if not matches:
                results.append({
                    'irl_did': irl_did,
                    'irl_smiles': irl_smiles,
                    'status': 'NO_SUBSTRUCTURE_MATCH',
                    'l4_fragments': []
                })
                continue
            
            # 为每个匹配位置生成L4片段
            l4_fragments = []
            for match in matches:
                try:
                    formatted_indices = [f'atom_{idx}' for idx in match]
                    fragment_atoms = self._bfs_expand_from_irl(formatted_indices)
                    fragment_smiles = self._get_substructure_smiles(fragment_atoms)
                    
                    if not fragment_smiles:
                        continue
                    
                    fragment_did = calculate_canonical_did(fragment_smiles)
                    
                    l4_fragments.append({
                        'fragment_smiles': fragment_smiles,
                        'fragment_did': fragment_did,
                        'atom_indices': list(fragment_atoms)
                    })
                except Exception as e:
                    # 跳过生成失败的片段
                    continue
            
            results.append({
                'irl_did': irl_did,
                'irl_smiles': irl_smiles,
                'status': 'SUCCESS',
                'l4_fragments': l4_fragments
            })
        
        self._extraction_results = results
        return results
    
    @property
    def fragments(self):
        """
        兼容原始MoleculeSplitter接口的属性
        返回平铺的L4片段列表，格式与原始方法相同
        """
        if self._fragments_cache is not None:
            return self._fragments_cache
        
        extraction_results = self.extract_l4_for_each_irl()
        
        # 转换格式：嵌套结构 -> 平铺列表
        flat_fragments = []
        for result in extraction_results:
            if result['status'] == 'SUCCESS':
                for l4 in result['l4_fragments']:
                    flat_fragments.append({
                        'fragment_smiles': l4['fragment_smiles'],
                        'fragment_DID': l4['fragment_did'],  # 键名转换
                        'fragment_IRL_did': result['irl_did'],
                        'fragment_IRL_smiles': result['irl_smiles'],
                        'fragment_atoms': l4.get('atom_indices', [])
                    })
        
        self._fragments_cache = flat_fragments
        return flat_fragments
    
    # 以下是兼容原始接口的辅助方法
    def get_all_IRL_DIDs(self):
        """兼容原始接口的方法"""
        return self.irls_in_ligand
    
    def mark_inactive(self, atom_idx):
        """兼容方法（空实现，因为独立版本不使用active状态）"""
        pass
    
    def get_active_atoms(self):
        """兼容方法：返回所有原子（因为不使用active状态）"""
        return list(self.marked_ligand_data.keys())
    
    def get_inactive_atoms(self):
        """兼容方法：返回空列表（因为不使用active状态）"""
        return []

