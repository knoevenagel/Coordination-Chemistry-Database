#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
from collections import deque
from rdkit import Chem
from rdkit.Chem.rdmolops import GetAdjacencyMatrix
from utils import load_elements_list
from tools.DID_calculate import calculate_canonical_did

_PT = Chem.GetPeriodicTable()

class ChemicalComplex:
    def __init__(self, pubchem_info_dict, metal_element_list, error_file="error.out"):
        self.cid = pubchem_info_dict['cid']
        self.smiles = pubchem_info_dict['smiles']
        self.corcomp = pubchem_info_dict['corcomp']
        self.inactive = 0
        self.error_file = error_file

        self._mol = Chem.MolFromSmiles(self.smiles)
        if self._mol is None:
            with open(self.error_file, "a") as f:
                f.write(f"Invalid SMILES at init '{self.smiles}'\n")
            self.inactive = 1
            self.did = None
            self.metal_element_list = metal_element_list
            self.complex_info = {
                "complex_info": {"DID": None, "complex_smiles": self.smiles, "corcomp": self.corcomp, "inactive": self.inactive},
                "central_metal_info": [],
                "neighbor_info": []
            }
            return

        adj = GetAdjacencyMatrix(self._mol)
        self._adj = [set() for _ in range(self._mol.GetNumAtoms())]
        for i in range(adj.shape[0]):
            nbrs = adj[i].nonzero()[0]
            self._adj[i].update(int(j) for j in nbrs)

        self.did = self.calculate_DID_from_fingerprint(self.smiles)
        self.metal_element_list = metal_element_list
        self.complex_info = {
            "complex_info": {"DID": self.did, "complex_smiles": self.smiles, "corcomp": self.corcomp, "inactive": 0},
            "central_metal_info": [],
            "neighbor_info": []
        }
        if self.did is None:
            return

        self.metal_info_list = self.extract_metal_content()
        self.metal_did_list = [self.generate_metal_did(symbol) for symbol in self.metal_info_list]

        self._metal_atom_indices = set()
        self._metal_atom_symbol = {}
        self._metal_atom_did = {}
        self._neighbor_cache = {}  # key: (metal_idx, neighbor_idx) -> dict

        self.process_metal_structure()
        self.process_complex()
        self.questionable_detection()
        self.is_metal_only(self.smiles)
        self.complex_info["complex_info"]["inactive"] = self.inactive

    def generate_metal_did(self, symbol):
        if not symbol:
            return None
        atomic_number = _PT.GetAtomicNumber(symbol)
        if not atomic_number:
            return None
        return f"D{atomic_number}E"

    def mark_smiles(self):
        if self._mol is None:
            return Chem.MolFromSmiles('')
        mol = Chem.Mol(self._mol)
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)
        return mol

    def extract_metal_content(self):
        pattern = re.compile(r'\[([^\]]+)\]')
        matches = pattern.findall(self.smiles)
        metal_matches = []
        for match in matches:
            clean_match = ''.join([char for char in match if char.isalpha()])
            metal_symbols = re.findall('[A-Z][a-z]*', clean_match)
            for metal_symbol in metal_symbols:
                if metal_symbol in self.metal_element_list:
                    metal_matches.append(metal_symbol)
                    break
        return metal_matches

    def questionable_detection(self, ligand_num_boundary=11):
        cc_smiles = self.smiles
        ligand_num = cc_smiles.count('.')
        if ligand_num > ligand_num_boundary:
            self.inactive = 3

    def calculate_DID_from_fingerprint(self, smiles):
        return calculate_canonical_did(smiles)

    def generate_ligand_smiles(self, coordinate_atom_index, central_metal_index, original_mol, central_metals_set):
        connected_atoms = set()
        visited = set()
        queue = deque([coordinate_atom_index])
        while queue:
            cur = queue.popleft()
            if cur in visited:
                continue
            visited.add(cur)
            if cur == central_metal_index:
                continue
            connected_atoms.add(cur)
            for nb in self._adj[cur]:
                if nb in central_metals_set or nb in visited:
                    continue
                queue.append(nb)
        ligand_index_list = list(connected_atoms)
        coordinating_atom_new_index = ligand_index_list.index(coordinate_atom_index)
        ligand_smiles = Chem.MolFragmentToSmiles(original_mol, atomsToUse=ligand_index_list, isomericSmiles=True, canonical=True)
        did = self.calculate_DID_from_fingerprint(ligand_smiles)
        return ligand_smiles, coordinating_atom_new_index, did

    def find_metals(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return
        metals_found = []
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol in self.metal_element_list:
                metals_found.append(symbol)
        return metals_found

    def get_neighbors_info(self, central_metals):
        mol = self._mol
        if not self._metal_atom_indices:
            tmp_indices = set()
            for cm in central_metals:
                for t in mol.GetSubstructMatches(cm):
                    tmp_indices.add(t[0])
            for idx in tmp_indices:
                sym = mol.GetAtomWithIdx(idx).GetSymbol()
                self._metal_atom_indices.add(idx)
                self._metal_atom_symbol[idx] = sym
                self._metal_atom_did[idx] = self.generate_metal_did(sym)

        for metal_idx in self._metal_atom_indices:
            for nb in self._adj[metal_idx]:
                key = (metal_idx, nb)
                if key in self._neighbor_cache:
                    self.complex_info["neighbor_info"].append(self._neighbor_cache[key])
                    continue
                ligand_smiles, coord_new_idx, ligand_did = self.generate_ligand_smiles(
                    nb, metal_idx, mol, self._metal_atom_indices
                )
                bond = mol.GetBondBetweenAtoms(metal_idx, nb)
                info = {
                    "coordinating_atom_index": coord_new_idx,
                    "central_atom": self._metal_atom_symbol[metal_idx],
                    "central_atom_did": self._metal_atom_did[metal_idx],
                    "bond_type": bond.GetBondType() if bond is not None else None,
                    "ligand_smiles": ligand_smiles,
                    "ligand_did": ligand_did
                }
                self._neighbor_cache[key] = info
                self.complex_info["neighbor_info"].append(info)

    def is_metal_only(self, smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                self.inactive = 1
                return
            for atom in mol.GetAtoms():
                if atom.GetSymbol() not in self.metal_element_list:
                    return
            self.inactive = 2
        except Exception:
            pass

    def calculate_valence(self, bonds, metal_atom_index):
        total_bond_state = 0
        for bond in bonds:
            a = bond.GetBeginAtom().GetIdx()
            b = bond.GetEndAtom().GetIdx()
            if a == metal_atom_index or b == metal_atom_index:
                total_bond_state += bond.GetBondTypeAsDouble()
        return total_bond_state

    def get_metal_info(self, central_metals):
        mol = self._mol
        bonds = mol.GetBonds()
        metal_atom_indices = set()
        for cm in central_metals:
            for t in mol.GetSubstructMatches(cm):
                metal_atom_indices.add(t[0])
        for idx in metal_atom_indices:
            sym = mol.GetAtomWithIdx(idx).GetSymbol()
            metal_did = self.generate_metal_did(sym)
            info = {"atom_index": idx, "central_metal": sym, "valence": 0, "metal_did": metal_did}
            info["valence"] = self.calculate_valence(bonds, idx)
            self.complex_info["central_metal_info"].append(info)

    def process_metal_structure(self):
        metal_symbols = self.metal_info_list
        if metal_symbols:
            metal_substructures = [Chem.MolFromSmarts('[' + symbol + ']') for symbol in metal_symbols]
            self.get_metal_info(metal_substructures)
            mol = self._mol
            tmp_indices = set()
            for cm in metal_substructures:
                for t in mol.GetSubstructMatches(cm):
                    tmp_indices.add(t[0])
            for idx in tmp_indices:
                sym = mol.GetAtomWithIdx(idx).GetSymbol()
                self._metal_atom_indices.add(idx)
                self._metal_atom_symbol[idx] = sym
                self._metal_atom_did[idx] = self.generate_metal_did(sym)

    def process_covalent_structure(self):
        metal_symbols = self.metal_info_list
        if metal_symbols:
            metal_substructures = [Chem.MolFromSmarts('[' + symbol + ']') for symbol in metal_symbols]
            self.get_neighbors_info(metal_substructures)

    def process_ion_structure(self):
        metal_symbols = self.metal_info_list
        metal_symbols_did = self.metal_did_list
        if metal_symbols:
            substructure_pattern = re.compile(r'[^\.]+')
            substructures = substructure_pattern.findall(self.smiles)
            for substructure in substructures:
                if any(metal in substructure for metal in metal_symbols):
                    self.process_covalent_structure()
                else:
                    did = self.calculate_DID_from_fingerprint(substructure)
                    ligand_info = {
                        'coordinating_atom_index': None,
                        'ligand_smiles': substructure,
                        'bond_type': 'IONIC',
                        'central_atom': metal_symbols,
                        'central_atom_did': metal_symbols_did,
                        'ligand_did': did
                    }
                    self.complex_info['neighbor_info'].append(ligand_info)

    def merge_neighbor_info(self):
        merged = {}
        for li in self.complex_info["neighbor_info"]:
            cad = li.get("central_atom_did", "")
            if isinstance(cad, list):
                cad = ",".join(map(str, cad))
            key = (li.get("coordinating_atom_index", ""), cad, li.get("ligand_did", ""))
            if key not in merged:
                merged[key] = li
            else:
                merged[key].update(li)
        out = self.complex_info.copy()
        out["neighbor_info"] = list(merged.values())

        # For ionic multi-molecule inputs like '....[Fe].[As]', build a simplified
        # central_metal_info as [{'Fe': '0'}, {'As': '0'}] in order of appearance.
        # This mirrors original behavior for downstream consumers that expect this summary.
        central_simple = []
        if self.is_multi_molecule():
            for part in self.smiles.split('.'):
                part = part.strip()
                m = re.fullmatch(r"\[([^\]]+)\]", part)
                if not m:
                    continue
                content = m.group(1)
                # Keep only letters and grab leading element symbol (e.g., Fe from Fe+2)
                letter_only = ''.join(ch for ch in content if ch.isalpha())
                sm = re.match(r"([A-Z][a-z]?)", letter_only)
                if sm:
                    symbol = sm.group(1)
                    central_simple.append({symbol: '0'})
        if central_simple:
            out["central_metal_info"] = central_simple
        return out

    def is_multi_molecule(self):
        return '.' in self.smiles

    def process_complex(self):
        if self.is_multi_molecule():
            return self.process_ion_structure()
        else:
            return self.process_covalent_structure()


if __name__ == "__main__":
    metal_path = '/Users/weilinzou/Code/DatabaseProject/ChemDB/Resources/metal_list.txt'
    metal_list = load_elements_list(metal_path)
    print('ok')
    test_info_dict = {'cid': '111111', 'smiles': 'CC1=C2CC(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(CC1OC(=O)[C@@H]([C@@H](C)C=C(C)C)O)O)OC(=O)C)(CO4)O)O)C.[Ac]','corcomp':True}
    complex_instance = ChemicalComplex(test_info_dict, metal_list)
    merged_dict = complex_instance.merge_neighbor_info()
    print(merged_dict)
    print(complex_instance.inactive)