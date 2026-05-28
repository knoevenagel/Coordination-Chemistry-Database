#!/usr/bin/env python3
"""
Common utility functions for ChemDB processing
"""

import hashlib
from typing import List
from rdkit import Chem
from tools.DID_calculate import calculate_canonical_did

def load_elements_list(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    return [line.strip() for line in lines if line.strip()]

def calculate_did(smiles: str) -> str:
    return calculate_canonical_did(smiles)

def is_coordination_complex(smiles: str, metal_elements: List[str]) -> bool:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in metal_elements:
                return True
        return False
    except:
        return False 