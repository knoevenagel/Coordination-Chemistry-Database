#!/usr/bin/env python3
"""
Fingerprint utility functions for ChemDB processing
"""

import base64
import logging
from typing import Tuple, Optional
from functools import lru_cache
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from utils import calculate_did

# Disable RDKit logging
RDLogger.DisableLog('rdApp.*')
RDLogger.DisableLog('rdMol.*')
RDLogger.DisableLog('rdDecomposition.*')

logger = logging.getLogger(__name__)

@lru_cache(maxsize=8)
def _get_morgan_generator(radius: int = 2, fp_size: int = 2048):
    """Cache Morgan fingerprint generator (only radius and fp_size matter)."""
    return rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=fp_size)

def generate_morgan_fingerprint(smiles: str) -> Tuple[Optional[str], int]:
    """
    Generate Morgan fingerprint for SMILES.
    Returns (base64_of_packed_bytes, on_bits_count) or (None, 0) if failed.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, 0

        fp = _get_morgan_generator().GetFingerprint(mol)
        packed_bytes = DataStructs.BitVectToBinaryText(fp)  # already packed bits
        fp_b64 = base64.b64encode(packed_bytes).decode("utf-8")
        return fp_b64, fp.GetNumOnBits()

    except Exception:
        return None, 0

def process_single_smiles(smiles: str) -> Tuple[Optional[str], Optional[str], int]:
    """
    Process single SMILES string and return (did, fingerprint, bit_count)
    
    Args:
        smiles: SMILES string
        calculate_did_func: Function to calculate DID
        
    Returns:
        Tuple of (did, fingerprint, bit_count) or (None, None, 0) if failed
    """
    try:
        if not smiles or smiles == 'nan' or smiles == '':
            return None, None, 0
        
        # Calculate DID
        did = calculate_did(smiles)
        if not did:
            return None, None, 0
        
        # Generate fingerprint
        fingerprint, bit_count = generate_morgan_fingerprint(smiles)
        if not fingerprint:
            return None, None, 0
        
        return did, fingerprint, bit_count
        
    except Exception as e:
        logger.error(f"处理SMILES失败: {smiles}, 错误: {e}")
        return None, None, 0
