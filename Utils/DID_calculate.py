# DID_calculate.py

from rdkit import Chem
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
import hashlib

# 初始化全局指纹生成器，只生成一次以节省资源
morgan_generator = GetMorganGenerator(radius=2, fpSize=4096)

def calculate_canonical_did(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Invalid SMILES: '{smiles}'")
            return None
        
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        canonical_mol = Chem.MolFromSmiles(canonical_smiles)

        # 使用新版 MorganGenerator 生成指纹
        fingerprint = morgan_generator.GetFingerprint(canonical_mol)
        fingerprint_bytes = fingerprint.ToBitString().encode()

        hash_value = hashlib.sha256(fingerprint_bytes).hexdigest()
        hash_decimal = int(hash_value, 16)
        truncated_hash_decimal = hash_decimal % (10 ** 10)
        truncated_hash_decimal_str = 'D' + str(truncated_hash_decimal)
        
        return truncated_hash_decimal_str
    except Exception as e:
        print(f"Error calculating DID for SMILES '{smiles}': {e}")
        return None

if __name__ == "__main__":
    smiles = "N=N"
    print(calculate_canonical_did(smiles))