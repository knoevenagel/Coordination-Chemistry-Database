"""
MolCLR FastAPI Server: Extract Molecular Embeddings via REST API

This API serves molecular embeddings from pre-trained MolCLR models.

Endpoint:
    POST /api/embeddings
    
Request Body:
    {
        "smiles": ["CCO", "c1ccccc1"],  # List of SMILES strings
        "model_type": "gin"              # Optional: "gin" or "gcn" (default: "gin")
    }

Response:
    {
        "embeddings_base64": "...",      # Base64-encoded numpy array (float32)
        "shape": [2, 512],               # Shape of the embedding array
        "dtype": "float32",              # Data type
        "valid_count": 2,                # Number of successfully processed molecules
        "invalid_smiles": []             # List of invalid SMILES (if any)
    }

Usage:
    uvicorn api:app --host 0.0.0.0 --port 8000
"""

import os
import sys
import base64
import shutil
import subprocess
from io import BytesIO
from typing import List, Optional

import numpy as np
import torch
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rdkit import Chem
from rdkit.Chem.rdchem import BondType as BT
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from torch_geometric.data import Data, Batch


# Constants for molecular featurization
ATOM_LIST = list(range(1, 119))
CHIRALITY_LIST = [
    Chem.rdchem.ChiralType.CHI_UNSPECIFIED,
    Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW,
    Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW,
    Chem.rdchem.ChiralType.CHI_OTHER
]
BOND_LIST = [BT.SINGLE, BT.DOUBLE, BT.TRIPLE, BT.AROMATIC]
BONDDIR_LIST = [
    Chem.rdchem.BondDir.NONE,
    Chem.rdchem.BondDir.ENDUPRIGHT,
    Chem.rdchem.BondDir.ENDDOWNRIGHT
]

# Global model cache
_models = {}


def smiles_to_graph(smiles: str) -> Optional[Data]:
    """Convert a SMILES string to a PyTorch Geometric Data object."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    mol = Chem.AddHs(mol)
    
    type_idx = []
    chirality_idx = []
    for atom in mol.GetAtoms():
        type_idx.append(ATOM_LIST.index(atom.GetAtomicNum()))
        chirality_idx.append(CHIRALITY_LIST.index(atom.GetChiralTag()))
    
    x1 = torch.tensor(type_idx, dtype=torch.long).view(-1, 1)
    x2 = torch.tensor(chirality_idx, dtype=torch.long).view(-1, 1)
    x = torch.cat([x1, x2], dim=-1)
    
    row, col, edge_feat = [], [], []
    for bond in mol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        row += [start, end]
        col += [end, start]
        bond_type = BOND_LIST.index(bond.GetBondType())
        bond_dir = BONDDIR_LIST.index(bond.GetBondDir())
        edge_feat.append([bond_type, bond_dir])
        edge_feat.append([bond_type, bond_dir])
    
    edge_index = torch.tensor([row, col], dtype=torch.long)
    edge_attr = torch.tensor(np.array(edge_feat), dtype=torch.long)
    
    return Data(x=x, edge_index=edge_index, edge_attr=edge_attr)


def check_model_files_exist(base_dir: str) -> bool:
    """Check if model checkpoint files exist."""
    tmp_dir = os.path.join(base_dir, 'tmp')
    ckpt_dir = os.path.join(tmp_dir, 'ckpt')
    # Check for at least one model checkpoint
    gin_ckpt = os.path.join(ckpt_dir, 'pretrained_gin', 'checkpoints', 'model.pth')
    gcn_ckpt = os.path.join(ckpt_dir, 'pretrained_gcn', 'checkpoints', 'model.pth')
    return os.path.exists(gin_ckpt) or os.path.exists(gcn_ckpt)


def prepare_model_files(base_dir: str) -> str:
    """Prepare model files by cloning MolCLR repo and copying ckpt folder.
    
    Returns:
        Path to ckpt directory
    """
    tmp_dir = os.path.join(base_dir, 'tmp')
    ckpt_dir = os.path.join(tmp_dir, 'ckpt')
    molclr_repo_dir = os.path.join(tmp_dir, 'MolCLR')
    
    # Check if ckpt already exists
    if os.path.exists(ckpt_dir):
        return ckpt_dir
    
    print("Model files not found. Preparing...")
    
    # Create tmp directory if not exists
    os.makedirs(tmp_dir, exist_ok=True)
    
    # Clone MolCLR repository if not exists
    if not os.path.exists(molclr_repo_dir):
        print("Cloning MolCLR repository...")
        subprocess.run(
            ['git', 'clone', 'https://github.com/yuyangw/MolCLR.git'],
            cwd=tmp_dir,
            check=True
        )
    
    # Copy ckpt folder from cloned repo to tmp
    src_ckpt = os.path.join(molclr_repo_dir, 'ckpt')
    if os.path.exists(src_ckpt):
        print(f"Copying ckpt folder from {src_ckpt} to {ckpt_dir}...")
        shutil.copytree(src_ckpt, ckpt_dir)
        print("Model files prepared successfully.")
    else:
        raise FileNotFoundError(f"ckpt folder not found in cloned repo: {src_ckpt}")
    
    return ckpt_dir


def get_model(model_type: str, device: str = 'cpu'):
    """Load and cache a pre-trained MolCLR model."""
    cache_key = f"{model_type}_{device}"
    
    if cache_key not in _models:
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        
        # Prepare model files if not exist
        ckpt_base = prepare_model_files(base_dir)
        
        # Add MolCLR directory to path for model imports
        molclr_dir = os.path.join(base_dir, 'tmp', 'MolCLR')
        if molclr_dir not in sys.path:
            sys.path.insert(0, molclr_dir)
        
        checkpoint_path = os.path.join(
            ckpt_base, f'pretrained_{model_type}', 'checkpoints', 'model.pth'
        )
        
        if not os.path.exists(checkpoint_path):
            raise FileNotFoundError(f"Checkpoint not found: {checkpoint_path}")
        
        if model_type == 'gin':
            from models.ginet_molclr import GINet
            model = GINet(num_layer=5, emb_dim=300, feat_dim=512, drop_ratio=0, pool='mean')
        elif model_type == 'gcn':
            from models.gcn_molclr import GCN
            model = GCN(num_layer=5, emb_dim=300, feat_dim=512, drop_ratio=0, pool='mean')
        else:
            raise ValueError(f"Unknown model type: {model_type}")
        
        state_dict = torch.load(checkpoint_path, map_location=device, weights_only=True)
        model.load_state_dict(state_dict)
        model = model.to(device)
        model.eval()
        
        _models[cache_key] = model
        # print(f"Loaded {model_type.upper()} model from {checkpoint_path}")
    
    return _models[cache_key]


def extract_embeddings(model, smiles_list: List[str], device: str = 'cpu'):
    """Extract molecular embeddings for a list of SMILES strings."""
    graphs = []
    valid_indices = []
    invalid_smiles = []
    
    for i, smiles in enumerate(smiles_list):
        graph = smiles_to_graph(smiles)
        if graph is not None:
            graphs.append(graph)
            valid_indices.append(i)
        else:
            invalid_smiles.append(smiles)
    
    if len(graphs) == 0:
        return np.array([], dtype=np.float32), valid_indices, invalid_smiles
    
    with torch.no_grad():
        batch = Batch.from_data_list(graphs).to(device)
        h, _ = model(batch)
        embeddings = h.cpu().numpy().astype(np.float32)
    
    return embeddings, valid_indices, invalid_smiles


def numpy_to_base64(arr: np.ndarray) -> str:
    """Encode a numpy array to base64 string."""
    buffer = BytesIO()
    np.save(buffer, arr)
    return base64.b64encode(buffer.getvalue()).decode('utf-8')


# FastAPI app
app = FastAPI(
    title="MolCLR Embedding API",
    description="Extract molecular embeddings from SMILES using pre-trained MolCLR models",
    version="1.0.0"
)


class EmbeddingRequest(BaseModel):
    smiles: List[str]
    model_type: str = "gin"
    
    class Config:
        json_schema_extra = {
            "example": {
                "smiles": ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"],
                "model_type": "gin"
            }
        }


class EmbeddingResponse(BaseModel):
    embeddings_base64: str
    shape: List[int]
    dtype: str
    valid_count: int
    invalid_smiles: List[str]


@app.get("/")
def root():
    """Health check endpoint."""
    return {"status": "ok", "message": "MolCLR Embedding API"}


@app.post("/api/embeddings", response_model=EmbeddingResponse)
def get_embeddings(request: EmbeddingRequest):
    """
    Extract molecular embeddings from SMILES strings.
    
    Returns embeddings as a base64-encoded numpy array.
    """
    if not request.smiles:
        raise HTTPException(status_code=400, detail="No SMILES provided")
    
    if request.model_type not in ["gin", "gcn"]:
        raise HTTPException(status_code=400, detail="model_type must be 'gin' or 'gcn'")
    
    try:
        model = get_model(request.model_type)
    except FileNotFoundError as e:
        raise HTTPException(status_code=500, detail=str(e))
    
    embeddings, valid_indices, invalid_smiles = extract_embeddings(
        model, request.smiles, device='cpu'
    )
    
    if len(embeddings) == 0:
        raise HTTPException(status_code=400, detail="No valid SMILES provided")
    
    return EmbeddingResponse(
        embeddings_base64=numpy_to_base64(embeddings),
        shape=list(embeddings.shape),
        dtype=str(embeddings.dtype),
        valid_count=len(valid_indices),
        invalid_smiles=invalid_smiles
    )


def ensure_models_ready():
    """Ensure model files are ready at startup. Downloads if needed and exits for restart."""
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    if not check_model_files_exist(base_dir):
        print("="*60)
        print("Model checkpoint files not found. Downloading...")
        print("="*60)
        try:
            prepare_model_files(base_dir)
            print("="*60)
            print("Model files downloaded successfully!")
            print("Please RESTART the server to load the models.")
            print("="*60)
            sys.exit(0)
        except Exception as e:
            print(f"ERROR: Failed to prepare model files: {e}")
            sys.exit(1)
    else:
        print("Model checkpoint files found. Ready to serve.")


if __name__ == "__main__":
    import uvicorn
    ensure_models_ready()
    uvicorn.run(app, host="0.0.0.0", port=3044)
