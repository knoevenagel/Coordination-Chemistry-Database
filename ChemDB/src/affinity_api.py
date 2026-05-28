#!/usr/bin/env python3
"""
Metal Affinity Calculation API

Calculates metal affinity scores for input SMILES using:
1. MolCLR embedding extraction
2. k-NN search on L3 ligand embeddings
3. Two-hop voting through M-L3 relationships

Endpoint:
    POST /api/affinity
    
Request Body:
    {
        "smiles": "CCO",
        "k": 10  # optional, default 10
    }

Response:
    {
        "smiles": "CCO",
        "affinities": [
            {"metal": "Fe", "score": 0.85},
            {"metal": "Cu", "score": 0.72},
            ...
        ],
        "k": 10,
        "neighbors": ["D123...", "D456...", ...]  # DIDs of k nearest neighbors
    }

Usage:
    python affinity_api.py
    # or
    uvicorn affinity_api:app --host 0.0.0.0 --port 3045
"""

import base64
import logging
import os
import sys
import time
from collections import defaultdict
from io import BytesIO
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from molclr_api import get_model, extract_embeddings, numpy_to_base64

# Configuration
BASE_DIR = Path(__file__).parent.parent
INPUT_DIR = BASE_DIR / "tmp"
DEFAULT_K = 10
DEFAULT_MODEL_TYPE = "gin"
API_PORT = 3045

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def base64_to_numpy(b64_str: str) -> np.ndarray:
    """Decode a base64 string back to numpy array."""
    buffer = BytesIO(base64.b64decode(b64_str))
    return np.load(buffer)


def get_most_free_gpu() -> str:
    """Select the GPU with most free memory using nvidia-smi."""
    import subprocess
    import torch
    
    if not torch.cuda.is_available():
        return 'cpu'
    
    try:
        result = subprocess.run(
            ['nvidia-smi', '--query-gpu=index,memory.free', '--format=csv,noheader,nounits'],
            capture_output=True, text=True, check=True
        )
        
        max_free = 0
        best_gpu = 0
        
        for line in result.stdout.strip().split('\n'):
            parts = [p.strip() for p in line.split(',')]
            gpu_idx = int(parts[0])
            free_mb = int(parts[1])
            
            if free_mb > max_free:
                max_free = free_mb
                best_gpu = gpu_idx
        
        return f'cuda:{best_gpu}'
        
    except Exception:
        return 'cuda:0'


class AffinityCalculator:
    """Calculates metal affinity scores using k-NN and two-hop voting."""
    
    def __init__(self, embeddings_file: Path, m_l3_file: Path, model_type: str = DEFAULT_MODEL_TYPE):
        self.model_type = model_type
        self.device = get_most_free_gpu()
        
        # Statistics
        self.stats = {
            'warmup_time': 0,
            'model_load_time': 0,
            'embeddings_load_time': 0,
            'relationships_load_time': 0,
            'embeddings_file_size_mb': 0,
            'relationships_file_size_mb': 0,
            'embeddings_memory_mb': 0,
            'total_memory_mb': 0,
            'device': self.device,
            'model_type': model_type,
            'num_ligands': 0,
            'num_metals': 0,
            'num_relationships': 0,
            'request_count': 0,
            'total_request_time': 0,
            'avg_request_time': 0
        }
        
        warmup_start = time.time()
        logger.info(f"Using device: {self.device}")
        
        # Get file sizes
        self.stats['embeddings_file_size_mb'] = round(embeddings_file.stat().st_size / (1024 * 1024), 2)
        self.stats['relationships_file_size_mb'] = round(m_l3_file.stat().st_size / (1024 * 1024), 2)
        logger.info(f"Embeddings file size: {self.stats['embeddings_file_size_mb']} MB")
        logger.info(f"Relationships file size: {self.stats['relationships_file_size_mb']} MB")
        
        # Load model
        model_start = time.time()
        logger.info(f"Loading {model_type} model...")
        self.model = get_model(model_type, device=self.device)
        self.stats['model_load_time'] = round(time.time() - model_start, 3)
        logger.info(f"Model loaded in {self.stats['model_load_time']}s")
        
        # Load embeddings
        emb_start = time.time()
        logger.info(f"Loading embeddings from {embeddings_file}...")
        self._load_embeddings(embeddings_file)
        self.stats['embeddings_load_time'] = round(time.time() - emb_start, 3)
        logger.info(f"Embeddings loaded in {self.stats['embeddings_load_time']}s")
        
        # Load M-L3 relationships
        rel_start = time.time()
        logger.info(f"Loading M-L3 relationships from {m_l3_file}...")
        self._load_m_l3_relationships(m_l3_file)
        self.stats['relationships_load_time'] = round(time.time() - rel_start, 3)
        logger.info(f"Relationships loaded in {self.stats['relationships_load_time']}s")
        
        # Calculate memory usage
        self.stats['embeddings_memory_mb'] = round(self.embeddings.nbytes / (1024 * 1024), 2)
        self.stats['total_memory_mb'] = round(
            (self.embeddings.nbytes + self.embeddings_normalized.nbytes + 
             sys.getsizeof(self.dids) + sys.getsizeof(self.smiles_list) +
             sys.getsizeof(self.did_to_metals)) / (1024 * 1024), 2
        )
        
        self.stats['warmup_time'] = round(time.time() - warmup_start, 3)
        self.stats['num_ligands'] = len(self.dids)
        self.stats['num_metals'] = len(self.all_metals)
        
        logger.info(f"="*50)
        logger.info(f"AffinityCalculator warmup complete!")
        logger.info(f"  Total warmup time: {self.stats['warmup_time']}s")
        logger.info(f"  Ligands: {self.stats['num_ligands']:,}")
        logger.info(f"  Metals: {self.stats['num_metals']}")
        logger.info(f"  Relationships: {self.stats['num_relationships']:,}")
        logger.info(f"  Embeddings memory: {self.stats['embeddings_memory_mb']} MB")
        logger.info(f"  Total memory: {self.stats['total_memory_mb']} MB")
        logger.info(f"="*50)
    
    def _load_embeddings(self, embeddings_file: Path):
        """Load L3 embeddings into memory."""
        df = pd.read_csv(embeddings_file, dtype=str)
        
        self.dids = df['did'].tolist()
        self.smiles_list = df['smiles'].tolist()
        
        # Decode base64 embeddings to numpy array
        embeddings_list = []
        for emb_b64 in df['embedding']:
            emb = base64_to_numpy(emb_b64)
            embeddings_list.append(emb)
        
        self.embeddings = np.vstack(embeddings_list).astype(np.float32)
        
        # Normalize embeddings for cosine similarity
        norms = np.linalg.norm(self.embeddings, axis=1, keepdims=True)
        self.embeddings_normalized = self.embeddings / (norms + 1e-8)
        
        # Create DID to index mapping
        self.did_to_idx = {did: i for i, did in enumerate(self.dids)}
        
        logger.info(f"Loaded {len(self.dids)} embeddings, shape: {self.embeddings.shape}")
    
    def _load_m_l3_relationships(self, m_l3_file: Path):
        """Load M-L3 relationships and build lookup structures."""
        df = pd.read_csv(m_l3_file, dtype=str)
        
        # Build DID -> metals mapping
        self.did_to_metals: Dict[str, List[str]] = defaultdict(list)
        for _, row in df.iterrows():
            did = row['did']
            metal = row['symbol']
            self.did_to_metals[did].append(metal)
        
        # Get all unique metals
        self.all_metals = sorted(df['symbol'].unique().tolist())
        
        self.stats['num_relationships'] = len(df)
        logger.info(f"Loaded {len(df)} M-L3 relationships, {len(self.all_metals)} unique metals")
    
    def get_embedding(self, smiles: str) -> Optional[np.ndarray]:
        """Extract embedding for a SMILES string."""
        embeddings, valid_indices, invalid_smiles = extract_embeddings(
            self.model, [smiles], device=self.device
        )
        
        if len(embeddings) == 0:
            return None
        
        return embeddings[0].astype(np.float32)
    
    def find_k_nearest(self, query_embedding: np.ndarray, k: int) -> List[tuple]:
        """Find k nearest neighbors using cosine similarity.
        
        Returns:
            List of (did, smiles, similarity_score) tuples
        """
        # Normalize query
        query_norm = query_embedding / (np.linalg.norm(query_embedding) + 1e-8)
        
        # Compute cosine similarities
        similarities = np.dot(self.embeddings_normalized, query_norm)
        
        # Get top-k indices
        top_k_indices = np.argsort(similarities)[-k:][::-1]
        
        results = []
        for idx in top_k_indices:
            results.append((
                self.dids[idx],
                self.smiles_list[idx],
                float(similarities[idx])
            ))
        
        return results
    
    def calculate_affinity(self, smiles: str, k: int = DEFAULT_K) -> Dict:
        """Calculate metal affinity scores for a SMILES.
        
        Uses two-hop voting:
        1. Find k nearest L3 ligands
        2. For each neighbor, get its associated metals from M-L3
        3. Vote for metals weighted by similarity score
        4. Normalize scores to 0-1
        """
        # Get query embedding
        query_embedding = self.get_embedding(smiles)
        if query_embedding is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        # Find k nearest neighbors
        neighbors = self.find_k_nearest(query_embedding, k)
        
        # Two-hop voting: accumulate metal scores
        metal_scores: Dict[str, float] = defaultdict(float)
        
        for did, neighbor_smiles, similarity in neighbors:
            # Get metals associated with this ligand
            metals = self.did_to_metals.get(did, [])
            
            for metal in metals:
                # Vote weighted by similarity
                metal_scores[metal] += similarity
        
        # Normalize to 0-1
        if metal_scores:
            max_score = max(metal_scores.values())
            min_score = min(metal_scores.values())
            score_range = max_score - min_score
            
            if score_range > 0:
                normalized_scores = {
                    metal: (score - min_score) / score_range
                    for metal, score in metal_scores.items()
                }
            else:
                # All scores are the same
                normalized_scores = {metal: 1.0 for metal in metal_scores}
        else:
            normalized_scores = {}
        
        # Sort by score descending
        sorted_affinities = sorted(
            [{"metal": m, "score": round(s, 4)} for m, s in normalized_scores.items()],
            key=lambda x: x["score"],
            reverse=True
        )
        
        return {
            "smiles": smiles,
            "affinities": sorted_affinities,
            "k": k,
            "neighbors": [did for did, _, _ in neighbors]
        }
    
    def update_request_stats(self, elapsed_time: float):
        """Update request statistics."""
        self.stats['request_count'] += 1
        self.stats['total_request_time'] += elapsed_time
        self.stats['avg_request_time'] = round(
            self.stats['total_request_time'] / self.stats['request_count'], 4
        )


# Global calculator instance
calculator: Optional[AffinityCalculator] = None


# FastAPI app
app = FastAPI(
    title="Metal Affinity API",
    description="Calculate metal affinity scores for molecules using k-NN and two-hop voting",
    version="1.0.0"
)


class AffinityRequest(BaseModel):
    smiles: str
    k: int = DEFAULT_K
    
    class Config:
        json_schema_extra = {
            "example": {
                "smiles": "CCO",
                "k": 10
            }
        }


class AffinityItem(BaseModel):
    metal: str
    score: float


class AffinityResponse(BaseModel):
    smiles: str
    affinities: List[AffinityItem]
    k: int
    neighbors: List[str]
    elapsed_ms: float


@app.on_event("startup")
async def startup_event():
    """Load data on startup."""
    global calculator
    
    embeddings_file = INPUT_DIR / "l3_embeddings.csv"
    m_l3_file = INPUT_DIR / "m_l3_relationships.csv"
    
    if not embeddings_file.exists():
        logger.error(f"Embeddings file not found: {embeddings_file}")
        raise FileNotFoundError(f"Embeddings file not found: {embeddings_file}")
    
    if not m_l3_file.exists():
        logger.error(f"M-L3 file not found: {m_l3_file}")
        raise FileNotFoundError(f"M-L3 file not found: {m_l3_file}")
    
    calculator = AffinityCalculator(embeddings_file, m_l3_file)
    logger.info("Affinity calculator ready")


@app.get("/")
def root():
    """Health check endpoint."""
    return {
        "status": "ok",
        "message": "Metal Affinity API",
        "ligands": len(calculator.dids) if calculator else 0,
        "metals": len(calculator.all_metals) if calculator else 0
    }


@app.get("/api/stats")
def get_stats():
    """Get server statistics."""
    if calculator is None:
        raise HTTPException(status_code=503, detail="Calculator not initialized")
    
    return calculator.stats


@app.post("/api/affinity", response_model=AffinityResponse)
def get_affinity(request: AffinityRequest):
    """Calculate metal affinity scores for a SMILES string."""
    if calculator is None:
        raise HTTPException(status_code=503, detail="Calculator not initialized")
    
    if not request.smiles:
        raise HTTPException(status_code=400, detail="SMILES is required")
    
    if request.k < 1:
        raise HTTPException(status_code=400, detail="k must be >= 1")
    
    try:
        start_time = time.time()
        result = calculator.calculate_affinity(request.smiles, request.k)
        elapsed = time.time() - start_time
        elapsed_ms = round(elapsed * 1000, 2)
        
        # Update stats
        calculator.update_request_stats(elapsed)
        
        logger.info(f"Affinity calculated for {request.smiles[:30]}... in {elapsed_ms}ms")
        result['elapsed_ms'] = elapsed_ms
        return result
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Error calculating affinity: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/api/metals")
def list_metals():
    """List all available metals."""
    if calculator is None:
        raise HTTPException(status_code=503, detail="Calculator not initialized")
    
    return {"metals": calculator.all_metals}


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=API_PORT)
