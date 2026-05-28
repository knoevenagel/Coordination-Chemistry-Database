# ChemDB API

---

## Metal Affinity API

**Port**: 3045

Calculate metal affinity scores for a molecule using k-NN and two-hop voting on ligand embeddings.

### POST /api/affinity

**Request**:
```json
{"smiles": "c1ccccc1", "k": 10}
```

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `smiles` | string | Yes | - | Molecule SMILES |
| `k` | int | No | 10 | k-NN neighbors |

**Response**:
```json
{
  "smiles": "c1ccccc1",
  "affinities": [
    {"metal": "Fe", "score": 1.0},
    {"metal": "Cu", "score": 0.8333}
  ],
  "k": 10,
  "neighbors": ["D249209940188324", "D007446817777545"],
  "elapsed_ms": 15.23
}
```

| Field | Description |
|-------|-------------|
| `affinities` | Metal scores (0-1), sorted descending |
| `neighbors` | DIDs of k nearest ligands |
| `elapsed_ms` | Processing time (ms) |

**Example**:
```bash
curl -X POST http://localhost:3045/api/affinity \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CCO", "k": 10}'
```

**Start Server**:
```bash
python src/affinity_api.py
```
