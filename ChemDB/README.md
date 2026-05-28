# ChemDB

A chemical database system for coordination complex analysis and ligand processing.

> **Documentation**: See the [reengineering report](doc/report.md) for detailed technical analysis.

## Overview

### Layer Definition

- M: Metal
- L1: Cordination Complex Space
- L2: Ligand Complexation Form
- L3: Ligand Standalone Form
- L4: Coordinative Fragment (Critical Subgraph)
- L5: Irreducible Ligand (Simplest Ligand-Metal Interaction)

### Processing Pipeline

![Processing Steps](doc/step_diagram.drawio.svg)

**Pipeline Steps:**
- **Step 1**: Initial data processing
- **Step 2**: Ligand repair
- **Step 4/5**: Generate IRL and GA
- **Step 6/7**: Ligand marking and L4 fragment generation
- **Step 8**: Generate CSV files for Neo4j bulk import
- **Step 9**: Generate SMILES fingerprint
- **Step 10**: Merge layers and generate combined database _(in progress)_
- **Step 11**: Generate embeddings of L3 and M-L3 relationships
- **Step 12**: Build precomputed indices from Step8 CSV outputs (no Neo4j required)
- **Step 13**: Batch candidate building with CSV output

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python src/main.py
```

## Project Structure

### Input Data
```
data/pubchem/*.csv           # Original PubChem data
data/metal_list.txt          # Metal list
data/p_elements_list.txt     # P elements
```

### Output Data
```
tmp/*.csv                    # Generated CSV files
tmp/*.json                   # Generation statistics
```

### Source Code
```
src/step*.py                 # Pipeline step scripts
src/proxy.py                 # Proxy server (UI/Tanimoto/Query)
src/server.py                # Query server
src/molclr_api.py            # MolCLR API server
src/affinity_api.py          # Affinity API server
src/main.py                  # Pipeline orchestration
src/comparison.py            # Legacy MySQL comparison
src/*_utils.py               # Utility modules
src/tools/*.py               # Migrated legacy tools
```

### Tools
- **tools/chemtool.py** - Calculate DID from SMILES
- **tools/neo4j/** - Neo4j Docker setup
- **tools/query.sh** - Query local CSV layer data
- **tools/printSummary.sh** - Print comparison summary

## Services & Ports

| Port | Service |
|------|----------|
| 3040 | Unified proxy |
| &nbsp;&nbsp;3041 | &nbsp;&nbsp;Query server |
| &nbsp;&nbsp;3042 | &nbsp;&nbsp;Tanimoto similarity server |
| &nbsp;&nbsp;3043 | &nbsp;&nbsp;UI |
| &nbsp;&nbsp;3044 | &nbsp;&nbsp;MolCLR API |
| &nbsp;&nbsp;3045 | &nbsp;&nbsp;Affinity API |
| 3074 | Neo4j browser |
| 3087 | Neo4j database |

## Step 12 & 13: Candidate Building Pipeline

Step 12 and Step 13 provide a **Neo4j-free** candidate building pipeline that works directly with Step 8 CSV outputs.

### Step 12: Build Precomputed Indices

Builds precomputed indices from Step 8 CSV outputs for efficient candidate searching.

**Outputs:**
- `l3_l5.json` - L3 → L5 mapping
- `l5_l3.json` - L5 → L3 inverted index
- `l5_freq_weight.json` - L5 frequency and IDF weights
- `l3_gac.json` - L3 GAC values
- `l5_l3_filtered_K{K}.json` - GAC-filtered L5→L3 index
- `m_l3_pairs.csv` - Metal-ligand pairs for candidate building

**Usage:**
```bash
cd src
python step12.py --input-dir ../tmp --output-dir ../tmp --K 30
```

### Step 13: Batch Candidate Building

Builds candidate sets for ligands and outputs directly to CSV format.

**Outputs:**
- `step13_targets.csv` - Target ligand metadata (DID, L5 sets, candidate count)
- `step13_candidates.csv` - Candidate relationships with distances

**Usage:**
```bash
cd src
python step13.py --input-dir ../tmp --output-dir ../tmp --num-samples 1000 --D 1000 --K 30
```

**Parameters:**
- `--num-samples`: Number of targets to process (-1 for all)
- `--D`: Rare/common L5 threshold (default: 1000)
- `--K`: GAC cutoff value (default: 30)
- `--seed`: Random seed for reproducibility (default: 42)
- `--no-distances`: Skip distance calculation for faster processing

## Neo4j Setup

1. **Start Neo4j container:**
   ```bash
   cd tools/neo4j
   docker-compose up -d
   ```

2. **Import CSV data:**
   ```bash
   cd tools/neo4j
   ./import.sh
   ```
