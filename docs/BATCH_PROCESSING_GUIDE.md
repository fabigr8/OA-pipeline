# Batch Processing Guide for Pocketeer → Boltz-2

## Overview

This guide shows how to process multiple protein targets from a CSV file, with organized outputs for each target.

## Directory Structure

After batch processing, you'll have:

```
results/
├── BATCH_SUMMARY.csv              # Summary of all targets
├── BATCH_SUMMARY.md               # Markdown summary
├── EGFR_1m17/                     # Target 1 (geneName_pdbID)
│   ├── TARGET_SUMMARY.md
│   ├── all_pockets.json
│   ├── pockets/
│   │   ├── pocket_0.json
│   │   ├── pocket_1.json
│   │   └── ...
│   └── boltz2_configs/
│       ├── orthosteric_pocket0.json
│       ├── allosteric_pocket1.json
│       ├── dual_ligand.json
│       └── comparison_pocket0.json
├── CDK2_1hck/                     # Target 2
│   └── ...
└── BRAF_3og7/                     # Target 3
    └── ...
```

## CSV Format

Your input CSV must have these columns:

| Column | Required | Description | Example |
|--------|----------|-------------|---------|
| geneName | Yes | Gene name | EGFR |
| PDB Filename | Yes | Path to PDB file | structures/1m17.pdb |
| uniprotID | No | UniProt ID | P00533 |
| OS PDB ID | No | Reference PDB for orthosteric site | 1M17 |
| OS Lig ID | No | Orthosteric ligand ID | AQ4 |
| AS PDB ID | No | Reference PDB for allosteric site | 3GT8 |
| AS Lig ID | No | Allosteric ligand ID | 03P |
| OS Comment | No | Notes on orthosteric site | ATP-competitive inhibitor |
| AS Comment | No | Notes on allosteric site | Allosteric antibody site |
| PDB Fasta | Yes | Protein sequence for Boltz2 | MRPSGTA... |
| OS SMILES | No | Orthosteric ligand SMILES | Nc1ncnc2... |
| AS SMILES | No | Allosteric ligand SMILES | CC(C)(C)c1... |

**See `targets_template.csv` for an example!**

## Quick Start

### 1. Prepare Your CSV

```bash
# Use the template
cp targets_template.csv my_targets.csv

# Edit with your data
# Make sure PDB files exist at specified paths
```

### 2. Run Batch Processing

```bash
# Process all targets
python batch_pocketeer_boltz2.py --csv my_targets.csv

# With custom results directory
python batch_pocketeer_boltz2.py --csv my_targets.csv --results-dir my_results

# Process only first 5 targets (for testing)
python batch_pocketeer_boltz2.py --csv my_targets.csv --max-targets 5

# Adjust distance cutoff
python batch_pocketeer_boltz2.py --csv my_targets.csv --distance-cutoff 6.0
```

### 3. Review Results

```bash
# Check batch summary
cat results/BATCH_SUMMARY.md

# Look at specific target
cd results/EGFR_1m17
cat TARGET_SUMMARY.md

# View Boltz2 config
cat boltz2_configs/dual_ligand.json
```

## Python API Usage

```python
from batch_pocketeer_boltz2 import BatchPocketeerBoltz2

# Initialize
batch = BatchPocketeerBoltz2(
    csv_path='my_targets.csv',
    results_base_dir='results'
)

# Process all targets
batch.process_all(distance_cutoff=5.0)

# Or process single target
row = batch.targets_df.iloc[0]
result = batch.process_target(row)

# Access results
for result in batch.results:
    if result['status'] == 'success':
        print(f"{result['gene_name']}: {result['n_pockets']} pockets")
```

## What Gets Generated

For each target, the pipeline generates:

### 1. Pocket Information
- `all_pockets.json` - All detected pockets
- `pockets/pocket_N.json` - Individual pocket details with residues

### 2. Boltz2 Configurations

Based on available SMILES in your CSV:

| File | When Created | Purpose |
|------|--------------|---------|
| `orthosteric_pocket0.json` | If OS SMILES provided | Orthosteric ligand in pocket 0 |
| `allosteric_pocket1.json` | If AS SMILES provided | Allosteric ligand in pocket 1 |
| `dual_ligand.json` | If both SMILES provided | Both ligands simultaneously |
| `comparison_pocket0-2.json` | If OS SMILES provided | Test ligand in multiple pockets |

### 3. Summaries
- `TARGET_SUMMARY.md` - Human-readable overview
- Target-specific metadata from CSV

## Common Workflows

### Workflow 1: Screen Multiple Ligands

```csv
geneName,PDB Filename,PDB Fasta,OS SMILES,AS SMILES
EGFR,1m17.pdb,MRPSG...,SMILES1,
EGFR,1m17.pdb,MRPSG...,SMILES2,
EGFR,1m17.pdb,MRPSG...,SMILES3,
```

Creates separate directories:
- `EGFR_1m17/` (first ligand)
- `EGFR_1m17_1/` (second ligand) 
- `EGFR_1m17_2/` (third ligand)

**Better approach**: Use same row, send different Boltz2 configs from `comparison_pocket*.json`

### Workflow 2: Allosteric Drug Discovery

```csv
geneName,PDB Filename,PDB Fasta,OS SMILES,AS SMILES
CDK2,1hck.pdb,MENFQ...,STU_SMILES,CANDIDATE1
CDK2,1hck.pdb,MENFQ...,STU_SMILES,CANDIDATE2
```

Each creates `dual_ligand.json` with:
- Known orthosteric (STU) in pocket 0
- Candidate allosteric in pocket 1

### Workflow 3: Pocket Validation

```csv
geneName,PDB Filename,PDB Fasta,OS SMILES,OS PDB ID,OS Lig ID
BRAF,3og7.pdb,MAALS...,VEMURAFENIB,3OG7,AZ4
```

Pipeline detects pockets, then you:
1. Compare pocket 0 location to known site (3OG7)
2. Validate pocket detection worked
3. Use for similar proteins without known ligand

## Sending to Boltz-2 NIM

After batch processing, send configs to Boltz-2:

### Option 1: Command Line

```bash
#!/bin/bash
# Script to send all configs to Boltz-2

ENDPOINT="https://your-nim-endpoint/biology/mit/boltz2/predict"
API_KEY="your-api-key"

for config in results/*/boltz2_configs/orthosteric_pocket0.json; do
    target_dir=$(dirname $(dirname $config))
    target_name=$(basename $target_dir)
    
    echo "Processing $target_name..."
    
    curl -X POST "$ENDPOINT" \
         -H "Content-Type: application/json" \
         -H "Authorization: Bearer $API_KEY" \
         -d @"$config" \
         -o "$target_dir/boltz2_result.json"
done
```

### Option 2: Python

```python
import requests
import json
from pathlib import Path

def send_to_boltz2(config_file, endpoint, api_key):
    """Send config to Boltz-2 and save result."""
    with open(config_file) as f:
        config = json.load(f)
    
    response = requests.post(
        f"{endpoint}/biology/mit/boltz2/predict",
        json=config,
        headers={"Authorization": f"Bearer {api_key}"}
    )
    
    return response.json()

# Process all targets
results_dir = Path('results')
endpoint = "https://your-nim-endpoint"
api_key = "your-key"

for target_dir in results_dir.glob("*/"):
    if target_dir.name.startswith('BATCH'):
        continue
    
    # Send orthosteric config
    config_file = target_dir / 'boltz2_configs/orthosteric_pocket0.json'
    if config_file.exists():
        print(f"Processing {target_dir.name}...")
        result = send_to_boltz2(config_file, endpoint, api_key)
        
        # Save result
        output = target_dir / 'boltz2_orthosteric_result.json'
        with open(output, 'w') as f:
            json.dump(result, f, indent=2)
        
        # Extract structure
        if 'structures' in result:
            structure = result['structures'][0]['structure']
            cif_file = target_dir / 'predicted_orthosteric.cif'
            with open(cif_file, 'w') as f:
                f.write(structure)
```

## Analyzing Results

### Compare Pocket Detection

```python
import pandas as pd
import json
from pathlib import Path

results = []

for target_dir in Path('results').glob("*/"):
    if target_dir.name.startswith('BATCH'):
        continue
    
    # Load pocket info
    pocket_file = target_dir / 'all_pockets.json'
    if not pocket_file.exists():
        continue
    
    with open(pocket_file) as f:
        pockets = json.load(f)
    
    # Get top pocket
    top_pocket = pockets[0]
    
    results.append({
        'target': target_dir.name,
        'n_pockets': len(pockets),
        'top_score': top_pocket['score'],
        'top_volume': top_pocket['volume'],
        'top_residues': top_pocket['n_residues']
    })

# Create summary dataframe
df = pd.DataFrame(results)
print(df)
df.to_csv('pocket_comparison.csv', index=False)
```

### Extract Boltz-2 Confidence Scores

```python
import json
from pathlib import Path

for target_dir in Path('results').glob("*/"):
    result_file = target_dir / 'boltz2_orthosteric_result.json'
    if result_file.exists():
        with open(result_file) as f:
            result = json.load(f)
        
        confidence = result['confidence_scores'][0]
        print(f"{target_dir.name}: {confidence:.3f}")
```

## Tips & Best Practices

### 1. Organize Your Input Files

```
project/
├── structures/          # All PDB files here
│   ├── 1m17.pdb
│   ├── 1hck.pdb
│   └── 3og7.pdb
├── targets.csv          # Reference as: structures/1m17.pdb
└── results/             # Output goes here
```

### 2. Test First

```bash
# Process just first target to verify setup
python batch_pocketeer_boltz2.py --csv targets.csv --max-targets 1
```

### 3. Handle Missing SMILES

If you have PDB but no SMILES yet:
- Leave OS SMILES / AS SMILES empty
- Pipeline still detects pockets
- Generate SMILES later, add to CSV, re-run

### 4. Naming Convention

The pipeline creates folders as `{geneName}_{pdb_stem}`:
- `EGFR_1m17`
- `CDK2_1hck` 
- `BRAF_3og7`

Use clear, unique gene names to avoid collisions.

### 5. Incremental Processing

```python
# Load existing results
batch = BatchPocketeerBoltz2('targets.csv', 'results')

# Only process targets without results
for idx, row in batch.targets_df.iterrows():
    target_dir = batch.get_target_dir(row)
    
    # Skip if already processed
    if (target_dir / 'all_pockets.json').exists():
        print(f"Skipping {row['geneName']} (already processed)")
        continue
    
    batch.process_target(row)
```

## Troubleshooting

### Error: PDB file not found

**Problem**: CSV has wrong path

**Solution**: 
```bash
# Check paths in CSV
cat targets.csv | cut -d',' -f2

# Verify files exist
ls structures/*.pdb
```

### Error: Invalid SMILES

**Problem**: Malformed SMILES string

**Solution**: Validate SMILES before adding to CSV:
```python
from rdkit import Chem

smiles = "CC(C)..."
mol = Chem.MolFromSmiles(smiles)
if mol is None:
    print("Invalid SMILES!")
```

### No pockets detected

**Problem**: Pocketeer finds no pockets

**Solution**: Check PDB file quality, adjust parameters:
```python
# In batch_pocketeer_boltz2.py, modify:
pockets = pt.find_pockets(
    atomarray,
    r_min=2.5,      # Lower threshold
    r_max=8.0,      # Higher threshold  
    min_spheres=20  # Fewer required
)
```

### Too many output files

**Problem**: Lots of comparison configs generated

**Solution**: Modify `_generate_boltz2_configs()` to limit:
```python
# Generate comparison configs for top 2 pockets only
for i in range(min(2, len(pockets))):  # Changed from 3
```

## Summary

✅ **CSV-based batch processing**  
✅ **Organized per-target directories**  
✅ **Automatic Boltz-2 config generation**  
✅ **Support for orthosteric + allosteric**  
✅ **Comprehensive summaries**

Start with the template CSV, add your targets, and run the batch script!
