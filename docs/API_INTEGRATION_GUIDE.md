# Batch Processing with Boltz-2 API Integration

## Overview

This script integrates NVIDIA Cloud Functions API calling directly into the batch processing pipeline. It:

1. Detects pockets with Pocketeer
2. Generates Boltz-2 configs
3. **Automatically calls Boltz-2 API**
4. Saves structures and confidence scores
5. Creates comprehensive summaries

## Setup

### Install Dependencies

```bash
pip install pocketeer httpx pandas numpy
```

### Set API Key

Option 1 - Environment variable (recommended):
```bash
export NVIDIA_API_KEY="nvapi-YOUR_KEY_HERE"
```

Option 2 - Command line:
```bash
python batch_pocketeer_boltz2_api.py --csv targets.csv --api-key "nvapi-YOUR_KEY_HERE"
```

## Usage

### Basic Usage

```bash
# Process all targets with API calls
python batch_pocketeer_boltz2_api.py --csv targets.csv
```

### Common Options

```bash
# Process first 3 targets (for testing)
python batch_pocketeer_boltz2_api.py --csv targets.csv --max-targets 3

# Custom results directory
python batch_pocketeer_boltz2_api.py --csv targets.csv --results-dir my_results

# Only generate configs, skip API calls
python batch_pocketeer_boltz2_api.py --csv targets.csv --skip-api

# Adjust distance cutoff
python batch_pocketeer_boltz2_api.py --csv targets.csv --distance-cutoff 6.0
```

## Output Structure

```
results/
├── BATCH_SUMMARY.csv
├── EGFR_1m17/
│   ├── TARGET_SUMMARY.md
│   ├── all_pockets.json
│   ├── pockets/
│   │   └── pocket_*.json
│   ├── boltz2_configs/
│   │   ├── orthosteric_pocket0.json
│   │   ├── allosteric_pocket1.json
│   │   └── dual_ligand.json
│   └── boltz2_results/              # NEW - API results
│       ├── orthosteric_result.json  # Full API response
│       ├── orthosteric_structure_0.mmcif  # Predicted structure
│       ├── allosteric_result.json
│       ├── allosteric_structure_0.mmcif
│       ├── dual_result.json
│       └── dual_structure_0.mmcif
└── ...
```

## API Call Details

### What Happens

For each target with SMILES data, the script:

1. **Generates configs** (orthosteric, allosteric, dual)
2. **Calls Boltz-2 API** for each config (except comparison configs)
3. **Handles long-polling** if request is queued (202 response)
4. **Saves results**:
   - Full JSON response: `{type}_result.json`
   - Structure files: `{type}_structure_0.mmcif`
5. **Updates summary** with confidence scores

### API Call Flow

```python
# 1. Load config
config = {
    "polymers": [...],
    "ligands": [...],
    "constraints": [pocket_constraint]  # From Pocketeer!
}

# 2. Call API
result = await make_boltz2_api_call(config, api_key)

# 3. Handle response
{
    "structures": [
        {
            "structure": "data_model\n_entry.id...",  # mmCIF content
            "format": "mmcif"
        }
    ],
    "confidence_scores": [0.87],
    "metrics": {...}
}

# 4. Save structure
with open('orthosteric_structure_0.mmcif', 'w') as f:
    f.write(result['structures'][0]['structure'])
```

### Which Configs Get Called

By default, the script calls Boltz-2 for:
- ✅ **Orthosteric** configs
- ✅ **Allosteric** configs
- ✅ **Dual ligand** configs
- ❌ **Comparison** configs (skipped to save API calls)

To include comparison configs, modify the script:
```python
# In _call_boltz2_for_configs method, remove this check:
if config_type == 'comparison':
    continue  # Remove or comment out
```

## Example Workflow

### 1. Prepare CSV

```csv
geneName,PDB Filename,PDB Fasta,OS SMILES,AS SMILES
EGFR,structures/1m17.pdb,MRPSGTA...,Nc1ncnc...,CC(C)(C)c1...
```

### 2. Run Pipeline

```bash
# Set API key
export NVIDIA_API_KEY="nvapi-YOUR_KEY"

# Run pipeline
python batch_pocketeer_boltz2_api.py \
    --csv targets.csv \
    --max-targets 1  # Start with 1 for testing
```

### 3. Monitor Progress

The script logs:
```
INFO - Processing: EGFR (1m17)
INFO - Step 1: Detecting pockets with Pocketeer...
INFO - Found 5 pockets
INFO - Step 2: Saving pocket information...
INFO - Step 3: Generating Boltz2 configurations...
INFO - Step 4: Calling Boltz-2 API...
INFO - Calling Boltz-2 for orthosteric_pocket0.json...
INFO - Task queued with ID: abc123, polling for results...
INFO - Task completed successfully
INFO - Saved result to orthosteric_result.json
INFO - Saved structure to orthosteric_structure_0.mmcif
INFO - ✓ Boltz-2 call successful for orthosteric
INFO - Confidence scores: [0.87]
```

### 4. Check Results

```bash
# View target summary
cat results/EGFR_1m17/TARGET_SUMMARY.md

# Check confidence scores
cat results/EGFR_1m17/boltz2_results/orthosteric_result.json | jq '.confidence_scores'

# View structure in PyMOL
pymol results/EGFR_1m17/boltz2_results/orthosteric_structure_0.mmcif
```

## Advanced Usage

### Programmatic Use

```python
import asyncio
from batch_pocketeer_boltz2_api import BatchPocketeerBoltz2WithAPI

async def run_pipeline():
    # Initialize
    batch = BatchPocketeerBoltz2WithAPI(
        csv_path='targets.csv',
        results_base_dir='results',
        api_key='nvapi-YOUR_KEY'
    )
    
    # Process all
    await batch.process_all(
        distance_cutoff=5.0,
        max_targets=None,
        call_boltz2=True
    )
    
    # Access results
    for result in batch.results:
        if result['status'] == 'success':
            print(f"{result['gene_name']}: {result['n_boltz2_results']} API calls")

asyncio.run(run_pipeline())
```

### Process Single Target

```python
import asyncio
from batch_pocketeer_boltz2_api import BatchPocketeerBoltz2WithAPI

async def process_one():
    batch = BatchPocketeerBoltz2WithAPI('targets.csv')
    
    # Get first row
    row = batch.targets_df.iloc[0]
    
    # Process it
    result = await batch.process_target(row, call_boltz2=True)
    
    print(f"Status: {result['status']}")
    print(f"Pockets: {result['n_pockets']}")
    print(f"API calls: {result['n_boltz2_results']}")

asyncio.run(process_one())
```

### Extract All Confidence Scores

```python
import json
from pathlib import Path

# Collect all confidence scores
scores = []

for target_dir in Path('results').glob('*/'):
    if target_dir.name.startswith('BATCH'):
        continue
    
    results_dir = target_dir / 'boltz2_results'
    if not results_dir.exists():
        continue
    
    for result_file in results_dir.glob('*_result.json'):
        with open(result_file) as f:
            data = json.load(f)
        
        config_type = result_file.stem.replace('_result', '')
        
        scores.append({
            'target': target_dir.name,
            'config_type': config_type,
            'confidence': data['confidence_scores'][0] if data['confidence_scores'] else None
        })

# Create DataFrame
import pandas as pd
df = pd.DataFrame(scores)
print(df)
df.to_csv('all_confidence_scores.csv', index=False)
```

## Rate Limiting & Performance

### API Calls Per Target

- Orthosteric only: **1 API call**
- Allosteric only: **1 API call**
- Both (dual): **3 API calls** (ortho + allo + dual)

### Typical Timing

- Pocket detection: ~5-30 seconds
- Config generation: <1 second
- Boltz-2 API call: 30-300 seconds (depends on queue)
- Total per target: ~1-10 minutes

### Optimization Tips

1. **Test first**:
   ```bash
   # Process 1 target to verify setup
   --max-targets 1
   ```

2. **Skip dual configs** if not needed:
   Comment out dual config generation in `_generate_boltz2_configs()`

3. **Generate configs first, call API later**:
   ```bash
   # Step 1: Generate all configs
   python batch_pocketeer_boltz2_api.py --csv targets.csv --skip-api
   
   # Step 2: Review configs, then call API
   python batch_pocketeer_boltz2_api.py --csv targets.csv
   ```

4. **Process in batches**:
   ```bash
   # Batch 1 (first 10)
   --max-targets 10
   
   # Batch 2 (next 10)
   # Modify CSV to skip processed ones
   ```

## Troubleshooting

### API Key Not Found

**Error**: "No API key provided"

**Solution**:
```bash
export NVIDIA_API_KEY="nvapi-YOUR_KEY"
# Or
--api-key "nvapi-YOUR_KEY"
```

### API Call Failed

**Error**: "API call failed with status 401"

**Solution**: Check your API key is valid and has correct permissions

### Task Timeout

**Error**: Request times out after 400 seconds

**Solution**: This is normal for queued requests. The script will keep polling. If it truly hangs, increase timeout:
```python
# In make_boltz2_api_call function
manual_timeout_seconds=600  # Increase from 400
```

### Too Many API Calls

**Problem**: Batch job makes hundreds of API calls

**Solution**:
```bash
# Only generate configs, review first
--skip-api

# Then run on subset
--max-targets 5
```

### No Structures Generated

**Problem**: API call succeeds but no structures saved

**Solution**: Check the result JSON:
```bash
cat results/TARGET/boltz2_results/orthosteric_result.json | jq '.structures'
```

If empty, check Boltz-2 API response for errors.

## Comparison: Old vs New

### Old Workflow (Manual API)

```bash
# 1. Generate configs
python batch_pocketeer_boltz2.py --csv targets.csv

# 2. Manually call API for each
for config in results/*/boltz2_configs/*.json; do
    curl -X POST ... -d @$config
done

# 3. Manually save results
# ... custom scripting ...
```

### New Workflow (Integrated API)

```bash
# 1. Everything in one command
python batch_pocketeer_boltz2_api.py --csv targets.csv

# Results automatically saved!
```

## Best Practices

1. **Always test with --max-targets 1 first**
2. **Set NVIDIA_API_KEY as environment variable**
3. **Review generated configs before large batches**
4. **Monitor logs for errors**
5. **Check confidence scores in summaries**
6. **Visualize structures in PyMOL/ChimeraX**
7. **Keep API key secure** (don't commit to git)

## Summary

✅ **Fully integrated API calling**  
✅ **Automatic structure saving**  
✅ **Long-polling support**  
✅ **Confidence score tracking**  
✅ **Async for efficiency**  
✅ **Comprehensive error handling**

Your pipeline is now completely automated from PDB → Pockets → Boltz-2 → Structures!
