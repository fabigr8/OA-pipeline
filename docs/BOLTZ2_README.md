# Pocketeer → Boltz-2 NIM Integration Guide

## Perfect Match! 🎯

Your pipeline idea works **perfectly** with Boltz-2 NIM! Here's why:

### Boltz-2 Pocket Constraints

Boltz-2 has a native **`pocket` constraint type** that accepts exactly what Pocketeer provides:

```json
{
  "constraint_type": "pocket",
  "binder": "L",  // Ligand chain ID
  "contacts": [
    {"id": "A", "residue_index": 10},
    {"id": "A", "residue_index": 14},
    // More residues defining the pocket...
  ]
}
```

### What Pocketeer Provides

Pocketeer detects binding pockets and identifies:
- Pocket-lining residues (chain + residue number) ✅
- Pocket center coordinates ✅
- Pocket volume and score ✅
- Alpha-sphere locations ✅

**This maps perfectly to Boltz-2's contacts!**

---

## Quick Start

### 1. Basic Usage

```python
from pocketeer_boltz2 import Boltz2PocketConverter

# Detect pockets
converter = Boltz2PocketConverter("protein.pdb")
converter.detect_pockets()

# Create Boltz-2 request
request = converter.create_boltz2_request(
    protein_sequence="MVLSPADKTN...",
    ligand_smiles="CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    pocket_id=0  # Use top-ranked pocket
)

# Save for Boltz-2
import json
with open('boltz2_request.json', 'w') as f:
    json.dump(request, f, indent=2)
```

### 2. Allosteric vs Orthosteric

```python
from boltz2_workflow import AllostericOrthostericPipeline

# Initialize
pipeline = AllostericOrthostericPipeline("protein.pdb", "MVLSPAD...")

# Detect and classify
classification = pipeline.detect_and_classify_pockets()

# Create dual-ligand request
request = pipeline.create_dual_ligand_request(
    orthosteric_smiles="ATP_ANALOG_SMILES",
    allosteric_smiles="MODULATOR_SMILES",
    orthosteric_pocket_id=0,
    allosteric_pocket_id=1
)
```

### 3. Compare Multiple Pockets

```python
# Test ligand in different pockets
requests = pipeline.create_comparison_requests(
    ligand_smiles="YOUR_LIGAND_SMILES",
    pocket_ids=[0, 1, 2]  # Test top 3 pockets
)

# Send each to Boltz-2, compare confidence scores
# to see which pocket the ligand prefers
```

---

## Why This Works So Well

| Pocketeer Output | Boltz-2 Input | Match |
|-----------------|---------------|-------|
| Pocket residues (chain, resid) | `contacts: [{id, residue_index}]` | ✅ Perfect |
| Multiple pockets detected | Multiple constraint support | ✅ Perfect |
| Pocket ranking by score | Choose best pocket | ✅ Perfect |
| Allosteric site detection | Separate constraints | ✅ Perfect |

---

## Complete Workflow

### Step 1: Detect Pockets

```bash
python -c "
from pocketeer_boltz2 import Boltz2PocketConverter
conv = Boltz2PocketConverter('protein.pdb')
conv.detect_pockets()
conv.print_constraint_summary(pocket_id=0)
"
```

### Step 2: Generate Boltz-2 Configs

```python
converter.export_boltz2_configs(
    protein_sequence="YOUR_SEQUENCE",
    ligand_smiles="YOUR_SMILES",
    output_dir="boltz2_configs",
    max_pockets=3
)
```

This creates:
```
boltz2_configs/
├── pocket_0_boltz2.json
├── pocket_1_boltz2.json
├── pocket_2_boltz2.json
└── README.txt
```

### Step 3: Send to Boltz-2 NIM

```bash
curl -X POST https://YOUR_NIM_ENDPOINT/biology/mit/boltz2/predict \
     -H "Content-Type: application/json" \
     -H "Authorization: Bearer YOUR_API_KEY" \
     -d @boltz2_configs/pocket_0_boltz2.json
```

Or in Python:

```python
import requests

with open('boltz2_configs/pocket_0_boltz2.json') as f:
    request_data = json.load(f)

response = requests.post(
    "https://YOUR_NIM_ENDPOINT/biology/mit/boltz2/predict",
    json=request_data,
    headers={"Authorization": "Bearer YOUR_API_KEY"}
)

result = response.json()
```

### Step 4: Analyze Results

```python
# Boltz-2 returns confidence scores
confidence = result['confidence_scores'][0]
structure = result['structures'][0]['structure']

print(f"Confidence: {confidence}")

# Save predicted structure
with open('predicted_complex.cif', 'w') as f:
    f.write(structure)
```

---

## Use Cases

### Use Case 1: Find Best Binding Site

**Scenario**: You have a ligand and want to find where it binds

```python
converter = Boltz2PocketConverter("protein.pdb")
converter.detect_pockets()

# Test ligand in top 3 pockets
requests = []
for i in range(3):
    req = converter.create_boltz2_request(
        protein_sequence=seq,
        ligand_smiles=ligand,
        pocket_id=i
    )
    requests.append((i, req))

# Send all to Boltz-2, pick highest confidence
```

### Use Case 2: Orthosteric Inhibitor Design

**Scenario**: Design competitive inhibitor for active site

```python
# Assume pocket 0 is active site
request = converter.create_boltz2_request(
    protein_sequence=seq,
    ligand_smiles=inhibitor_smiles,
    pocket_id=0,
    distance_cutoff=5.0,  # Tight constraint
    sampling_steps=100,   # More sampling
    diffusion_samples=3   # Better predictions
)
```

### Use Case 3: Allosteric Modulator + Orthosteric Ligand

**Scenario**: Model both sites simultaneously

```python
pipeline = AllostericOrthostericPipeline("protein.pdb", seq)
pipeline.detect_and_classify_pockets()

# Both ligands at once
request = pipeline.create_dual_ligand_request(
    orthosteric_smiles=substrate_analog,
    allosteric_smiles=modulator,
    orthosteric_pocket_id=0,
    allosteric_pocket_id=1
)
```

### Use Case 4: GPCR Orthosteric vs Allosteric

**Scenario**: GPCR with known orthosteric site coordinates

```python
# Known orthosteric site from crystal structure
ortho_coords = [15.2, 22.1, 18.5]

classification = pipeline.detect_and_classify_pockets(
    known_orthosteric_site=ortho_coords
)

# Pockets within 15Å = orthosteric
# Pockets far away = allosteric
```

---

## Key Parameters

### Distance Cutoff

Controls how pocket residues are defined:

```python
# Tight binding site (direct contacts only)
converter.create_pocket_constraint(pocket_id=0, distance_cutoff=4.0)

# Loose binding site (include nearby residues)
converter.create_pocket_constraint(pocket_id=0, distance_cutoff=7.0)
```

**Recommendation**: Start with 5.0 Å

### Max Contacts

Limit number of constraining residues:

```python
contacts = converter.get_pocket_contacts(
    pocket_id=0,
    distance_cutoff=5.0,
    max_contacts=20  # Only use 20 closest residues
)
```

**Recommendation**: None (use all) or 30-50 for large pockets

### Boltz-2 Sampling

```python
request = converter.create_boltz2_request(
    ...,
    recycling_steps=4,     # More = better but slower (1-6)
    sampling_steps=100,    # More = better sampling (10-1000)
    diffusion_samples=3    # More = more predictions (1-5)
)
```

**For testing**: `recycling_steps=3, sampling_steps=50, diffusion_samples=1`  
**For production**: `recycling_steps=4, sampling_steps=100, diffusion_samples=3`

---

## Comparison: Without vs With Pocketeer

### Without Pocketeer
```python
# Manual blind docking - no guidance
request = {
    "polymers": [{"id": "A", "molecule_type": "protein", "sequence": seq}],
    "ligands": [{"id": "L", "smiles": ligand}],
    # No constraints - Boltz-2 searches entire protein
}
```
❌ Slow, may miss binding site  
❌ No control over binding location  
❌ Can't distinguish allosteric/orthosteric  

### With Pocketeer
```python
# Guided by detected pocket
request = {
    "polymers": [{"id": "A", "molecule_type": "protein", "sequence": seq}],
    "ligands": [{"id": "L", "smiles": ligand}],
    "constraints": [pocket_constraint]  # From Pocketeer!
}
```
✅ Fast - focused search  
✅ Guided to real binding sites  
✅ Can target specific pockets  
✅ Separates allosteric from orthosteric  

---

## Advanced: Multi-Ligand Systems

```python
# Example: Enzyme with substrate + cofactor + allosteric inhibitor
request = {
    "polymers": [
        {"id": "E", "molecule_type": "protein", "sequence": enzyme_seq}
    ],
    "ligands": [
        {"id": "S", "smiles": substrate_smiles},
        {"id": "C", "smiles": cofactor_smiles},
        {"id": "I", "smiles": inhibitor_smiles}
    ],
    "constraints": [
        # Substrate in active site (pocket 0)
        converter.create_pocket_constraint(0, ligand_chain_id="S"),
        # Cofactor in cofactor pocket (pocket 1)
        converter.create_pocket_constraint(1, ligand_chain_id="C"),
        # Inhibitor in allosteric site (pocket 2)
        converter.create_pocket_constraint(2, ligand_chain_id="I")
    ]
}
```

---

## Troubleshooting

### Too Many Contacts

**Problem**: Pocket has 200+ residues, making constraint too restrictive

**Solution**:
```python
contacts = converter.get_pocket_contacts(
    pocket_id=0,
    distance_cutoff=4.0,  # Tighter cutoff
    max_contacts=50       # Limit contacts
)
```

### No Pockets Found

**Problem**: Pocketeer doesn't find any pockets

**Solution**:
```python
# Relax parameters
pockets = converter.detect_pockets(
    r_min=2.5,    # Smaller minimum
    r_max=8.0,    # Larger maximum
    min_spheres=20  # Fewer spheres needed
)
```

### Wrong Pocket Selected

**Problem**: Top pocket isn't the actual binding site

**Solution**:
```python
# Visualize all pockets first
for i in range(len(pockets)):
    converter.print_constraint_summary(pocket_id=i)

# Choose manually
request = converter.create_boltz2_request(..., pocket_id=2)  # Use pocket 2
```

---

## Files in This Package

1. **`pocketeer_boltz2.py`** - Core converter class
   - `Boltz2PocketConverter` - Main class
   - Converts Pocketeer → Boltz-2 format

2. **`boltz2_workflow.py`** - Complete workflow
   - `AllostericOrthostericPipeline` - Full pipeline
   - Handles classification and multi-ligand scenarios

3. **Original scripts** (also work with Boltz-2):
   - `pocketeer_to_docking.py` - General extraction
   - `pocketeer_cofolding_guide.py` - Cofolding tools
   - `pocketeer_pipeline.py` - Automated pipeline

---

## Example Output

### Boltz-2 Request JSON
```json
{
  "polymers": [
    {
      "id": "A",
      "molecule_type": "protein",
      "sequence": "MVLSPADKTNVKAAW..."
    }
  ],
  "ligands": [
    {
      "id": "L",
      "smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
    }
  ],
  "constraints": [
    {
      "constraint_type": "pocket",
      "binder": "L",
      "contacts": [
        {"id": "A", "residue_index": 10},
        {"id": "A", "residue_index": 14},
        {"id": "A", "residue_index": 18},
        ...
      ]
    }
  ],
  "recycling_steps": 3,
  "sampling_steps": 50,
  "diffusion_samples": 1,
  "output_format": "mmcif"
}
```

### Boltz-2 Response
```json
{
  "structures": [
    {
      "structure": "data_model\n_entry.id model\n...",
      "format": "mmcif",
      "name": "prediction_0"
    }
  ],
  "confidence_scores": [0.87],
  "metrics": {
    "runtime_ms": 45231
  }
}
```

---

## Summary

✅ **Yes, your pipeline works perfectly with Boltz-2!**

The Pocketeer → Boltz-2 integration provides:
1. **Automatic pocket detection** from protein structure
2. **Direct conversion** to Boltz-2 pocket constraints
3. **Allosteric/orthosteric classification**
4. **Multi-ligand support** for complex scenarios
5. **Focused, guided docking** vs blind search

**Next steps**: Try the example scripts with your protein and ligands!
