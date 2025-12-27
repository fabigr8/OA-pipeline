# Pocketeer → Docking/Cofolding Pipeline Scripts

This collection of scripts helps you extract useful information from Pocketeer's binding pocket detection and convert it into formats suitable for molecular docking and protein-ligand cofolding.

## Overview

**Pocketeer** detects binding pockets using alpha-spheres (similar to fpocket). These scripts convert that information into:
- Docking box coordinates for AutoDock Vina, Glide, etc.
- Pocket residue lists for cofolding constraints
- Visualization scripts for PyMOL
- JSON files with structured pocket data

## Installation

```bash
pip install pocketeer
pip install numpy  # Usually included with pocketeer
```

## Scripts

### 1. `pocketeer_to_docking.py` - Basic Extraction

**What it does:**
- Extracts pocket centers and sphere coordinates
- Calculates docking boxes with customizable padding
- Identifies pocket-lining residues
- Exports AutoDock Vina configs and PyMOL scripts

**Usage:**
```python
from pocketeer_to_docking import PocketExtractor

# Initialize
extractor = PocketExtractor("your_protein.pdb")

# Find pockets
pockets = extractor.find_pockets(r_min=3.0, r_max=6.0, min_spheres=30)

# Get pocket center (for pocket 0)
center = extractor.get_pocket_center(pocket_id=0)
print(f"Pocket center: {center}")

# Get docking box
box = extractor.get_docking_box(pocket_id=0, padding=5.0)
print(f"Box center: ({box['center_x']}, {box['center_y']}, {box['center_z']})")
print(f"Box size: ({box['size_x']}, {box['size_y']}, {box['size_z']})")

# Get pocket residues
residues = extractor.get_pocket_residues(pocket_id=0, distance_cutoff=5.0)

# Export Vina config
extractor.export_vina_config(pocket_id=0, output_file='vina_config.txt')

# Export PyMOL visualization
extractor.export_pymol_selection(pocket_id=0, output_file='pymol_pocket.pml')
```

### 2. `pocketeer_cofolding_guide.py` - Cofolding Integration

**What it does:**
- Classifies pockets as orthosteric or allosteric
- Creates spatial constraints for AlphaFold/RoseTTAFold
- Generates residue masks for focused cofolding
- Exports structured JSON for automated pipelines

**Usage:**
```python
from pocketeer_cofolding_guide import CofoldingPocketGuide

# Initialize
guide = CofoldingPocketGuide("your_protein.pdb")

# Detect pockets
guide.detect_pockets(r_min=3.0, r_max=6.0, min_spheres=30)

# Create constraints
# Option 1: Auto-classify (top pocket = orthosteric)
guide.create_constraints(distance_cutoff=5.0)

# Option 2: With known active site coordinates
known_site = [10.5, 20.3, 15.7]  # Your known coordinates
guide.create_constraints(distance_cutoff=5.0, known_active_site=known_site)

# Export for AlphaFold
guide.export_alphafold_constraints('af_constraints.json')

# Export strategy guide
guide.export_cofolding_plan('cofolding_plan.md')

# Get residue mask for focused cofolding
residue_mask = guide.create_residue_mask(pocket_id=0, expand_shells=2)
```

### 3. `pocketeer_pipeline.py` - Complete Automated Pipeline

**What it does:**
- Runs complete workflow from detection to output generation
- Creates configs for multiple docking tools
- Generates visualization scripts
- Exports everything in organized directory structure

**Command-line usage:**
```bash
# Basic usage
python pocketeer_pipeline.py --protein receptor.pdb

# With custom parameters
python pocketeer_pipeline.py \
    --protein receptor.pdb \
    --output my_results \
    --min-volume 250 \
    --top-n 3
```

**Python usage:**
```python
from pocketeer_pipeline import PocketeerPipeline

# Initialize
pipeline = PocketeerPipeline("receptor.pdb", output_dir="results")

# Run complete pipeline
pipeline.run_complete_pipeline(min_volume=200, top_n=5)
```

**Outputs:**
```
results/
├── pocket_summary.txt              # Human-readable summary
├── cofolding_constraints.json      # For AlphaFold/cofolding
├── visualize_pockets.pml           # PyMOL visualization
├── vina_configs/
│   ├── pocket_0_vina.txt          # Vina config for pocket 0
│   ├── pocket_1_vina.txt          # Vina config for pocket 1
│   └── run_all_vina.sh            # Batch script
└── glide_grids/
    ├── pocket_0_grid.in           # Glide grid for pocket 0
    └── pocket_1_grid.in           # Glide grid for pocket 1
```

## Use Cases

### Use Case 1: Find Best Docking Site

```python
from pocketeer_to_docking import PocketExtractor

extractor = PocketExtractor("kinase.pdb")
pockets = extractor.find_pockets()

# Get top pocket
box = extractor.get_docking_box(pocket_id=0, padding=5.0)

# Use in your docking workflow
print(f"Dock ligand at: ({box['center_x']:.2f}, {box['center_y']:.2f}, {box['center_z']:.2f})")
print(f"Search box size: ({box['size_x']:.1f}, {box['size_y']:.1f}, {box['size_z']:.1f})")
```

### Use Case 2: Separate Allosteric and Orthosteric Sites

```python
from pocketeer_cofolding_guide import CofoldingPocketGuide

guide = CofoldingPocketGuide("gpcr.pdb")
guide.detect_pockets()

# If you know the orthosteric site location
orthosteric_site = [15.2, 22.1, 18.5]  # From crystal structure
guide.create_constraints(known_active_site=orthosteric_site)

# Now pockets are classified
guide.print_summary()
# Output shows which are orthosteric vs allosteric

# Export separate configs
guide.export_alphafold_constraints('all_sites.json', include_allosteric=True)
```

### Use Case 3: Multi-Ligand Cofolding

```python
from pocketeer_pipeline import PocketeerPipeline

pipeline = PocketeerPipeline("protein.pdb")
pipeline.run_complete_pipeline(top_n=3)

# Now use cofolding_constraints.json:
import json
with open('pocketeer_output/cofolding_constraints.json') as f:
    constraints = json.load(f)

# Run cofolding separately for each pocket
for pocket in constraints['pockets']:
    if pocket['rank'] <= 2:  # Top 2 pockets
        center = pocket['center']
        print(f"Pocket {pocket['pocket_id']}: Center at ({center['x']}, {center['y']}, {center['z']})")
        # Use in your AlphaFold pipeline
```

## Key Features

### 1. Sphere Coordinates Access
All scripts provide access to the raw alpha-sphere coordinates from Pocketeer:

```python
# Access sphere data directly
pocket = pockets[0]
for sphere in pocket.spheres:
    print(f"Sphere center: {sphere.coord}, radius: {sphere.radius}")
```

### 2. Flexible Distance Cutoffs
Control how pocket residues are defined:

```python
# Tight definition (5Å from spheres)
residues_tight = extractor.get_pocket_residues(pocket_id=0, distance_cutoff=5.0)

# Loose definition (8Å from spheres)
residues_loose = extractor.get_pocket_residues(pocket_id=0, distance_cutoff=8.0)
```

### 3. Multiple Export Formats
- **Vina**: `.txt` config files
- **Glide**: `.in` grid files
- **PyMOL**: `.pml` scripts
- **JSON**: Structured data for custom pipelines
- **Markdown**: Human-readable strategies

## Integration Examples

### With AutoDock Vina

```bash
# 1. Generate config
python pocketeer_pipeline.py --protein receptor.pdb

# 2. Run Vina
cd pocketeer_output/vina_configs
bash run_all_vina.sh
```

### With AlphaFold

```python
import json
from alphafold import ...  # Your AlphaFold imports

# Load pocket constraints
with open('cofolding_constraints.json') as f:
    pockets = json.load(f)

# Use top pocket for ligand placement
top_pocket = pockets['pockets'][0]
ligand_center = top_pocket['center']
search_radius = top_pocket['radius']

# Feed to AlphaFold with spatial constraints
# ... your AlphaFold code ...
```

## Tips

1. **Volume filtering**: Pockets < 200 Å³ are often too small for drug-like molecules
2. **Top N selection**: Usually top 3-5 pockets by score are most relevant
3. **Padding**: Add 5-10 Å padding to docking boxes for flexibility
4. **Distance cutoff**: 5 Å is good for pocket residues; increase for broader regions
5. **Visualization**: Always visualize pockets in PyMOL before running expensive docking

## Troubleshooting

**Q: Pocketeer finds too many pockets**
- Increase `min_spheres` parameter (default: 30)
- Increase `min_volume` filter
- Keep only top N by score

**Q: Pocket coordinates seem wrong**
- Check your PDB file has correct atom positions
- Visualize with PyMOL to verify
- Try adjusting `r_min` and `r_max` parameters

**Q: Want to focus on specific region**
- Use `known_active_site` parameter in cofolding guide
- Filter pockets by distance to known site
- Manually select pocket by visual inspection

## References

- Pocketeer: https://github.com/cch1999/pocketeer
- fpocket algorithm: Le Guilloux et al. (2009) BMC Bioinformatics
- AutoDock Vina: Trott & Olson (2010) J. Comput. Chem.

## License

These scripts are provided as examples. Modify as needed for your workflow!
