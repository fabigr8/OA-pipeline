"""
Complete Pipeline: Pocketeer → Docking/Cofolding

This script demonstrates a complete workflow for:
1. Detecting pockets with Pocketeer
2. Ranking and filtering pockets
3. Preparing inputs for different docking tools
4. Setting up cofolding with AlphaFold or similar

Usage:
    python pocketeer_pipeline.py --protein receptor.pdb --ligand ligand.sdf
"""

import argparse
import json
from pathlib import Path
import numpy as np
import pocketeer as pt


class PocketeerPipeline:
    """Complete pipeline from pocket detection to docking/cofolding prep."""
    
    def __init__(self, protein_path, output_dir='pocketeer_output'):
        self.protein_path = Path(protein_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Load protein
        print(f"Loading protein: {protein_path}")
        self.atomarray = pt.load_structure(str(protein_path))
        self.pockets = None
        
    def run_pocket_detection(self, min_volume=200, top_n=5):
        """
        Run pocket detection and filter results.
        
        Args:
            min_volume: Minimum pocket volume (Å³)
            top_n: Maximum number of pockets to keep
        """
        print("\n" + "="*70)
        print("STEP 1: POCKET DETECTION")
        print("="*70)
        
        # Run Pocketeer
        all_pockets = pt.find_pockets(
            self.atomarray,
            r_min=3.0,
            r_max=6.0,
            min_spheres=30,
            merge_distance=2.5,
            sasa_threshold=25.0
        )
        
        print(f"\nInitial pockets found: {len(all_pockets)}")
        
        # Filter by volume
        self.pockets = [p for p in all_pockets if p.volume >= min_volume]
        print(f"After volume filter (≥{min_volume} Å³): {len(self.pockets)}")
        
        # Keep top N by score
        self.pockets = sorted(self.pockets, 
                             key=lambda p: p.score, 
                             reverse=True)[:top_n]
        print(f"Top {top_n} pockets by score: {len(self.pockets)}")
        
        # Save summary
        self._save_pocket_summary()
        
        return self.pockets
    
    def _save_pocket_summary(self):
        """Save pocket summary to file."""
        summary_file = self.output_dir / 'pocket_summary.txt'
        
        with open(summary_file, 'w') as f:
            f.write("POCKET DETECTION SUMMARY\n")
            f.write("=" * 70 + "\n\n")
            
            for i, pocket in enumerate(self.pockets):
                sphere_coords = np.array([s.coord for s in pocket.spheres])
                center = np.mean(sphere_coords, axis=0)
                
                f.write(f"Pocket {i}:\n")
                f.write(f"  Rank:        {i + 1}\n")
                f.write(f"  Score:       {pocket.score:.3f}\n")
                f.write(f"  Volume:      {pocket.volume:.1f} Å³\n")
                f.write(f"  Spheres:     {pocket.n_spheres}\n")
                f.write(f"  Center (x,y,z): {center[0]:8.3f} {center[1]:8.3f} {center[2]:8.3f}\n")
                f.write("\n")
        
        print(f"Summary saved: {summary_file}")
    
    def prepare_vina_docking(self, padding=5.0):
        """
        Prepare AutoDock Vina configuration files for each pocket.
        
        Args:
            padding: Extra space around pocket (Å)
        """
        print("\n" + "="*70)
        print("STEP 2: PREPARE VINA DOCKING CONFIGS")
        print("="*70)
        
        vina_dir = self.output_dir / 'vina_configs'
        vina_dir.mkdir(exist_ok=True)
        
        for i, pocket in enumerate(self.pockets):
            # Calculate bounding box
            sphere_coords = np.array([s.coord for s in pocket.spheres])
            min_coords = np.min(sphere_coords, axis=0) - padding
            max_coords = np.max(sphere_coords, axis=0) + padding
            center = (min_coords + max_coords) / 2
            size = max_coords - min_coords
            
            # Write config file
            config_file = vina_dir / f'pocket_{i}_vina.txt'
            
            with open(config_file, 'w') as f:
                f.write(f"# AutoDock Vina config for Pocket {i}\n")
                f.write(f"# Score: {pocket.score:.3f}, Volume: {pocket.volume:.1f} Å³\n\n")
                f.write(f"receptor = receptor.pdbqt\n")
                f.write(f"ligand = ligand.pdbqt\n\n")
                f.write(f"center_x = {center[0]:.3f}\n")
                f.write(f"center_y = {center[1]:.3f}\n")
                f.write(f"center_z = {center[2]:.3f}\n\n")
                f.write(f"size_x = {size[0]:.1f}\n")
                f.write(f"size_y = {size[1]:.1f}\n")
                f.write(f"size_z = {size[2]:.1f}\n\n")
                f.write(f"exhaustiveness = 8\n")
                f.write(f"num_modes = 10\n")
            
            print(f"  Created: {config_file.name}")
        
        # Create batch script
        batch_file = vina_dir / 'run_all_vina.sh'
        with open(batch_file, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# Batch script to run Vina docking on all pockets\n\n")
            for i in range(len(self.pockets)):
                f.write(f"vina --config pocket_{i}_vina.txt --out pocket_{i}_result.pdbqt\n")
        
        batch_file.chmod(0o755)
        print(f"\n  Batch script: {batch_file.name}")
    
    def prepare_glide_grids(self, inner_box_ratio=0.5):
        """
        Prepare Schrödinger Glide grid generation files.
        
        Args:
            inner_box_ratio: Ratio of inner box to outer box
        """
        print("\n" + "="*70)
        print("STEP 3: PREPARE GLIDE GRID FILES")
        print("="*70)
        
        glide_dir = self.output_dir / 'glide_grids'
        glide_dir.mkdir(exist_ok=True)
        
        for i, pocket in enumerate(self.pockets):
            sphere_coords = np.array([s.coord for s in pocket.spheres])
            center = np.mean(sphere_coords, axis=0)
            
            # Calculate box size
            distances = np.linalg.norm(sphere_coords - center, axis=0)
            radius = np.max(distances) + 5.0
            outer_size = radius * 2
            inner_size = outer_size * inner_box_ratio
            
            # Write Glide grid input file
            grid_file = glide_dir / f'pocket_{i}_grid.in'
            
            with open(grid_file, 'w') as f:
                f.write(f"# Glide grid for Pocket {i}\n")
                f.write(f"GRID_CENTER {center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}\n")
                f.write(f"GRIDFILE pocket_{i}_grid.zip\n")
                f.write(f"INNERBOX {inner_size:.1f}, {inner_size:.1f}, {inner_size:.1f}\n")
                f.write(f"OUTERBOX {outer_size:.1f}, {outer_size:.1f}, {outer_size:.1f}\n")
                f.write(f"RECEP_FILE receptor.mae\n")
            
            print(f"  Created: {grid_file.name}")
    
    def prepare_cofolding_json(self, distance_cutoff=5.0):
        """
        Prepare JSON file with pocket information for cofolding.
        
        Args:
            distance_cutoff: Distance to define pocket residues (Å)
        """
        print("\n" + "="*70)
        print("STEP 4: PREPARE COFOLDING CONSTRAINTS")
        print("="*70)
        
        cofolding_data = {
            'protein': str(self.protein_path),
            'n_pockets': len(self.pockets),
            'pockets': []
        }
        
        for i, pocket in enumerate(self.pockets):
            # Calculate pocket center and radius
            sphere_coords = np.array([s.coord for s in pocket.spheres])
            center = np.mean(sphere_coords, axis=0)
            distances = np.linalg.norm(sphere_coords - center, axis=0)
            radius = np.max(distances)
            
            # Get pocket residues
            residues = []
            atom_coords = self.atomarray.coord
            
            for sphere_center in sphere_coords:
                dists = np.linalg.norm(atom_coords - sphere_center, axis=1)
                nearby = np.where(dists <= distance_cutoff)[0]
                
                for idx in nearby:
                    res_info = {
                        'chain': self.atomarray.chain_id[idx],
                        'resid': int(self.atomarray.res_id[idx]),
                        'resname': self.atomarray.res_name[idx]
                    }
                    if res_info not in residues:
                        residues.append(res_info)
            
            # Sort residues
            residues = sorted(residues, key=lambda x: (x['chain'], x['resid']))
            
            pocket_data = {
                'pocket_id': i,
                'rank': i + 1,
                'score': float(pocket.score),
                'volume': float(pocket.volume),
                'n_spheres': pocket.n_spheres,
                'center': {
                    'x': float(center[0]),
                    'y': float(center[1]),
                    'z': float(center[2])
                },
                'radius': float(radius),
                'search_box': {
                    'center': center.tolist(),
                    'size': float(radius * 2 + 10)
                },
                'residues': residues,
                'n_residues': len(residues),
                'cofolding_strategy': self._suggest_strategy(i, pocket)
            }
            
            cofolding_data['pockets'].append(pocket_data)
        
        # Save to JSON
        json_file = self.output_dir / 'cofolding_constraints.json'
        with open(json_file, 'w') as f:
            json.dump(cofolding_data, f, indent=2)
        
        print(f"  Created: {json_file.name}")
        print(f"  Total pockets: {len(self.pockets)}")
        print(f"  Pocket residues: {[len(p['residues']) for p in cofolding_data['pockets']]}")
    
    def _suggest_strategy(self, pocket_id, pocket):
        """Suggest cofolding strategy for a pocket."""
        if pocket_id == 0:
            return {
                'type': 'orthosteric',
                'priority': 'high',
                'recommendation': 'Primary target for ligand placement. Use as main constraint.'
            }
        elif pocket.volume > 500:
            return {
                'type': 'potential_allosteric',
                'priority': 'medium',
                'recommendation': 'Large pocket. Could be allosteric site or secondary binding site.'
            }
        else:
            return {
                'type': 'potential_allosteric',
                'priority': 'low',
                'recommendation': 'Smaller pocket. Consider for specialized ligands or allosteric modulators.'
            }
    
    def create_pymol_visualization(self):
        """Create PyMOL script to visualize all pockets."""
        print("\n" + "="*70)
        print("STEP 5: CREATE VISUALIZATION SCRIPT")
        print("="*70)
        
        pymol_file = self.output_dir / 'visualize_pockets.pml'
        
        with open(pymol_file, 'w') as f:
            f.write("# PyMOL script to visualize detected pockets\n\n")
            f.write(f"load {self.protein_path}\n")
            f.write("hide everything\n")
            f.write("show cartoon\n")
            f.write("color gray80\n\n")
            
            colors = ['red', 'blue', 'green', 'yellow', 'orange', 'purple']
            
            for i, pocket in enumerate(self.pockets):
                sphere_coords = np.array([s.coord for s in pocket.spheres])
                center = np.mean(sphere_coords, axis=0)
                color = colors[i % len(colors)]
                
                f.write(f"# Pocket {i} (Score: {pocket.score:.2f})\n")
                f.write(f"pseudoatom pocket_{i}_center, pos=[{center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}]\n")
                f.write(f"show spheres, pocket_{i}_center\n")
                f.write(f"color {color}, pocket_{i}_center\n")
                f.write(f"set sphere_scale, 2.0, pocket_{i}_center\n\n")
            
            f.write("# Zoom to show all pockets\n")
            f.write("zoom\n")
        
        print(f"  Created: {pymol_file.name}")
    
    def run_complete_pipeline(self, min_volume=200, top_n=5):
        """Run the complete pipeline."""
        print("\n" + "="*80)
        print(" POCKETEER → DOCKING/COFOLDING PIPELINE")
        print("="*80)
        
        # Step 1: Detect pockets
        self.run_pocket_detection(min_volume=min_volume, top_n=top_n)
        
        # Step 2: Prepare Vina configs
        self.prepare_vina_docking(padding=5.0)
        
        # Step 3: Prepare Glide grids
        self.prepare_glide_grids()
        
        # Step 4: Prepare cofolding constraints
        self.prepare_cofolding_json(distance_cutoff=5.0)
        
        # Step 5: Create visualization
        self.create_pymol_visualization()
        
        # Final summary
        print("\n" + "="*80)
        print(" PIPELINE COMPLETE!")
        print("="*80)
        print(f"\nOutput directory: {self.output_dir}")
        print("\nGenerated files:")
        print(f"  ├── pocket_summary.txt")
        print(f"  ├── cofolding_constraints.json")
        print(f"  ├── visualize_pockets.pml")
        print(f"  ├── vina_configs/")
        print(f"  │   ├── pocket_0_vina.txt")
        print(f"  │   ├── pocket_1_vina.txt")
        print(f"  │   └── run_all_vina.sh")
        print(f"  └── glide_grids/")
        print(f"      ├── pocket_0_grid.in")
        print(f"      └── pocket_1_grid.in")
        print("\nNext steps:")
        print("  1. Review pocket_summary.txt to select best pockets")
        print("  2. For docking: Use configs in vina_configs/ or glide_grids/")
        print("  3. For cofolding: Use constraints in cofolding_constraints.json")
        print("  4. Visualize in PyMOL: pymol visualize_pockets.pml")
        print()


def main():
    """Command-line interface."""
    parser = argparse.ArgumentParser(
        description='Run Pocketeer pocket detection and prepare docking/cofolding inputs'
    )
    parser.add_argument(
        '--protein',
        type=str,
        required=True,
        help='Input protein structure (PDB or CIF)'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='pocketeer_output',
        help='Output directory (default: pocketeer_output)'
    )
    parser.add_argument(
        '--min-volume',
        type=float,
        default=200.0,
        help='Minimum pocket volume in Å³ (default: 200)'
    )
    parser.add_argument(
        '--top-n',
        type=int,
        default=5,
        help='Maximum number of top pockets to keep (default: 5)'
    )
    
    args = parser.parse_args()
    
    # Run pipeline
    pipeline = PocketeerPipeline(args.protein, args.output)
    pipeline.run_complete_pipeline(
        min_volume=args.min_volume,
        top_n=args.top_n
    )


if __name__ == '__main__':
    main()
