"""
Batch Pocketeer → Boltz-2 Pipeline from CSV

Processes multiple targets from a CSV file with proper organization.
Each target gets its own subdirectory in the results folder.

CSV Format:
- geneName: Gene name
- PDB Filename: Path to PDB file for Pocketeer
- uniprotID: UniProt identifier
- OS PDB ID: Orthosteric site PDB reference
- OS Lig ID: Orthosteric ligand ID
- AS PDB ID: Allosteric site PDB reference  
- AS Lig ID: Allosteric ligand ID
- OS Comment: Notes on orthosteric site
- AS Comment: Notes on allosteric site
- PDB Fasta: Protein sequence for Boltz2
- OS SMILES: Orthosteric ligand SMILES
- AS SMILES: Allosteric ligand SMILES
"""

import pandas as pd
import json
from pathlib import Path
import numpy as np
from datetime import datetime
from typing import Optional
import pocketeer as pt


class BatchPocketeerBoltz2:
    """Batch process multiple targets from CSV."""
    
    def __init__(self, csv_path, results_base_dir='results'):
        """
        Initialize batch processor.
        
        Args:
            csv_path: Path to CSV file with targets
            results_base_dir: Base directory for all results
        """
        self.csv_path = Path(csv_path)
        self.results_base = Path(results_base_dir)
        self.results_base.mkdir(exist_ok=True)
        
        # Load targets
        self.targets_df = pd.read_csv(csv_path)
        print(f"Loaded {len(self.targets_df)} targets from {csv_path}")
        
        # Results tracking
        self.results = []
        
    def get_target_dir(self, row):
        """
        Create organized directory name for target.
        
        Args:
            row: DataFrame row with target info
            
        Returns:
            Path: Directory path for this target
        """
        # Use gene name and PDB ID for folder name
        pdb_stem = Path(row['PDB Filename']).stem
        gene_name = row['geneName'].replace('/', '_').replace(' ', '_')
        
        # Format: geneName_pdbID (e.g., EGFR_1m17, CDK2_1hck)
        dir_name = f"{gene_name}_{pdb_stem}"
        
        target_dir = self.results_base / dir_name
        target_dir.mkdir(exist_ok=True)
        
        return target_dir
    
    def process_target(self, row, distance_cutoff=5.0):
        """
        Process a single target from the CSV.
        
        Args:
            row: DataFrame row with target data
            distance_cutoff: Distance for pocket contacts
            
        Returns:
            dict: Processing results
        """
        gene_name = row['geneName']
        pdb_file = row['PDB Filename']
        
        print(f"\n{'='*80}")
        print(f"Processing: {gene_name} ({Path(pdb_file).stem})")
        print(f"{'='*80}")
        
        # Create target directory
        target_dir = self.get_target_dir(row)
        print(f"Output directory: {target_dir}")
        
        # Check if PDB file exists
        if not Path(pdb_file).exists():
            error_msg = f"PDB file not found: {pdb_file}"
            print(f"ERROR: {error_msg}")
            return {
                'gene_name': gene_name,
                'pdb_file': pdb_file,
                'status': 'error',
                'error': error_msg
            }
        
        try:
            # Step 1: Run Pocketeer
            print("\nStep 1: Detecting pockets with Pocketeer...")
            atomarray = pt.load_structure(pdb_file)
            pockets = pt.find_pockets(
                atomarray,
                r_min=3.0,
                r_max=6.0,
                min_spheres=30,
                merge_distance=2.5
            )
            print(f"Found {len(pockets)} pockets")
            
            # Step 2: Save pocket information
            print("\nStep 2: Saving pocket information...")
            pocket_info = self._save_pocket_info(pockets, atomarray, target_dir, 
                                                 distance_cutoff)
            
            # Step 3: Generate Boltz2 configs
            print("\nStep 3: Generating Boltz2 configurations...")
            boltz_configs = self._generate_boltz2_configs(
                row, pockets, atomarray, target_dir, distance_cutoff
            )
            
            # Step 4: Create target summary
            print("\nStep 4: Creating target summary...")
            self._create_target_summary(row, pockets, target_dir, boltz_configs)
            
            result = {
                'gene_name': gene_name,
                'pdb_file': pdb_file,
                'target_dir': str(target_dir),
                'status': 'success',
                'n_pockets': len(pockets),
                'n_boltz_configs': len(boltz_configs),
                'top_pocket_score': float(pockets[0].score) if pockets else None,
                'top_pocket_volume': float(pockets[0].volume) if pockets else None
            }
            
            print(f"\n✓ Successfully processed {gene_name}")
            return result
            
        except Exception as e:
            error_msg = f"Error processing {gene_name}: {str(e)}"
            print(f"ERROR: {error_msg}")
            import traceback
            traceback.print_exc()
            
            return {
                'gene_name': gene_name,
                'pdb_file': pdb_file,
                'status': 'error',
                'error': error_msg
            }
    
    def _save_pocket_info(self, pockets, atomarray, target_dir, distance_cutoff):
        """Save detailed pocket information."""
        pockets_dir = target_dir / 'pockets'
        pockets_dir.mkdir(exist_ok=True)
        
        pocket_info = []
        
        for i, pocket in enumerate(pockets):
            sphere_coords = np.array([s.coord for s in pocket.spheres])
            center = np.mean(sphere_coords, axis=0)
            radius = np.max(np.linalg.norm(sphere_coords - center, axis=0))
            
            # Get pocket residues
            residues = self._get_pocket_residues(
                pocket, atomarray, distance_cutoff
            )
            
            info = {
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
                'residues': residues,
                'n_residues': len(residues)
            }
            
            pocket_info.append(info)
            
            # Save individual pocket JSON
            pocket_file = pockets_dir / f'pocket_{i}.json'
            with open(pocket_file, 'w') as f:
                json.dump(info, f, indent=2)
        
        # Save all pockets JSON
        all_pockets_file = target_dir / 'all_pockets.json'
        with open(all_pockets_file, 'w') as f:
            json.dump(pocket_info, f, indent=2)
        
        print(f"  Saved pocket info: {all_pockets_file}")
        
        return pocket_info
    
    def _get_pocket_residues(self, pocket, atomarray, distance_cutoff):
        """Extract residues near pocket."""
        sphere_coords = np.array([s.coord for s in pocket.spheres])
        atom_coords = atomarray.coord
        
        residues = set()
        
        for sphere_center in sphere_coords:
            distances = np.linalg.norm(atom_coords - sphere_center, axis=1)
            nearby = np.where(distances <= distance_cutoff)[0]
            
            for idx in nearby:
                residues.add((
                    atomarray.chain_id[idx],
                    int(atomarray.res_id[idx]),
                    atomarray.res_name[idx]
                ))
        
        # Convert to list of dicts and sort
        residues_list = [
            {'chain': r[0], 'resid': r[1], 'resname': r[2]}
            for r in sorted(residues, key=lambda x: (x[0], x[1]))
        ]
        
        return residues_list
    
    def _generate_boltz2_configs(self, row, pockets, atomarray, 
                                 target_dir, distance_cutoff):
        """Generate Boltz2 configuration files."""
        boltz_dir = target_dir / 'boltz2_configs'
        boltz_dir.mkdir(exist_ok=True)
        
        configs = []
        
        # Get sequences and SMILES from row
        protein_seq = row['PDB Fasta']
        os_smiles = row.get('OS SMILES', '')
        as_smiles = row.get('AS SMILES', '')
        
        # Generate configs based on available data
        if pd.notna(os_smiles) and os_smiles:
            # Orthosteric ligand config (pocket 0)
            config = self._create_boltz2_config(
                protein_seq, os_smiles, pockets[0], atomarray,
                ligand_chain_id='O', ligand_type='orthosteric',
                distance_cutoff=distance_cutoff
            )
            
            config_file = boltz_dir / 'orthosteric_pocket0.json'
            with open(config_file, 'w') as f:
                json.dump(config, f, indent=2)
            
            configs.append({
                'type': 'orthosteric',
                'pocket_id': 0,
                'file': str(config_file.name)
            })
            print(f"  Created: {config_file.name}")
        
        if pd.notna(as_smiles) and as_smiles and len(pockets) > 1:
            # Allosteric ligand config (pocket 1)
            config = self._create_boltz2_config(
                protein_seq, as_smiles, pockets[1], atomarray,
                ligand_chain_id='A', ligand_type='allosteric',
                distance_cutoff=distance_cutoff
            )
            
            config_file = boltz_dir / 'allosteric_pocket1.json'
            with open(config_file, 'w') as f:
                json.dump(config, f, indent=2)
            
            configs.append({
                'type': 'allosteric',
                'pocket_id': 1,
                'file': str(config_file.name)
            })
            print(f"  Created: {config_file.name}")
        
        # Dual ligand config if both available
        if (pd.notna(os_smiles) and os_smiles and 
            pd.notna(as_smiles) and as_smiles and len(pockets) > 1):
            
            dual_config = self._create_dual_ligand_config(
                protein_seq, os_smiles, as_smiles, 
                pockets[0], pockets[1], atomarray,
                distance_cutoff
            )
            
            config_file = boltz_dir / 'dual_ligand.json'
            with open(config_file, 'w') as f:
                json.dump(dual_config, f, indent=2)
            
            configs.append({
                'type': 'dual',
                'pocket_id': [0, 1],
                'file': str(config_file.name)
            })
            print(f"  Created: {config_file.name}")
        
        # Generate comparison configs (test ligand in multiple pockets)
        if pd.notna(os_smiles) and os_smiles:
            for i in range(min(3, len(pockets))):
                config = self._create_boltz2_config(
                    protein_seq, os_smiles, pockets[i], atomarray,
                    ligand_chain_id='L', ligand_type='comparison',
                    distance_cutoff=distance_cutoff
                )
                
                config_file = boltz_dir / f'comparison_pocket{i}.json'
                with open(config_file, 'w') as f:
                    json.dump(config, f, indent=2)
                
                configs.append({
                    'type': 'comparison',
                    'pocket_id': i,
                    'file': str(config_file.name)
                })
        
        return configs
    
    def _create_boltz2_config(self, protein_seq, ligand_smiles, pocket, 
                             atomarray, ligand_chain_id, ligand_type,
                             distance_cutoff):
        """Create single Boltz2 config."""
        # Get pocket contacts
        contacts = self._get_boltz2_contacts(pocket, atomarray, distance_cutoff)
        
        config = {
            "_metadata": {
                "ligand_type": ligand_type,
                "pocket_score": float(pocket.score),
                "pocket_volume": float(pocket.volume)
            },
            "polymers": [
                {
                    "id": "P",
                    "molecule_type": "protein",
                    "sequence": protein_seq
                }
            ],
            "ligands": [
                {
                    "id": ligand_chain_id,
                    "smiles": ligand_smiles
                }
            ],
            "constraints": [
                {
                    "constraint_type": "pocket",
                    "binder": ligand_chain_id,
                    "contacts": contacts
                }
            ],
            "recycling_steps": 3,
            "sampling_steps": 50,
            "diffusion_samples": 1,
            "output_format": "mmcif"
        }
        
        return config
    
    def _create_dual_ligand_config(self, protein_seq, os_smiles, as_smiles,
                                   os_pocket, as_pocket, atomarray, distance_cutoff):
        """Create dual ligand Boltz2 config."""
        os_contacts = self._get_boltz2_contacts(os_pocket, atomarray, distance_cutoff)
        as_contacts = self._get_boltz2_contacts(as_pocket, atomarray, distance_cutoff)
        
        config = {
            "_metadata": {
                "ligand_type": "dual",
                "orthosteric_score": float(os_pocket.score),
                "allosteric_score": float(as_pocket.score)
            },
            "polymers": [
                {
                    "id": "P",
                    "molecule_type": "protein",
                    "sequence": protein_seq
                }
            ],
            "ligands": [
                {
                    "id": "O",
                    "smiles": os_smiles
                },
                {
                    "id": "A",
                    "smiles": as_smiles
                }
            ],
            "constraints": [
                {
                    "constraint_type": "pocket",
                    "binder": "O",
                    "contacts": os_contacts
                },
                {
                    "constraint_type": "pocket",
                    "binder": "A",
                    "contacts": as_contacts
                }
            ],
            "recycling_steps": 4,
            "sampling_steps": 100,
            "diffusion_samples": 2,
            "output_format": "mmcif"
        }
        
        return config
    
    def _get_boltz2_contacts(self, pocket, atomarray, distance_cutoff):
        """Get pocket contacts in Boltz2 format."""
        sphere_coords = np.array([s.coord for s in pocket.spheres])
        atom_coords = atomarray.coord
        
        contacts = set()
        
        for sphere_center in sphere_coords:
            distances = np.linalg.norm(atom_coords - sphere_center, axis=1)
            nearby = np.where(distances <= distance_cutoff)[0]
            
            for idx in nearby:
                contacts.add((
                    atomarray.chain_id[idx],
                    int(atomarray.res_id[idx])
                ))
        
        # Convert to Boltz2 format
        boltz_contacts = [
            {"id": chain, "residue_index": resid}
            for chain, resid in sorted(contacts, key=lambda x: (x[0], x[1]))
        ]
        
        return boltz_contacts
    
    def _create_target_summary(self, row, pockets, target_dir, boltz_configs):
        """Create summary file for target."""
        summary_file = target_dir / 'TARGET_SUMMARY.md'
        
        lines = [
            f"# Target Summary: {row['geneName']}",
            "",
            "## Target Information",
            f"- **Gene Name**: {row['geneName']}",
            f"- **UniProt ID**: {row.get('uniprotID', 'N/A')}",
            f"- **PDB File**: {row['PDB Filename']}",
            "",
            "## Orthosteric Site",
            f"- **Reference PDB**: {row.get('OS PDB ID', 'N/A')}",
            f"- **Ligand ID**: {row.get('OS Lig ID', 'N/A')}",
            f"- **SMILES**: {row.get('OS SMILES', 'N/A')}",
            f"- **Comment**: {row.get('OS Comment', 'N/A')}",
            "",
            "## Allosteric Site",
            f"- **Reference PDB**: {row.get('AS PDB ID', 'N/A')}",
            f"- **Ligand ID**: {row.get('AS Lig ID', 'N/A')}",
            f"- **SMILES**: {row.get('AS SMILES', 'N/A')}",
            f"- **Comment**: {row.get('AS Comment', 'N/A')}",
            "",
            "## Detected Pockets",
            ""
        ]
        
        for i, pocket in enumerate(pockets[:5]):  # Top 5
            sphere_coords = np.array([s.coord for s in pocket.spheres])
            center = np.mean(sphere_coords, axis=0)
            
            lines.extend([
                f"### Pocket {i}",
                f"- Score: {pocket.score:.3f}",
                f"- Volume: {pocket.volume:.1f} Å³",
                f"- Center: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})",
                f"- Spheres: {pocket.n_spheres}",
                ""
            ])
        
        lines.extend([
            "## Generated Boltz2 Configs",
            ""
        ])
        
        for config in boltz_configs:
            lines.append(f"- **{config['type']}**: `boltz2_configs/{config['file']}`")
        
        lines.extend([
            "",
            "## Files",
            "```",
            f"{target_dir.name}/",
            "├── TARGET_SUMMARY.md          # This file",
            "├── all_pockets.json           # All pocket data",
            "├── pockets/                   # Individual pocket JSONs",
            "│   ├── pocket_0.json",
            "│   └── ...",
            "└── boltz2_configs/            # Boltz2 configuration files",
            "    ├── orthosteric_pocket0.json",
            "    ├── allosteric_pocket1.json",
            "    ├── dual_ligand.json",
            "    └── comparison_pocket*.json",
            "```",
            ""
        ])
        
        with open(summary_file, 'w') as f:
            f.write('\n'.join(lines))
        
        print(f"  Summary: {summary_file}")
    
    def process_all(self, distance_cutoff=5.0, max_targets=None):
        """
        Process all targets in the CSV.
        
        Args:
            distance_cutoff: Distance for pocket contacts
            max_targets: Optional limit on number to process
        """
        print(f"\n{'='*80}")
        print(f"BATCH PROCESSING: {len(self.targets_df)} TARGETS")
        print(f"{'='*80}")
        
        n_targets = min(len(self.targets_df), max_targets) if max_targets else len(self.targets_df)
        
        for idx, row in self.targets_df.iterrows():
            if max_targets and idx >= max_targets:
                break
            
            result = self.process_target(row, distance_cutoff)
            self.results.append(result)
        
        # Create batch summary
        self._create_batch_summary()
        
        print(f"\n{'='*80}")
        print("BATCH PROCESSING COMPLETE")
        print(f"{'='*80}")
        print(f"Processed: {len(self.results)} targets")
        print(f"Success: {sum(1 for r in self.results if r['status'] == 'success')}")
        print(f"Errors: {sum(1 for r in self.results if r['status'] == 'error')}")
        print(f"\nResults directory: {self.results_base}")
        print(f"Batch summary: {self.results_base / 'BATCH_SUMMARY.csv'}")
    
    def _create_batch_summary(self):
        """Create summary CSV of all results."""
        summary_df = pd.DataFrame(self.results)
        summary_file = self.results_base / 'BATCH_SUMMARY.csv'
        summary_df.to_csv(summary_file, index=False)
        
        print(f"\nBatch summary saved: {summary_file}")
        
        # Also create markdown summary
        md_file = self.results_base / 'BATCH_SUMMARY.md'
        lines = [
            "# Batch Processing Summary",
            "",
            f"**Date**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"**Total Targets**: {len(self.results)}",
            f"**Successful**: {sum(1 for r in self.results if r['status'] == 'success')}",
            f"**Errors**: {sum(1 for r in self.results if r['status'] == 'error')}",
            "",
            "## Results",
            ""
        ]
        
        for result in self.results:
            status_icon = "✓" if result['status'] == 'success' else "✗"
            lines.append(f"- {status_icon} **{result['gene_name']}**")
            if result['status'] == 'success':
                lines.append(f"  - Directory: `{Path(result['target_dir']).name}`")
                lines.append(f"  - Pockets: {result['n_pockets']}")
                lines.append(f"  - Boltz configs: {result['n_boltz_configs']}")
            else:
                lines.append(f"  - Error: {result.get('error', 'Unknown')}")
            lines.append("")
        
        with open(md_file, 'w') as f:
            f.write('\n'.join(lines))


def main():
    """Example usage."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Batch process targets from CSV for Pocketeer → Boltz2'
    )
    parser.add_argument(
        '--csv',
        type=str,
        required=True,
        help='Path to CSV file with targets'
    )
    parser.add_argument(
        '--results-dir',
        type=str,
        default='results',
        help='Base directory for results (default: results)'
    )
    parser.add_argument(
        '--distance-cutoff',
        type=float,
        default=5.0,
        help='Distance cutoff for pocket contacts in Å (default: 5.0)'
    )
    parser.add_argument(
        '--max-targets',
        type=int,
        default=None,
        help='Maximum number of targets to process (default: all)'
    )
    
    args = parser.parse_args()
    
    # Run batch processing
    batch = BatchPocketeerBoltz2(args.csv, args.results_dir)
    batch.process_all(
        distance_cutoff=args.distance_cutoff,
        max_targets=args.max_targets
    )


if __name__ == '__main__':
    main()
