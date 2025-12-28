"""
Batch Pocketeer → Boltz-2 Pipeline with Two Execution Modes

Mode 1: Orthosteric ligand only with pocket constraint
Mode 2: Both ligands (OS + AS) without pocket constraints

Reads FASTA from files and SMILES from CSV.
"""

import asyncio
import pandas as pd
import json
from pathlib import Path
import numpy as np
from datetime import datetime
from typing import Optional, Dict, Any, List
import logging
import sys
import os

# Check for required dependencies
missing_deps = []
try:
    import httpx
except ImportError:
    missing_deps.append("httpx")
try:
    import pocketeer as pt
except ImportError:
    missing_deps.append("pocketeer")

if missing_deps:
    print("Error: Missing required dependencies. Please install them using:")
    print(f"pip install {' '.join(missing_deps)}")
    sys.exit(1)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# NVIDIA Cloud Functions URLs
PUBLIC_URL = "https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"
STATUS_URL = "https://api.nvcf.nvidia.com/v2/nvcf/pexec/status/{task_id}"


async def make_boltz2_api_call(
    data: Dict[str, Any],
    api_key: str,
    nvcf_poll_seconds: int = 300,
    manual_timeout_seconds: int = 400
) -> Dict:
    """Make a call to Boltz-2 via NVIDIA Cloud Functions API."""
    async with httpx.AsyncClient() as client:
        headers = {
            "Authorization": f"Bearer {api_key}",
            "NVCF-POLL-SECONDS": f"{nvcf_poll_seconds}",
            "Content-Type": "application/json"
        }
        
        logger.debug(f"Making Boltz-2 API call to {PUBLIC_URL}")
        
        response = await client.post(
            PUBLIC_URL,
            json=data,
            headers=headers,
            timeout=manual_timeout_seconds
        )
        
        logger.debug(f"Response status: {response.status_code}")
        
        if response.status_code == 202:
            # Handle 202 Accepted - need to poll for results
            task_id = response.headers.get("nvcf-reqid")
            logger.info(f"Task queued with ID: {task_id}, polling for results...")
            
            while True:
                await asyncio.sleep(5)
                
                status_response = await client.get(
                    STATUS_URL.format(task_id=task_id),
                    headers=headers,
                    timeout=manual_timeout_seconds
                )
                
                if status_response.status_code == 200:
                    logger.info("Task completed successfully")
                    return status_response.json()
                elif status_response.status_code in [400, 401, 404, 422, 500]:
                    error_msg = f"Error while polling task: {status_response.text}"
                    logger.error(error_msg)
                    raise Exception(error_msg)
                    
        elif response.status_code == 200:
            logger.info("Request completed immediately")
            return response.json()
        else:
            error_msg = f"API call failed with status {response.status_code}: {response.text}"
            logger.error(error_msg)
            raise Exception(error_msg)


def read_fasta_file(fasta_path):
    """
    Read sequence from FASTA file.
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        str: Protein sequence
    """
    fasta_file = Path(fasta_path)
    
    if not fasta_file.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    
    sequence_lines = []
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            # Skip header lines
            if line.startswith('>'):
                continue
            # Collect sequence lines
            if line:
                sequence_lines.append(line)
    
    sequence = ''.join(sequence_lines)
    
    if not sequence:
        raise ValueError(f"No sequence found in FASTA file: {fasta_path}")
    
    return sequence


class BatchPocketeerBoltz2Pipeline:
    """
    Batch processing pipeline with two execution modes:
    Mode 1: OS ligand only with pocket constraint
    Mode 2: Both ligands without pocket constraints
    """
    
    def __init__(self, csv_path, results_base_dir='results', api_key=None):
        """
        Initialize batch processor.
        
        Args:
            csv_path: Path to CSV file with targets
            results_base_dir: Base directory for all results
            api_key: NVIDIA API key (or set NVIDIA_API_KEY env var)
        """
        self.csv_path = Path(csv_path)
        self.results_base = Path(results_base_dir)
        self.results_base.mkdir(exist_ok=True)
        
        # Get API key
        self.api_key = api_key or os.getenv('NVIDIA_API_KEY')
        if not self.api_key:
            logger.warning("No API key provided. Boltz-2 API calls will fail.")
            logger.warning("Set NVIDIA_API_KEY environment variable or pass api_key parameter")
        
        # Load targets
        self.targets_df = pd.read_csv(csv_path)
        logger.info(f"Loaded {len(self.targets_df)} targets from {csv_path}")
        
        # Results tracking
        self.results = []
        
    def get_target_dir(self, row):
        """Create organized directory name for target."""
        pdb_stem = Path(row['PDB Filename']).stem
        gene_name = row['geneName'].replace('/', '_').replace(' ', '_')
        
        dir_name = f"{gene_name}_{pdb_stem}"
        target_dir = self.results_base / dir_name
        target_dir.mkdir(exist_ok=True)
        
        return target_dir
    
    async def process_target(self, row, mode=1, distance_cutoff=5.0, call_boltz2=True):
        """
        Process a single target from the CSV.
        
        Args:
            row: DataFrame row with target data
            mode: Execution mode (1 or 2)
            distance_cutoff: Distance for pocket contacts
            call_boltz2: Whether to actually call Boltz-2 API
            
        Returns:
            dict: Processing results
        """
        gene_name = row['geneName']
        pdb_file = row['PDB Filename']
        
        logger.info(f"\n{'='*80}")
        logger.info(f"Processing: {gene_name} ({Path(pdb_file).stem}) - MODE {mode}")
        logger.info(f"{'='*80}")
        
        # Create target directory
        target_dir = self.get_target_dir(row)
        logger.info(f"Output directory: {target_dir}")
        
        # Check if PDB file exists
        if not Path(pdb_file).exists():
            error_msg = f"PDB file not found: {pdb_file}"
            logger.error(error_msg)
            return {
                'gene_name': gene_name,
                'pdb_file': pdb_file,
                'status': 'error',
                'error': error_msg
            }
        
        try:
            # Step 1: Run Pocketeer
            logger.info("Step 1: Detecting pockets with Pocketeer...")
            atomarray = pt.load_structure(pdb_file)
            pockets = pt.find_pockets(
                atomarray,
                r_min=3.0,
                r_max=6.0,
                min_spheres=30,
                merge_distance=2.5
            )
            logger.info(f"Found {len(pockets)} pockets")
            
            # Step 2: Save pocket information
            logger.info("Step 2: Saving pocket information...")
            pocket_info = self._save_pocket_info(pockets, atomarray, target_dir, 
                                                 distance_cutoff)
            
            # Step 3: Read FASTA sequence
            logger.info("Step 3: Reading FASTA sequence...")
            fasta_file = row['Fasta File']
            protein_sequence = read_fasta_file(fasta_file)
            logger.info(f"Loaded sequence: {len(protein_sequence)} residues")
            
            # Step 4: Generate Boltz2 configs based on mode
            logger.info(f"Step 4: Generating Boltz2 configurations (MODE {mode})...")
            boltz_configs = self._generate_boltz2_configs_by_mode(
                row, pockets, atomarray, target_dir, mode, 
                protein_sequence, distance_cutoff
            )
            
            # Step 5: Call Boltz-2 API if requested
            boltz2_results = []
            if call_boltz2 and self.api_key:
                logger.info("Step 5: Calling Boltz-2 API...")
                boltz2_results = await self._call_boltz2_for_configs(
                    target_dir, boltz_configs
                )
            else:
                logger.info("Step 5: Skipping Boltz-2 API calls")
            
            # Step 6: Create target summary
            logger.info("Step 6: Creating target summary...")
            self._create_target_summary(row, pockets, target_dir, boltz_configs, 
                                       boltz2_results, mode)
            
            result = {
                'gene_name': gene_name,
                'pdb_file': pdb_file,
                'target_dir': str(target_dir),
                'mode': mode,
                'status': 'success',
                'n_pockets': len(pockets),
                'n_boltz_configs': len(boltz_configs),
                'n_boltz2_results': len(boltz2_results),
                'top_pocket_score': float(pockets[0].score) if pockets else None,
                'top_pocket_volume': float(pockets[0].volume) if pockets else None
            }
            
            logger.info(f"✓ Successfully processed {gene_name}")
            return result
            
        except Exception as e:
            error_msg = f"Error processing {gene_name}: {str(e)}"
            logger.error(error_msg)
            import traceback
            traceback.print_exc()
            
            return {
                'gene_name': gene_name,
                'pdb_file': pdb_file,
                'status': 'error',
                'error': error_msg
            }
    
    def _generate_boltz2_configs_by_mode(self, row, pockets, atomarray, 
                                        target_dir, mode, protein_sequence, 
                                        distance_cutoff):
        """
        Generate Boltz2 configs based on execution mode.
        
        Mode 1: OS ligand only with pocket constraint (orthosteric pocket)
        Mode 2: Both ligands without pocket constraints
        """
        boltz_dir = target_dir / 'boltz2_configs'
        boltz_dir.mkdir(exist_ok=True)
        
        configs = []
        
        # Get SMILES from CSV
        os_smiles = row.get('OS SMILES', '')
        as_smiles = row.get('AS SMILES', '')
        
        if mode == 1:
            # Mode 1: OS ligand only with orthosteric pocket constraint
            logger.info("MODE 1: Orthosteric ligand with pocket constraint")
            
            if pd.notna(os_smiles) and os_smiles and len(pockets) > 0:
                # Use pocket 0 (top-ranked, assumed orthosteric)
                config = self._create_boltz2_config_with_pocket(
                    protein_sequence, 
                    os_smiles, 
                    pockets[0], 
                    atomarray,
                    ligand_chain_id='O',
                    distance_cutoff=distance_cutoff
                )
                
                config_file = boltz_dir / 'mode1_os_with_pocket.json'
                with open(config_file, 'w') as f:
                    json.dump(config, f, indent=2)
                
                configs.append({
                    'type': 'mode1_os_pocket',
                    'pocket_id': 0,
                    'file': config_file.name,
                    'description': 'Orthosteric ligand with pocket 0 constraint'
                })
                logger.info(f"Created: {config_file.name}")
            else:
                logger.warning("No OS SMILES or no pockets found for Mode 1")
        
        elif mode == 2:
            # Mode 2: Both ligands without pocket constraints
            logger.info("MODE 2: Both ligands without pocket constraints")
            
            if pd.notna(os_smiles) and os_smiles and pd.notna(as_smiles) and as_smiles:
                config = self._create_boltz2_config_no_pockets(
                    protein_sequence,
                    os_smiles,
                    as_smiles
                )
                
                config_file = boltz_dir / 'mode2_both_no_pockets.json'
                with open(config_file, 'w') as f:
                    json.dump(config, f, indent=2)
                
                configs.append({
                    'type': 'mode2_both_free',
                    'pocket_id': None,
                    'file': config_file.name,
                    'description': 'Both ligands without pocket constraints (free docking)'
                })
                logger.info(f"Created: {config_file.name}")
            else:
                logger.warning("Missing OS or AS SMILES for Mode 2")
        
        else:
            raise ValueError(f"Invalid mode: {mode}. Must be 1 or 2.")
        
        return configs
    
    def _create_boltz2_config_with_pocket(self, protein_seq, ligand_smiles, 
                                         pocket, atomarray, ligand_chain_id,
                                         distance_cutoff):
        """Create Boltz2 config with pocket constraint (Mode 1)."""
        # Get pocket contacts
        contacts = self._get_boltz2_contacts(pocket, atomarray, distance_cutoff)
        
        config = {
            "_metadata": {
                "mode": 1,
                "description": "Orthosteric ligand with pocket constraint",
                "pocket_score": float(pocket.score),
                "pocket_volume": float(pocket.volume),
                "n_contacts": len(contacts)
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
    
    def _create_boltz2_config_no_pockets(self, protein_seq, os_smiles, as_smiles):
        """Create Boltz2 config without pocket constraints (Mode 2)."""
        config = {
            "_metadata": {
                "mode": 2,
                "description": "Both ligands without pocket constraints (free docking)"
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
            # NO constraints - free docking
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
    
    async def _call_boltz2_for_configs(self, target_dir, boltz_configs):
        """Call Boltz-2 API for all configs in a target directory."""
        boltz_dir = target_dir / 'boltz2_configs'
        results_dir = target_dir / 'boltz2_results'
        results_dir.mkdir(exist_ok=True)
        
        boltz2_results = []
        
        for config_info in boltz_configs:
            config_file = boltz_dir / config_info['file']
            config_type = config_info['type']
            
            logger.info(f"Calling Boltz-2 for {config_file.name}...")
            
            try:
                # Load config
                with open(config_file) as f:
                    config_data = json.load(f)
                
                # Remove metadata if present
                if '_metadata' in config_data:
                    metadata = config_data.pop('_metadata')
                else:
                    metadata = {}
                
                # Call API
                result = await make_boltz2_api_call(
                    data=config_data,
                    api_key=self.api_key
                )
                
                # Save full result
                result_file = results_dir / f"{config_type}_result.json"
                with open(result_file, 'w') as f:
                    json.dump(result, f, indent=2)
                
                logger.info(f"Saved result to {result_file.name}")
                
                # Extract and save structures
                if 'structures' in result:
                    for i, structure_data in enumerate(result['structures']):
                        structure_content = structure_data['structure']
                        structure_format = structure_data['format']
                        
                        # Save structure file
                        structure_file = results_dir / f"{config_type}_structure_{i}.{structure_format}"
                        with open(structure_file, 'w') as f:
                            f.write(structure_content)
                        
                        logger.info(f"Saved structure to {structure_file.name}")
                
                # Store result info
                boltz2_results.append({
                    'config_type': config_type,
                    'result_file': str(result_file.name),
                    'confidence_scores': result.get('confidence_scores', []),
                    'n_structures': len(result.get('structures', []))
                })
                
                logger.info(f"✓ Boltz-2 call successful for {config_type}")
                logger.info(f"  Confidence scores: {result.get('confidence_scores', [])}")
                
            except Exception as e:
                error_msg = f"Error calling Boltz-2 for {config_file.name}: {str(e)}"
                logger.error(error_msg)
                boltz2_results.append({
                    'config_type': config_type,
                    'error': error_msg
                })
        
        return boltz2_results
    
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
        
        logger.info(f"Saved pocket info: {all_pockets_file.name}")
        
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
    
    def _create_target_summary(self, row, pockets, target_dir, boltz_configs, 
                              boltz2_results, mode):
        """Create summary file for target."""
        summary_file = target_dir / 'TARGET_SUMMARY.md'
        
        lines = [
            f"# Target Summary: {row['geneName']}",
            "",
            f"**Execution Mode**: {mode}",
            "",
            "## Target Information",
            f"- **Gene Name**: {row['geneName']}",
            f"- **UniProt ID**: {row.get('uniprotID', 'N/A')}",
            f"- **PDB File**: {row['PDB Filename']}",
            f"- **FASTA File**: {row['Fasta File']}",
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
            "## Execution Mode Details",
            ""
        ]
        
        if mode == 1:
            lines.extend([
                "**Mode 1**: Orthosteric ligand with pocket constraint",
                "- Ligand: Orthosteric only",
                "- Pocket constraint: Yes (pocket 0)",
                "- Purpose: Guided docking to orthosteric site",
                ""
            ])
        else:
            lines.extend([
                "**Mode 2**: Both ligands without pocket constraints",
                "- Ligands: Orthosteric + Allosteric",
                "- Pocket constraint: None",
                "- Purpose: Free docking of both ligands",
                ""
            ])
        
        lines.extend([
            "## Detected Pockets",
            ""
        ])
        
        for i, pocket in enumerate(pockets[:5]):
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
            lines.append(f"  - {config['description']}")
        
        lines.append("")
        
        # Add Boltz2 results if available
        if boltz2_results:
            lines.extend([
                "## Boltz-2 API Results",
                ""
            ])
            
            for result in boltz2_results:
                if 'error' not in result:
                    lines.extend([
                        f"### {result['config_type']}",
                        f"- **Result file**: `boltz2_results/{result['result_file']}`",
                        f"- **Structures**: {result['n_structures']}",
                        f"- **Confidence scores**: {result['confidence_scores']}",
                        ""
                    ])
                else:
                    lines.extend([
                        f"### {result['config_type']}",
                        f"- **Error**: {result['error']}",
                        ""
                    ])
        
        with open(summary_file, 'w') as f:
            f.write('\n'.join(lines))
        
        logger.info(f"Summary: {summary_file.name}")
    
    async def process_all(self, mode=1, distance_cutoff=5.0, max_targets=None, 
                         call_boltz2=True):
        """
        Process all targets in the CSV.
        
        Args:
            mode: Execution mode (1 or 2)
            distance_cutoff: Distance for pocket contacts
            max_targets: Optional limit on number to process
            call_boltz2: Whether to call Boltz-2 API
        """
        logger.info(f"\n{'='*80}")
        logger.info(f"BATCH PROCESSING: {len(self.targets_df)} TARGETS - MODE {mode}")
        logger.info(f"{'='*80}")
        
        if mode == 1:
            logger.info("Mode 1: Orthosteric ligand with pocket constraint")
        else:
            logger.info("Mode 2: Both ligands without pocket constraints")
        
        if not call_boltz2:
            logger.warning("Boltz-2 API calls disabled - only generating configs")
        
        n_targets = min(len(self.targets_df), max_targets) if max_targets else len(self.targets_df)
        
        for idx, row in self.targets_df.iterrows():
            if max_targets and idx >= max_targets:
                break
            
            result = await self.process_target(row, mode, distance_cutoff, call_boltz2)
            self.results.append(result)
        
        # Create batch summary
        self._create_batch_summary(mode)
        
        logger.info(f"\n{'='*80}")
        logger.info("BATCH PROCESSING COMPLETE")
        logger.info(f"{'='*80}")
        logger.info(f"Processed: {len(self.results)} targets")
        logger.info(f"Success: {sum(1 for r in self.results if r['status'] == 'success')}")
        logger.info(f"Errors: {sum(1 for r in self.results if r['status'] == 'error')}")
        logger.info(f"\nResults directory: {self.results_base}")
        logger.info(f"Batch summary: {self.results_base / f'BATCH_SUMMARY_MODE{mode}.csv'}")
    
    def _create_batch_summary(self, mode):
        """Create summary CSV of all results."""
        summary_df = pd.DataFrame(self.results)
        summary_file = self.results_base / f'BATCH_SUMMARY_MODE{mode}.csv'
        summary_df.to_csv(summary_file, index=False)
        
        logger.info(f"\nBatch summary saved: {summary_file}")


async def main():
    """Main entry point with CLI."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Batch process targets with Pocketeer → Boltz2 (two modes)'
    )
    parser.add_argument(
        '--csv',
        type=str,
        required=True,
        help='Path to CSV file with targets'
    )
    parser.add_argument(
        '--mode',
        type=int,
        choices=[1, 2],
        required=True,
        help='Execution mode: 1=OS with pocket, 2=Both without pockets'
    )
    parser.add_argument(
        '--results-dir',
        type=str,
        default='results',
        help='Base directory for results (default: results)'
    )
    parser.add_argument(
        '--api-key',
        type=str,
        help='NVIDIA API key (or set NVIDIA_API_KEY env var)'
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
    parser.add_argument(
        '--skip-api',
        action='store_true',
        help='Skip Boltz-2 API calls (only generate configs)'
    )
    
    args = parser.parse_args()
    
    # Run batch processing
    batch = BatchPocketeerBoltz2Pipeline(args.csv, args.results_dir, args.api_key)
    await batch.process_all(
        mode=args.mode,
        distance_cutoff=args.distance_cutoff,
        max_targets=args.max_targets,
        call_boltz2=not args.skip_api
    )


if __name__ == '__main__':
    asyncio.run(main())
