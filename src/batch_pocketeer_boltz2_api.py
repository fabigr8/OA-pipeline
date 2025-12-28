"""
Batch Pocketeer → Boltz-2 Pipeline with API Integration

Processes multiple targets from CSV and calls Boltz-2 NIM via NVIDIA Cloud Functions API.
Handles async API calls, long-polling, and saves results with structures.
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
    """
    Make a call to Boltz-2 via NVIDIA Cloud Functions API.
    
    Args:
        data: Request payload for Boltz-2
        api_key: NVIDIA API key
        nvcf_poll_seconds: Polling interval
        manual_timeout_seconds: Overall timeout
        
    Returns:
        dict: Boltz-2 response with structures and scores
    """
    async with httpx.AsyncClient() as client:
        headers = {
            "Authorization": f"Bearer {api_key}",
            "NVCF-POLL-SECONDS": f"{nvcf_poll_seconds}",
            "Content-Type": "application/json"
        }
        
        logger.debug(f"Making Boltz-2 API call to {PUBLIC_URL}")
        logger.debug(f"Request data: {json.dumps(data, indent=2)}")
        
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
                await asyncio.sleep(5)  # Poll every 5 seconds
                
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


class BatchPocketeerBoltz2WithAPI:
    """Batch process with integrated Boltz-2 API calling."""
    
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
    
    async def process_target(self, row, distance_cutoff=5.0, call_boltz2=True):
        """
        Process a single target from the CSV.
        
        Args:
            row: DataFrame row with target data
            distance_cutoff: Distance for pocket contacts
            call_boltz2: Whether to actually call Boltz-2 API
            
        Returns:
            dict: Processing results
        """
        gene_name = row['geneName']
        pdb_file = row['PDB Filename']
        
        logger.info(f"\n{'='*80}")
        logger.info(f"Processing: {gene_name} ({Path(pdb_file).stem})")
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
            
            # Step 3: Generate Boltz2 configs
            logger.info("Step 3: Generating Boltz2 configurations...")
            boltz_configs = self._generate_boltz2_configs(
                row, pockets, atomarray, target_dir, distance_cutoff
            )
            
            # Step 4: Call Boltz-2 API if requested
            boltz2_results = []
            if call_boltz2 and self.api_key:
                logger.info("Step 4: Calling Boltz-2 API...")
                boltz2_results = await self._call_boltz2_for_configs(
                    target_dir, boltz_configs
                )
            else:
                logger.info("Step 4: Skipping Boltz-2 API calls")
            
            # Step 5: Create target summary
            logger.info("Step 5: Creating target summary...")
            self._create_target_summary(row, pockets, target_dir, boltz_configs, 
                                       boltz2_results)
            
            result = {
                'gene_name': gene_name,
                'pdb_file': pdb_file,
                'target_dir': str(target_dir),
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
    
    async def _call_boltz2_for_configs(self, target_dir, boltz_configs):
        """Call Boltz-2 API for all configs in a target directory."""
        boltz_dir = target_dir / 'boltz2_configs'
        results_dir = target_dir / 'boltz2_results'
        results_dir.mkdir(exist_ok=True)
        
        boltz2_results = []
        
        for config_info in boltz_configs:
            config_file = boltz_dir / config_info['file']
            config_type = config_info['type']
            
            # Skip comparison configs to avoid too many API calls
            if config_type == 'comparison':
                logger.info(f"Skipping comparison config: {config_file.name}")
                continue
            
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
        if pd.notna(os_smiles) and os_smiles and len(pockets) > 0:
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
                'file': config_file.name
            })
            logger.info(f"Created: {config_file.name}")
        
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
                'file': config_file.name
            })
            logger.info(f"Created: {config_file.name}")
        
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
                'file': config_file.name
            })
            logger.info(f"Created: {config_file.name}")
        
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
    
    def _create_target_summary(self, row, pockets, target_dir, boltz_configs, 
                              boltz2_results):
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
                        f"### {result['config_type'].title()}",
                        f"- **Result file**: `boltz2_results/{result['result_file']}`",
                        f"- **Structures**: {result['n_structures']}",
                        f"- **Confidence scores**: {result['confidence_scores']}",
                        ""
                    ])
                else:
                    lines.extend([
                        f"### {result['config_type'].title()}",
                        f"- **Error**: {result['error']}",
                        ""
                    ])
        
        with open(summary_file, 'w') as f:
            f.write('\n'.join(lines))
        
        logger.info(f"Summary: {summary_file.name}")
    
    async def process_all(self, distance_cutoff=5.0, max_targets=None, call_boltz2=True):
        """
        Process all targets in the CSV.
        
        Args:
            distance_cutoff: Distance for pocket contacts
            max_targets: Optional limit on number to process
            call_boltz2: Whether to call Boltz-2 API
        """
        logger.info(f"\n{'='*80}")
        logger.info(f"BATCH PROCESSING: {len(self.targets_df)} TARGETS")
        logger.info(f"{'='*80}")
        
        if not call_boltz2:
            logger.warning("Boltz-2 API calls disabled - only generating configs")
        
        n_targets = min(len(self.targets_df), max_targets) if max_targets else len(self.targets_df)
        
        for idx, row in self.targets_df.iterrows():
            if max_targets and idx >= max_targets:
                break
            
            result = await self.process_target(row, distance_cutoff, call_boltz2)
            self.results.append(result)
        
        # Create batch summary
        self._create_batch_summary()
        
        logger.info(f"\n{'='*80}")
        logger.info("BATCH PROCESSING COMPLETE")
        logger.info(f"{'='*80}")
        logger.info(f"Processed: {len(self.results)} targets")
        logger.info(f"Success: {sum(1 for r in self.results if r['status'] == 'success')}")
        logger.info(f"Errors: {sum(1 for r in self.results if r['status'] == 'error')}")
        logger.info(f"\nResults directory: {self.results_base}")
        logger.info(f"Batch summary: {self.results_base / 'BATCH_SUMMARY.csv'}")
    
    def _create_batch_summary(self):
        """Create summary CSV of all results."""
        summary_df = pd.DataFrame(self.results)
        summary_file = self.results_base / 'BATCH_SUMMARY.csv'
        summary_df.to_csv(summary_file, index=False)
        
        logger.info(f"\nBatch summary saved: {summary_file}")


async def main():
    """Main entry point with CLI."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Batch process targets from CSV for Pocketeer → Boltz2 with API'
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
    batch = BatchPocketeerBoltz2WithAPI(args.csv, args.results_dir, args.api_key)
    await batch.process_all(
        distance_cutoff=args.distance_cutoff,
        max_targets=args.max_targets,
        call_boltz2=not args.skip_api
    )


if __name__ == '__main__':
    asyncio.run(main())
