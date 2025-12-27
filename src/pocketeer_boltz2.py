"""
Pocketeer → Boltz-2 NIM Integration

This script converts Pocketeer pocket detection output into Boltz-2 constraint format
for guided protein-ligand cofolding.

Boltz-2 supports "pocket constraints" which specify binding site interactions.
This is perfect for using Pocketeer-detected pockets to guide ligand placement.
"""

import json
import numpy as np
import pocketeer as pt
from typing import List, Dict, Optional, Tuple


class Boltz2PocketConverter:
    """Convert Pocketeer pockets to Boltz-2 constraint format."""
    
    def __init__(self, protein_path):
        """
        Initialize with protein structure.
        
        Args:
            protein_path: Path to PDB or CIF file
        """
        self.protein_path = protein_path
        self.atomarray = pt.load_structure(protein_path)
        self.pockets = None
        
    def detect_pockets(self, **kwargs):
        """
        Run Pocketeer pocket detection.
        
        Args:
            **kwargs: Parameters for pt.find_pockets()
        """
        print(f"Detecting pockets in {self.protein_path}...")
        self.pockets = pt.find_pockets(self.atomarray, **kwargs)
        print(f"Found {len(self.pockets)} pockets")
        return self.pockets
    
    def get_pocket_contacts(self, pocket_id=0, distance_cutoff=5.0, 
                           max_contacts=None):
        """
        Get pocket residues as Boltz-2 contact format.
        
        Args:
            pocket_id: Index of pocket to extract
            distance_cutoff: Distance from pocket spheres to define contacts (Å)
            max_contacts: Optional limit on number of contacts (None = all)
            
        Returns:
            list: List of contact dictionaries for Boltz-2
        """
        if self.pockets is None:
            raise ValueError("Run detect_pockets() first")
        
        pocket = self.pockets[pocket_id]
        sphere_coords = np.array([s.coord for s in pocket.spheres])
        atom_coords = self.atomarray.coord
        
        # Find all residues near pocket
        contact_residues = set()
        
        for sphere_center in sphere_coords:
            distances = np.linalg.norm(atom_coords - sphere_center, axis=1)
            nearby_indices = np.where(distances <= distance_cutoff)[0]
            
            for idx in nearby_indices:
                chain_id = self.atomarray.chain_id[idx]
                res_id = int(self.atomarray.res_id[idx])
                contact_residues.add((chain_id, res_id))
        
        # Sort by chain and residue number
        sorted_contacts = sorted(list(contact_residues), key=lambda x: (x[0], x[1]))
        
        # Limit if requested
        if max_contacts is not None:
            sorted_contacts = sorted_contacts[:max_contacts]
        
        # Convert to Boltz-2 format
        boltz2_contacts = []
        for chain_id, res_id in sorted_contacts:
            contact = {
                "id": chain_id,
                "residue_index": res_id
            }
            boltz2_contacts.append(contact)
        
        return boltz2_contacts
    
    def create_pocket_constraint(self, pocket_id=0, ligand_chain_id="L",
                                distance_cutoff=5.0, max_contacts=None):
        """
        Create a single Boltz-2 pocket constraint.
        
        Args:
            pocket_id: Index of pocket
            ligand_chain_id: Chain ID for the ligand in Boltz-2 (default: "L")
            distance_cutoff: Distance to define pocket contacts (Å)
            max_contacts: Optional limit on contacts
            
        Returns:
            dict: Boltz-2 pocket constraint
        """
        contacts = self.get_pocket_contacts(pocket_id, distance_cutoff, max_contacts)
        
        constraint = {
            "constraint_type": "pocket",
            "binder": ligand_chain_id,
            "contacts": contacts
        }
        
        return constraint
    
    def create_boltz2_request(self, 
                             protein_sequence: str,
                             ligand_smiles: str,
                             pocket_id: int = 0,
                             protein_chain_id: str = "A",
                             ligand_chain_id: str = "L",
                             distance_cutoff: float = 5.0,
                             recycling_steps: int = 3,
                             sampling_steps: int = 50,
                             diffusion_samples: int = 1):
        """
        Create a complete Boltz-2 API request with pocket constraints.
        
        Args:
            protein_sequence: Protein amino acid sequence
            ligand_smiles: SMILES string for ligand
            pocket_id: Which pocket to use as constraint
            protein_chain_id: Chain ID for protein
            ligand_chain_id: Chain ID for ligand
            distance_cutoff: Distance to define pocket
            recycling_steps: Boltz-2 recycling steps (1-6)
            sampling_steps: Boltz-2 sampling steps (10-1000)
            diffusion_samples: Number of diffusion samples (1-5)
            
        Returns:
            dict: Complete Boltz-2 API request
        """
        # Create pocket constraint
        pocket_constraint = self.create_pocket_constraint(
            pocket_id=pocket_id,
            ligand_chain_id=ligand_chain_id,
            distance_cutoff=distance_cutoff
        )
        
        # Build Boltz-2 request
        request = {
            "polymers": [
                {
                    "id": protein_chain_id,
                    "molecule_type": "protein",
                    "sequence": protein_sequence
                }
            ],
            "ligands": [
                {
                    "id": ligand_chain_id,
                    "smiles": ligand_smiles
                }
            ],
            "constraints": [pocket_constraint],
            "recycling_steps": recycling_steps,
            "sampling_steps": sampling_steps,
            "diffusion_samples": diffusion_samples,
            "output_format": "mmcif"
        }
        
        return request
    
    def create_multi_pocket_requests(self,
                                    protein_sequence: str,
                                    ligand_smiles: str,
                                    pocket_ids: List[int],
                                    protein_chain_id: str = "A",
                                    ligand_chain_id: str = "L",
                                    distance_cutoff: float = 5.0,
                                    **boltz_params):
        """
        Create multiple Boltz-2 requests for different pockets.
        Useful for testing orthosteric vs allosteric binding.
        
        Args:
            protein_sequence: Protein sequence
            ligand_smiles: Ligand SMILES
            pocket_ids: List of pocket indices to test
            protein_chain_id: Protein chain ID
            ligand_chain_id: Ligand chain ID  
            distance_cutoff: Distance for pocket definition
            **boltz_params: Additional Boltz-2 parameters
            
        Returns:
            list: List of Boltz-2 requests, one per pocket
        """
        requests = []
        
        for pocket_id in pocket_ids:
            request = self.create_boltz2_request(
                protein_sequence=protein_sequence,
                ligand_smiles=ligand_smiles,
                pocket_id=pocket_id,
                protein_chain_id=protein_chain_id,
                ligand_chain_id=ligand_chain_id,
                distance_cutoff=distance_cutoff,
                **boltz_params
            )
            
            # Add metadata
            request['_metadata'] = {
                'pocket_id': pocket_id,
                'pocket_score': float(self.pockets[pocket_id].score),
                'pocket_volume': float(self.pockets[pocket_id].volume)
            }
            
            requests.append(request)
        
        return requests
    
    def export_boltz2_configs(self, 
                             protein_sequence: str,
                             ligand_smiles: str,
                             output_dir: str = "boltz2_configs",
                             max_pockets: int = 3,
                             distance_cutoff: float = 5.0):
        """
        Export Boltz-2 configuration files for top pockets.
        
        Args:
            protein_sequence: Protein sequence
            ligand_smiles: Ligand SMILES
            output_dir: Output directory
            max_pockets: Maximum number of pockets to export
            distance_cutoff: Distance for pocket contacts
        """
        from pathlib import Path
        
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        n_pockets = min(len(self.pockets), max_pockets)
        
        for i in range(n_pockets):
            request = self.create_boltz2_request(
                protein_sequence=protein_sequence,
                ligand_smiles=ligand_smiles,
                pocket_id=i,
                distance_cutoff=distance_cutoff
            )
            
            # Save to JSON file
            output_file = output_path / f'pocket_{i}_boltz2.json'
            with open(output_file, 'w') as f:
                json.dump(request, f, indent=2)
            
            print(f"Created: {output_file}")
        
        # Create summary
        self._create_boltz2_summary(output_path, n_pockets, distance_cutoff)
    
    def _create_boltz2_summary(self, output_path, n_pockets, distance_cutoff):
        """Create summary file for Boltz-2 configs."""
        summary_file = output_path / 'README.txt'
        
        with open(summary_file, 'w') as f:
            f.write("BOLTZ-2 POCKET CONSTRAINT CONFIGURATIONS\n")
            f.write("=" * 70 + "\n\n")
            f.write(f"Generated {n_pockets} configurations from Pocketeer pockets\n")
            f.write(f"Distance cutoff: {distance_cutoff} Å\n\n")
            
            for i in range(n_pockets):
                pocket = self.pockets[i]
                sphere_coords = np.array([s.coord for s in pocket.spheres])
                center = np.mean(sphere_coords, axis=0)
                
                f.write(f"Pocket {i}:\n")
                f.write(f"  File: pocket_{i}_boltz2.json\n")
                f.write(f"  Score: {pocket.score:.3f}\n")
                f.write(f"  Volume: {pocket.volume:.1f} Å³\n")
                f.write(f"  Center: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})\n")
                f.write(f"  Spheres: {pocket.n_spheres}\n\n")
            
            f.write("\nUSAGE:\n")
            f.write("1. Choose a pocket configuration (pocket_0_boltz2.json recommended)\n")
            f.write("2. Send as POST request to Boltz-2 NIM endpoint\n")
            f.write("3. Compare results from different pockets if testing multiple sites\n")
        
        print(f"Summary: {summary_file}")
    
    def print_constraint_summary(self, pocket_id=0, distance_cutoff=5.0):
        """Print summary of pocket constraint."""
        contacts = self.get_pocket_contacts(pocket_id, distance_cutoff)
        pocket = self.pockets[pocket_id]
        sphere_coords = np.array([s.coord for s in pocket.spheres])
        center = np.mean(sphere_coords, axis=0)
        
        print(f"\nPOCKET {pocket_id} CONSTRAINT SUMMARY")
        print("=" * 70)
        print(f"Score:           {pocket.score:.3f}")
        print(f"Volume:          {pocket.volume:.1f} Å³")
        print(f"Center:          ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
        print(f"Contact residues: {len(contacts)}")
        print(f"\nFirst 10 contacts:")
        for contact in contacts[:10]:
            print(f"  Chain {contact['id']}, Residue {contact['residue_index']}")
        
        if len(contacts) > 10:
            print(f"  ... and {len(contacts) - 10} more")


def example_usage():
    """Example of using Boltz-2 converter."""
    
    print("=" * 80)
    print("POCKETEER → BOLTZ-2 NIM INTEGRATION")
    print("=" * 80)
    
    # 1. Initialize and detect pockets
    converter = Boltz2PocketConverter("protein.pdb")
    converter.detect_pockets(
        r_min=3.0,
        r_max=6.0,
        min_spheres=30
    )
    
    # 2. Example protein sequence and ligand
    # Replace with your actual sequence and SMILES
    protein_seq = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTT"
    ligand_smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"  # Ibuprofen example
    
    # 3. Create single pocket constraint (orthosteric)
    print("\n" + "=" * 80)
    print("OPTION 1: Single Pocket (Orthosteric Site)")
    print("=" * 80)
    
    request_orthosteric = converter.create_boltz2_request(
        protein_sequence=protein_seq,
        ligand_smiles=ligand_smiles,
        pocket_id=0,  # Top-ranked pocket
        protein_chain_id="A",
        ligand_chain_id="L",
        distance_cutoff=5.0,
        recycling_steps=3,
        sampling_steps=50,
        diffusion_samples=1
    )
    
    print(f"\nRequest for pocket 0 (orthosteric):")
    print(f"  Constraints: {len(request_orthosteric['constraints'])}")
    print(f"  Contacts in pocket: {len(request_orthosteric['constraints'][0]['contacts'])}")
    
    # 4. Create multiple requests for allosteric vs orthosteric
    print("\n" + "=" * 80)
    print("OPTION 2: Multiple Pockets (Allosteric vs Orthosteric)")
    print("=" * 80)
    
    if len(converter.pockets) >= 2:
        requests = converter.create_multi_pocket_requests(
            protein_sequence=protein_seq,
            ligand_smiles=ligand_smiles,
            pocket_ids=[0, 1],  # Test top 2 pockets
            distance_cutoff=5.0,
            recycling_steps=3,
            sampling_steps=50,
            diffusion_samples=1
        )
        
        print(f"\nCreated {len(requests)} requests for different pockets:")
        for i, req in enumerate(requests):
            meta = req['_metadata']
            print(f"  Pocket {meta['pocket_id']}: Score={meta['pocket_score']:.2f}, "
                  f"Volume={meta['pocket_volume']:.1f} Å³")
    
    # 5. Export all configs
    print("\n" + "=" * 80)
    print("OPTION 3: Export All Configurations")
    print("=" * 80)
    
    converter.export_boltz2_configs(
        protein_sequence=protein_seq,
        ligand_smiles=ligand_smiles,
        output_dir="boltz2_configs",
        max_pockets=3,
        distance_cutoff=5.0
    )
    
    # 6. Print detailed summary
    print("\n" + "=" * 80)
    print("POCKET CONSTRAINT DETAILS")
    print("=" * 80)
    
    converter.print_constraint_summary(pocket_id=0, distance_cutoff=5.0)
    
    print("\n" + "=" * 80)
    print("NEXT STEPS")
    print("=" * 80)
    print("""
1. Review boltz2_configs/ directory
2. Choose best pocket configuration
3. Send JSON to Boltz-2 NIM endpoint:
   
   POST /biology/mit/boltz2/predict
   
   with the JSON from pocket_0_boltz2.json (or other pocket)
   
4. Compare results from different pockets if testing multiple sites
5. Boltz-2 will use the pocket contacts to guide ligand placement
    """)


def create_sample_request():
    """Create a minimal working example."""
    
    sample_request = {
        "polymers": [
            {
                "id": "A",
                "molecule_type": "protein",
                "sequence": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTT"
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
                    # Add more contacts from Pocketeer output
                ]
            }
        ],
        "recycling_steps": 3,
        "sampling_steps": 50,
        "diffusion_samples": 1,
        "output_format": "mmcif"
    }
    
    print("\nSAMPLE BOLTZ-2 REQUEST:")
    print(json.dumps(sample_request, indent=2))
    
    return sample_request


if __name__ == "__main__":
    # Run example
    example_usage()
    
    # Show minimal sample
    print("\n" + "=" * 80)
    create_sample_request()
