"""
Complete Workflow: Pocketeer → Boltz-2 NIM for Allosteric/Orthosteric Docking

This script demonstrates the full pipeline:
1. Detect pockets with Pocketeer
2. Classify as orthosteric/allosteric
3. Generate Boltz-2 requests with pocket constraints
4. Handle multiple ligands for different sites

Example use case: 
- Orthosteric ligand (competitive inhibitor)
- Allosteric ligand (modulator)
"""

import json
import requests
from pathlib import Path
from pocketeer_boltz2 import Boltz2PocketConverter


class AllostericOrthostericPipeline:
    """Pipeline for handling multiple binding sites with Boltz-2."""
    
    def __init__(self, protein_pdb, protein_sequence):
        """
        Initialize pipeline.
        
        Args:
            protein_pdb: Path to protein PDB file
            protein_sequence: Full protein amino acid sequence
        """
        self.converter = Boltz2PocketConverter(protein_pdb)
        self.protein_sequence = protein_sequence
        self.pockets = None
        
    def detect_and_classify_pockets(self, known_orthosteric_site=None):
        """
        Detect pockets and classify as orthosteric/allosteric.
        
        Args:
            known_orthosteric_site: Optional (x,y,z) coordinates of known active site
            
        Returns:
            dict: Classification results
        """
        import numpy as np
        
        # Detect pockets
        self.pockets = self.converter.detect_pockets(
            r_min=3.0,
            r_max=6.0,
            min_spheres=30,
            merge_distance=2.5
        )
        
        # Classify pockets
        classification = {
            'orthosteric': [],
            'allosteric': []
        }
        
        if known_orthosteric_site is None:
            # Assume top-ranked is orthosteric
            classification['orthosteric'].append(0)
            classification['allosteric'].extend(range(1, len(self.pockets)))
        else:
            # Classify by distance to known site
            known_site = np.array(known_orthosteric_site)
            
            for i, pocket in enumerate(self.pockets):
                sphere_coords = np.array([s.coord for s in pocket.spheres])
                pocket_center = np.mean(sphere_coords, axis=0)
                distance = np.linalg.norm(pocket_center - known_site)
                
                if distance < 15.0:  # Within 15Å = orthosteric
                    classification['orthosteric'].append(i)
                else:
                    classification['allosteric'].append(i)
        
        return classification
    
    def create_dual_ligand_request(self,
                                   orthosteric_smiles,
                                   allosteric_smiles,
                                   orthosteric_pocket_id=0,
                                   allosteric_pocket_id=1,
                                   distance_cutoff=5.0):
        """
        Create Boltz-2 request with both orthosteric and allosteric ligands.
        
        Args:
            orthosteric_smiles: SMILES for orthosteric ligand
            allosteric_smiles: SMILES for allosteric ligand
            orthosteric_pocket_id: Pocket ID for orthosteric site
            allosteric_pocket_id: Pocket ID for allosteric site
            distance_cutoff: Distance for pocket contacts
            
        Returns:
            dict: Boltz-2 request with dual constraints
        """
        # Create constraints for both pockets
        ortho_constraint = self.converter.create_pocket_constraint(
            pocket_id=orthosteric_pocket_id,
            ligand_chain_id="O",  # Orthosteric ligand
            distance_cutoff=distance_cutoff
        )
        
        allo_constraint = self.converter.create_pocket_constraint(
            pocket_id=allosteric_pocket_id,
            ligand_chain_id="A",  # Allosteric ligand
            distance_cutoff=distance_cutoff
        )
        
        # Build request
        request = {
            "polymers": [
                {
                    "id": "P",
                    "molecule_type": "protein",
                    "sequence": self.protein_sequence
                }
            ],
            "ligands": [
                {
                    "id": "O",
                    "smiles": orthosteric_smiles
                },
                {
                    "id": "A",
                    "smiles": allosteric_smiles
                }
            ],
            "constraints": [
                ortho_constraint,
                allo_constraint
            ],
            "recycling_steps": 4,  # More recycling for dual ligands
            "sampling_steps": 100,  # More sampling for complexity
            "diffusion_samples": 2,
            "output_format": "mmcif"
        }
        
        return request
    
    def create_comparison_requests(self, ligand_smiles, pocket_ids=None,
                                  distance_cutoff=5.0):
        """
        Create multiple requests to test ligand in different pockets.
        
        Args:
            ligand_smiles: SMILES string for ligand
            pocket_ids: List of pocket IDs to test (None = all)
            distance_cutoff: Distance for pocket contacts
            
        Returns:
            list: Multiple Boltz-2 requests
        """
        if pocket_ids is None:
            pocket_ids = list(range(min(3, len(self.pockets))))
        
        requests = []
        
        for pocket_id in pocket_ids:
            request = self.converter.create_boltz2_request(
                protein_sequence=self.protein_sequence,
                ligand_smiles=ligand_smiles,
                pocket_id=pocket_id,
                protein_chain_id="P",
                ligand_chain_id="L",
                distance_cutoff=distance_cutoff,
                recycling_steps=3,
                sampling_steps=50,
                diffusion_samples=1
            )
            
            requests.append({
                'pocket_id': pocket_id,
                'request': request,
                'pocket_info': {
                    'score': float(self.pockets[pocket_id].score),
                    'volume': float(self.pockets[pocket_id].volume)
                }
            })
        
        return requests
    
    def export_workflow(self, output_dir="boltz2_workflow"):
        """Export complete workflow with configurations."""
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        # Create workflow guide
        guide_path = output_path / "workflow_guide.md"
        with open(guide_path, 'w') as f:
            f.write(self._generate_workflow_guide())
        
        print(f"Workflow guide: {guide_path}")
        
    def _generate_workflow_guide(self):
        """Generate markdown workflow guide."""
        import numpy as np
        
        guide = [
            "# Pocketeer → Boltz-2 Workflow for Allosteric/Orthosteric Sites",
            "",
            "## Detected Pockets",
            ""
        ]
        
        for i, pocket in enumerate(self.pockets):
            sphere_coords = np.array([s.coord for s in pocket.spheres])
            center = np.mean(sphere_coords, axis=0)
            
            guide.extend([
                f"### Pocket {i}",
                f"- **Score**: {pocket.score:.3f}",
                f"- **Volume**: {pocket.volume:.1f} Å³",
                f"- **Center**: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})",
                f"- **Spheres**: {pocket.n_spheres}",
                ""
            ])
        
        guide.extend([
            "## Workflow Options",
            "",
            "### Option 1: Test Single Ligand in Multiple Pockets",
            "```python",
            "# Test if ligand prefers orthosteric or allosteric site",
            "requests = pipeline.create_comparison_requests(",
            "    ligand_smiles='YOUR_SMILES',",
            "    pocket_ids=[0, 1]  # Test top 2 pockets",
            ")",
            "",
            "# Send each request to Boltz-2",
            "# Compare confidence scores to determine preferred site",
            "```",
            "",
            "### Option 2: Dual Ligand Binding",
            "```python",
            "# Model both orthosteric and allosteric ligands simultaneously",
            "request = pipeline.create_dual_ligand_request(",
            "    orthosteric_smiles='COMPETITIVE_INHIBITOR_SMILES',",
            "    allosteric_smiles='MODULATOR_SMILES',",
            "    orthosteric_pocket_id=0,",
            "    allosteric_pocket_id=1",
            ")",
            "",
            "# Send to Boltz-2 to model both ligands at once",
            "```",
            "",
            "### Option 3: Single Site Focus",
            "```python",
            "# Focus only on one pocket (e.g., validated orthosteric site)",
            "request = converter.create_boltz2_request(",
            "    protein_sequence=sequence,",
            "    ligand_smiles='YOUR_SMILES',",
            "    pocket_id=0  # Just orthosteric",
            ")",
            "```",
            ""
        ])
        
        return '\n'.join(guide)


def example_orthosteric_allosteric():
    """Example: Different ligands for orthosteric vs allosteric sites."""
    
    print("=" * 80)
    print("EXAMPLE: ORTHOSTERIC + ALLOSTERIC LIGAND DOCKING")
    print("=" * 80)
    
    # Example protein (replace with your actual data)
    protein_pdb = "kinase.pdb"
    protein_sequence = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTT"
    
    # Example ligands
    atp_competitive = "CC1=C2N=CN(C2=C(NC1=O)N)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O"  # ATP analog
    allosteric_modulator = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"  # Example allosteric compound
    
    # Initialize pipeline
    pipeline = AllostericOrthostericPipeline(protein_pdb, protein_sequence)
    
    # Detect and classify pockets
    print("\nStep 1: Detecting pockets...")
    classification = pipeline.detect_and_classify_pockets()
    
    print(f"\nClassification results:")
    print(f"  Orthosteric pockets: {classification['orthosteric']}")
    print(f"  Allosteric pockets:  {classification['allosteric']}")
    
    # Create dual ligand request
    print("\nStep 2: Creating dual-ligand Boltz-2 request...")
    if classification['orthosteric'] and classification['allosteric']:
        ortho_id = classification['orthosteric'][0]
        allo_id = classification['allosteric'][0]
        
        dual_request = pipeline.create_dual_ligand_request(
            orthosteric_smiles=atp_competitive,
            allosteric_smiles=allosteric_modulator,
            orthosteric_pocket_id=ortho_id,
            allosteric_pocket_id=allo_id
        )
        
        print(f"\nDual-ligand request created:")
        print(f"  Orthosteric constraint: Pocket {ortho_id}, "
              f"{len(dual_request['constraints'][0]['contacts'])} contacts")
        print(f"  Allosteric constraint:  Pocket {allo_id}, "
              f"{len(dual_request['constraints'][1]['contacts'])} contacts")
        
        # Save request
        with open('dual_ligand_boltz2.json', 'w') as f:
            json.dump(dual_request, f, indent=2)
        print("\nSaved to: dual_ligand_boltz2.json")
    
    # Create comparison requests
    print("\nStep 3: Creating comparison requests...")
    comparison_requests = pipeline.create_comparison_requests(
        ligand_smiles=atp_competitive,
        pocket_ids=[0, 1]
    )
    
    for req in comparison_requests:
        filename = f"pocket_{req['pocket_id']}_comparison.json"
        with open(filename, 'w') as f:
            json.dump(req['request'], f, indent=2)
        print(f"  Created: {filename} (Score: {req['pocket_info']['score']:.2f})")
    
    # Export workflow
    print("\nStep 4: Exporting workflow guide...")
    pipeline.export_workflow()
    
    print("\n" + "=" * 80)
    print("PIPELINE COMPLETE")
    print("=" * 80)
    print("""
Next steps:
1. Review generated JSON files
2. Send requests to Boltz-2 NIM endpoint:
   
   curl -X POST https://YOUR_NIM_ENDPOINT/biology/mit/boltz2/predict \\
        -H "Content-Type: application/json" \\
        -d @dual_ligand_boltz2.json
   
3. For comparison mode:
   - Run both pocket_0_comparison.json and pocket_1_comparison.json
   - Compare confidence_scores to see which site ligand prefers
   
4. For dual ligand mode:
   - Use dual_ligand_boltz2.json
   - Boltz-2 will model both ligands simultaneously
    """)


def send_to_boltz2(request_data, nim_endpoint, api_key=None):
    """
    Send request to Boltz-2 NIM (example function).
    
    Args:
        request_data: Boltz-2 request dictionary
        nim_endpoint: URL of Boltz-2 NIM endpoint
        api_key: Optional API key
        
    Returns:
        dict: Boltz-2 response
    """
    headers = {"Content-Type": "application/json"}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"
    
    response = requests.post(
        f"{nim_endpoint}/biology/mit/boltz2/predict",
        json=request_data,
        headers=headers
    )
    
    response.raise_for_status()
    return response.json()


if __name__ == "__main__":
    example_orthosteric_allosteric()
