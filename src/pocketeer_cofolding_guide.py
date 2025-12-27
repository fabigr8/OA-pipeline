"""
Integration script for using Pocketeer pocket detection with cofolding workflows.

This script shows how to:
1. Identify binding pockets before cofolding
2. Create spatial constraints for ligand placement
3. Generate input for AlphaFold3 or other cofolding tools
4. Handle multiple pockets (allosteric/orthosteric)
"""

import json
import numpy as np
import pocketeer as pt
from dataclasses import dataclass
from typing import List, Tuple, Optional


@dataclass
class PocketConstraint:
    """Constraint information for a binding pocket."""
    pocket_id: int
    center: np.ndarray
    radius: float
    residues: List[Tuple[str, int, str]]
    score: float
    volume: float
    pocket_type: str  # 'orthosteric', 'allosteric', or 'unknown'


class CofoldingPocketGuide:
    """Guide cofolding using Pocketeer-detected pockets."""
    
    def __init__(self, protein_path):
        self.protein_path = protein_path
        self.atomarray = pt.load_structure(protein_path)
        self.pockets = None
        self.constraints = []
        
    def detect_pockets(self, **kwargs):
        """Run pocket detection."""
        print(f"Detecting pockets in {self.protein_path}...")
        self.pockets = pt.find_pockets(self.atomarray, **kwargs)
        print(f"Detected {len(self.pockets)} pockets")
        return self.pockets
    
    def classify_pockets(self, known_active_site=None, distance_threshold=15.0):
        """
        Classify pockets as orthosteric or allosteric.
        
        Args:
            known_active_site: (x, y, z) coordinates of known active site
            distance_threshold: Distance in Å to classify as orthosteric
            
        Returns:
            list: Classified pocket types
        """
        if known_active_site is None:
            # If no known site, assume top-ranked is orthosteric
            return ['orthosteric'] + ['allosteric'] * (len(self.pockets) - 1)
        
        classifications = []
        known_site = np.array(known_active_site)
        
        for i, pocket in enumerate(self.pockets):
            # Get pocket center
            sphere_coords = np.array([s.coord for s in pocket.spheres])
            pocket_center = np.mean(sphere_coords, axis=0)
            
            # Calculate distance to known active site
            distance = np.linalg.norm(pocket_center - known_site)
            
            if distance < distance_threshold:
                classifications.append('orthosteric')
            else:
                classifications.append('allosteric')
        
        return classifications
    
    def create_constraints(self, distance_cutoff=5.0, known_active_site=None):
        """
        Create cofolding constraints from detected pockets.
        
        Args:
            distance_cutoff: Distance to define pocket residues
            known_active_site: Optional known active site coordinates
        """
        classifications = self.classify_pockets(known_active_site)
        
        for i, pocket in enumerate(self.pockets):
            # Get pocket information
            sphere_coords = np.array([s.coord for s in pocket.spheres])
            center = np.mean(sphere_coords, axis=0)
            
            # Calculate effective radius
            distances = np.linalg.norm(sphere_coords - center, axis=0)
            radius = np.max(distances)
            
            # Get pocket residues
            residues = self._get_pocket_residues(pocket, distance_cutoff)
            
            # Create constraint
            constraint = PocketConstraint(
                pocket_id=i,
                center=center,
                radius=radius,
                residues=residues,
                score=pocket.score,
                volume=pocket.volume,
                pocket_type=classifications[i]
            )
            
            self.constraints.append(constraint)
        
        return self.constraints
    
    def _get_pocket_residues(self, pocket, distance_cutoff):
        """Helper to extract residues near pocket."""
        sphere_coords = np.array([s.coord for s in pocket.spheres])
        atom_coords = self.atomarray.coord
        
        pocket_residues = set()
        
        for sphere_center in sphere_coords:
            distances = np.linalg.norm(atom_coords - sphere_center, axis=1)
            nearby_indices = np.where(distances <= distance_cutoff)[0]
            
            for idx in nearby_indices:
                chain_id = self.atomarray.chain_id[idx]
                res_num = self.atomarray.res_id[idx]
                res_name = self.atomarray.res_name[idx]
                pocket_residues.add((chain_id, res_num, res_name))
        
        return sorted(list(pocket_residues), key=lambda x: (x[0], x[1]))
    
    def export_alphafold_constraints(self, output_file='af_constraints.json',
                                    ligand_name='LIG', include_allosteric=True):
        """
        Export constraints for AlphaFold3 or similar tools.
        
        Args:
            output_file: Output JSON file
            ligand_name: Name/identifier for ligand
            include_allosteric: Whether to include allosteric pockets
        """
        af_constraints = {
            'pockets': [],
            'recommended_strategy': {
                'description': 'Multi-pocket cofolding strategy',
                'steps': []
            }
        }
        
        for constraint in self.constraints:
            # Skip allosteric if requested
            if not include_allosteric and constraint.pocket_type == 'allosteric':
                continue
            
            pocket_data = {
                'pocket_id': constraint.pocket_id,
                'pocket_type': constraint.pocket_type,
                'center': constraint.center.tolist(),
                'radius': float(constraint.radius),
                'score': float(constraint.score),
                'volume': float(constraint.volume),
                'residues': [
                    {
                        'chain': r[0],
                        'resid': int(r[1]),
                        'resname': r[2]
                    }
                    for r in constraint.residues
                ],
                'ligand_placement_hint': {
                    'center_x': float(constraint.center[0]),
                    'center_y': float(constraint.center[1]),
                    'center_z': float(constraint.center[2]),
                    'search_radius': float(constraint.radius + 5.0)
                }
            }
            
            af_constraints['pockets'].append(pocket_data)
        
        # Add strategy recommendations
        orthosteric = [c for c in self.constraints if c.pocket_type == 'orthosteric']
        allosteric = [c for c in self.constraints if c.pocket_type == 'allosteric']
        
        if orthosteric:
            af_constraints['recommended_strategy']['steps'].append({
                'step': 1,
                'description': f'Cofold with ligand in orthosteric pocket (pocket {orthosteric[0].pocket_id})',
                'pocket_id': orthosteric[0].pocket_id,
                'ligand_constraint': {
                    'center': orthosteric[0].center.tolist(),
                    'radius': float(orthosteric[0].radius)
                }
            })
        
        if allosteric and include_allosteric:
            for i, allost in enumerate(allosteric[:2]):  # Limit to top 2 allosteric
                af_constraints['recommended_strategy']['steps'].append({
                    'step': 2 + i,
                    'description': f'Optionally cofold with ligand in allosteric pocket (pocket {allost.pocket_id})',
                    'pocket_id': allost.pocket_id,
                    'ligand_constraint': {
                        'center': allost.center.tolist(),
                        'radius': float(allost.radius)
                    }
                })
        
        with open(output_file, 'w') as f:
            json.dump(af_constraints, f, indent=2)
        
        print(f"AlphaFold constraints written to {output_file}")
    
    def create_residue_mask(self, pocket_id=0, expand_shells=2):
        """
        Create a residue mask for focused cofolding.
        
        Args:
            pocket_id: Which pocket to create mask for
            expand_shells: Number of neighbor shells to include
            
        Returns:
            list: Residue IDs to include in cofolding
        """
        constraint = self.constraints[pocket_id]
        core_residues = set(constraint.residues)
        
        # Expand to neighboring residues
        all_residues = core_residues.copy()
        
        for _ in range(expand_shells):
            new_residues = set()
            for chain, resnum, resname in all_residues:
                # Add sequential neighbors
                new_residues.add((chain, resnum - 1, None))
                new_residues.add((chain, resnum + 1, None))
            
            all_residues.update(new_residues)
        
        # Filter to only existing residues
        existing_res = set(
            (self.atomarray.chain_id[i], self.atomarray.res_id[i])
            for i in range(len(self.atomarray))
        )
        
        final_residues = [
            (chain, resnum) for chain, resnum, _ in all_residues
            if (chain, resnum) in existing_res
        ]
        
        return sorted(final_residues, key=lambda x: (x[0], x[1]))
    
    def export_cofolding_plan(self, output_file='cofolding_plan.md'):
        """
        Export a markdown file with cofolding strategy.
        
        Args:
            output_file: Output markdown file
        """
        lines = [
            "# Cofolding Plan Generated from Pocketeer",
            "",
            f"**Protein:** {self.protein_path}",
            f"**Pockets Detected:** {len(self.constraints)}",
            "",
            "---",
            ""
        ]
        
        # Orthosteric pockets
        orthosteric = [c for c in self.constraints if c.pocket_type == 'orthosteric']
        if orthosteric:
            lines.extend([
                "## Orthosteric Pocket(s)",
                "",
                "*Primary binding site(s) - highest priority for cofolding*",
                ""
            ])
            
            for constraint in orthosteric:
                lines.extend([
                    f"### Pocket {constraint.pocket_id}",
                    f"- **Score:** {constraint.score:.2f}",
                    f"- **Volume:** {constraint.volume:.1f} Å³",
                    f"- **Center:** ({constraint.center[0]:.2f}, {constraint.center[1]:.2f}, {constraint.center[2]:.2f})",
                    f"- **Radius:** {constraint.radius:.2f} Å",
                    f"- **Lining Residues:** {len(constraint.residues)}",
                    "",
                    "**Recommended approach:**",
                    "1. Use this pocket for primary ligand placement",
                    "2. Center ligand search around the coordinates above",
                    f"3. Constrain ligand within {constraint.radius + 5:.1f} Å radius",
                    ""
                ])
        
        # Allosteric pockets
        allosteric = [c for c in self.constraints if c.pocket_type == 'allosteric']
        if allosteric:
            lines.extend([
                "---",
                "",
                "## Allosteric Pocket(s)",
                "",
                "*Secondary binding sites - useful for multi-ligand systems*",
                ""
            ])
            
            for constraint in allosteric[:3]:  # Show top 3
                lines.extend([
                    f"### Pocket {constraint.pocket_id}",
                    f"- **Score:** {constraint.score:.2f}",
                    f"- **Volume:** {constraint.volume:.1f} Å³",
                    f"- **Center:** ({constraint.center[0]:.2f}, {constraint.center[1]:.2f}, {constraint.center[2]:.2f})",
                    f"- **Radius:** {constraint.radius:.2f} Å",
                    f"- **Lining Residues:** {len(constraint.residues)}",
                    ""
                ])
        
        # General strategy
        lines.extend([
            "---",
            "",
            "## Suggested Cofolding Strategy",
            "",
            "### Option 1: Single Ligand (Orthosteric)",
            "1. Focus on the top-ranked orthosteric pocket",
            "2. Constrain ligand placement to pocket center ± radius",
            "3. Include pocket-lining residues + 2 shells of neighbors",
            "",
            "### Option 2: Multi-Ligand (Orthosteric + Allosteric)",
            "1. Run separate cofolding for each pocket",
            "2. Compare binding energies/scores",
            "3. Consider simultaneous binding if pockets are distant",
            "",
            "### Option 3: Full Protein",
            "1. Use pocket information to bias search",
            "2. Allow flexibility for novel binding modes",
            "3. Validate against detected pockets",
            ""
        ])
        
        with open(output_file, 'w') as f:
            f.write('\n'.join(lines))
        
        print(f"Cofolding plan written to {output_file}")
    
    def print_summary(self):
        """Print summary of pockets and constraints."""
        print("\n" + "="*70)
        print("COFOLDING POCKET GUIDE SUMMARY")
        print("="*70)
        
        orthosteric = [c for c in self.constraints if c.pocket_type == 'orthosteric']
        allosteric = [c for c in self.constraints if c.pocket_type == 'allosteric']
        
        print(f"\nTotal pockets: {len(self.constraints)}")
        print(f"  Orthosteric: {len(orthosteric)}")
        print(f"  Allosteric:  {len(allosteric)}")
        
        print("\nPocket Details:")
        for constraint in self.constraints:
            print(f"\n  Pocket {constraint.pocket_id} [{constraint.pocket_type.upper()}]")
            print(f"    Score:     {constraint.score:.2f}")
            print(f"    Volume:    {constraint.volume:.1f} Å³")
            print(f"    Center:    ({constraint.center[0]:.1f}, {constraint.center[1]:.1f}, {constraint.center[2]:.1f})")
            print(f"    Radius:    {constraint.radius:.1f} Å")
            print(f"    Residues:  {len(constraint.residues)}")


def main():
    """Example workflow."""
    
    # Initialize
    protein_file = "protein.pdb"
    guide = CofoldingPocketGuide(protein_file)
    
    # Detect pockets
    guide.detect_pockets(
        r_min=3.0,
        r_max=6.0,
        min_spheres=30
    )
    
    # Create constraints
    # Option 1: No known active site (auto-classify)
    guide.create_constraints(distance_cutoff=5.0)
    
    # Option 2: With known active site
    # known_site = [10.5, 20.3, 15.7]  # Replace with actual coordinates
    # guide.create_constraints(distance_cutoff=5.0, known_active_site=known_site)
    
    # Print summary
    guide.print_summary()
    
    # Export files
    print("\n" + "="*70)
    print("EXPORTING COFOLDING GUIDES")
    print("="*70 + "\n")
    
    guide.export_alphafold_constraints('alphafold_constraints.json')
    guide.export_cofolding_plan('cofolding_plan.md')
    
    # Example: Get residue mask for focused cofolding
    if len(guide.constraints) > 0:
        residue_mask = guide.create_residue_mask(pocket_id=0, expand_shells=2)
        print(f"\nResidue mask for pocket 0: {len(residue_mask)} residues")
        print(f"  Example residues: {residue_mask[:10]}")
    
    print("\n✓ Cofolding guides generated successfully!")


if __name__ == "__main__":
    main()
