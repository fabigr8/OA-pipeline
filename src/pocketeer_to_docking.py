"""
Extract binding pocket information from Pocketeer for docking/cofolding pipelines.

This script demonstrates how to:
1. Run Pocketeer to detect binding pockets
2. Extract sphere coordinates and pocket information
3. Convert to docking box coordinates
4. Identify pocket residues for cofolding constraints
5. Export in various formats
"""

import json
import numpy as np
import pocketeer as pt
from pathlib import Path


class PocketExtractor:
    """Extract and process pocket information from Pocketeer results."""
    
    def __init__(self, protein_path):
        """
        Initialize with protein structure.
        
        Args:
            protein_path: Path to PDB or CIF file
        """
        self.protein_path = protein_path
        self.atomarray = pt.load_structure(protein_path)
        self.pockets = None
        
    def find_pockets(self, **kwargs):
        """
        Run Pocketeer pocket detection.
        
        Args:
            **kwargs: Parameters to pass to pt.find_pockets()
                     (r_min, r_max, min_spheres, etc.)
        """
        print(f"Finding pockets in {self.protein_path}...")
        self.pockets = pt.find_pockets(self.atomarray, **kwargs)
        print(f"Found {len(self.pockets)} pockets")
        return self.pockets
    
    def get_pocket_center(self, pocket_id=0):
        """
        Get the geometric center of a pocket.
        
        Args:
            pocket_id: Index of pocket (0 = top-ranked)
            
        Returns:
            np.array: [x, y, z] coordinates of pocket center
        """
        pocket = self.pockets[pocket_id]
        
        # Get all sphere centers
        sphere_coords = np.array([sphere.coord for sphere in pocket.spheres])
        
        # Calculate geometric center
        center = np.mean(sphere_coords, axis=0)
        
        return center
    
    def get_docking_box(self, pocket_id=0, padding=5.0):
        """
        Get docking box coordinates for AutoDock Vina, Glide, etc.
        
        Args:
            pocket_id: Index of pocket
            padding: Extra space around pocket (Angstroms)
            
        Returns:
            dict: Docking box parameters
        """
        pocket = self.pockets[pocket_id]
        sphere_coords = np.array([sphere.coord for sphere in pocket.spheres])
        
        # Calculate bounding box
        min_coords = np.min(sphere_coords, axis=0) - padding
        max_coords = np.max(sphere_coords, axis=0) + padding
        
        center = (min_coords + max_coords) / 2
        size = max_coords - min_coords
        
        return {
            'center_x': float(center[0]),
            'center_y': float(center[1]),
            'center_z': float(center[2]),
            'size_x': float(size[0]),
            'size_y': float(size[1]),
            'size_z': float(size[2]),
            'min_x': float(min_coords[0]),
            'min_y': float(min_coords[1]),
            'min_z': float(min_coords[2]),
            'max_x': float(max_coords[0]),
            'max_y': float(max_coords[1]),
            'max_z': float(max_coords[2])
        }
    
    def get_pocket_residues(self, pocket_id=0, distance_cutoff=5.0):
        """
        Get residues lining the pocket (for cofolding constraints).
        
        Args:
            pocket_id: Index of pocket
            distance_cutoff: Distance from sphere centers to include residues (Å)
            
        Returns:
            list: List of (chain_id, res_num, res_name) tuples
        """
        pocket = self.pockets[pocket_id]
        sphere_coords = np.array([sphere.coord for sphere in pocket.spheres])
        
        # Get protein atom coordinates
        atom_coords = self.atomarray.coord
        
        # Find atoms within cutoff of any sphere center
        pocket_residues = set()
        
        for sphere_center in sphere_coords:
            # Calculate distances to all atoms
            distances = np.linalg.norm(atom_coords - sphere_center, axis=1)
            
            # Get atoms within cutoff
            nearby_indices = np.where(distances <= distance_cutoff)[0]
            
            # Extract residue information
            for idx in nearby_indices:
                chain_id = self.atomarray.chain_id[idx]
                res_num = self.atomarray.res_id[idx]
                res_name = self.atomarray.res_name[idx]
                pocket_residues.add((chain_id, res_num, res_name))
        
        # Sort by chain and residue number
        sorted_residues = sorted(list(pocket_residues), 
                                key=lambda x: (x[0], x[1]))
        
        return sorted_residues
    
    def get_pocket_info(self, pocket_id=0):
        """
        Get comprehensive pocket information.
        
        Args:
            pocket_id: Index of pocket
            
        Returns:
            dict: All pocket information
        """
        pocket = self.pockets[pocket_id]
        
        return {
            'pocket_id': pocket_id,
            'rank': pocket_id + 1,
            'score': float(pocket.score),
            'volume': float(pocket.volume),
            'n_spheres': pocket.n_spheres,
            'center': self.get_pocket_center(pocket_id).tolist(),
            'docking_box': self.get_docking_box(pocket_id),
            'residues': [
                {'chain': r[0], 'number': int(r[1]), 'name': r[2]}
                for r in self.get_pocket_residues(pocket_id)
            ]
        }
    
    def export_vina_config(self, pocket_id=0, output_file='vina_config.txt', 
                          exhaustiveness=8):
        """
        Export AutoDock Vina configuration file.
        
        Args:
            pocket_id: Index of pocket
            output_file: Output filename
            exhaustiveness: Vina exhaustiveness parameter
        """
        box = self.get_docking_box(pocket_id)
        
        config = f"""# AutoDock Vina configuration for pocket {pocket_id}
# Generated from Pocketeer output

receptor = receptor.pdbqt
ligand = ligand.pdbqt

center_x = {box['center_x']:.3f}
center_y = {box['center_y']:.3f}
center_z = {box['center_z']:.3f}

size_x = {box['size_x']:.1f}
size_y = {box['size_y']:.1f}
size_z = {box['size_z']:.1f}

exhaustiveness = {exhaustiveness}
"""
        
        with open(output_file, 'w') as f:
            f.write(config)
        
        print(f"Vina config written to {output_file}")
    
    def export_pymol_selection(self, pocket_id=0, output_file='pymol_pocket.pml'):
        """
        Export PyMOL script to visualize pocket residues.
        
        Args:
            pocket_id: Index of pocket
            output_file: Output filename
        """
        residues = self.get_pocket_residues(pocket_id)
        center = self.get_pocket_center(pocket_id)
        
        # Build selection string
        selections = []
        for chain, resnum, _ in residues:
            selections.append(f"(chain {chain} and resi {resnum})")
        
        selection_str = " or ".join(selections)
        
        script = f"""# PyMOL script to visualize pocket {pocket_id}
# Generated from Pocketeer output

# Select pocket residues
select pocket_{pocket_id}, {selection_str}

# Show as sticks
show sticks, pocket_{pocket_id}
color cyan, pocket_{pocket_id}

# Show pocket center
pseudoatom pocket_center_{pocket_id}, pos=[{center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}]
show spheres, pocket_center_{pocket_id}
color red, pocket_center_{pocket_id}

# Zoom to pocket
zoom pocket_{pocket_id}
"""
        
        with open(output_file, 'w') as f:
            f.write(script)
        
        print(f"PyMOL script written to {output_file}")
    
    def export_all_pockets_json(self, output_file='pockets_info.json'):
        """
        Export all pocket information to JSON.
        
        Args:
            output_file: Output filename
        """
        all_pockets = []
        
        for i in range(len(self.pockets)):
            all_pockets.append(self.get_pocket_info(i))
        
        with open(output_file, 'w') as f:
            json.dump(all_pockets, f, indent=2)
        
        print(f"All pocket info written to {output_file}")
    
    def print_summary(self):
        """Print a summary of all detected pockets."""
        print("\n" + "="*70)
        print("POCKET DETECTION SUMMARY")
        print("="*70)
        
        for i, pocket in enumerate(self.pockets):
            center = self.get_pocket_center(i)
            n_residues = len(self.get_pocket_residues(i))
            
            print(f"\nPocket {i} (Rank {i+1}):")
            print(f"  Score:      {pocket.score:.2f}")
            print(f"  Volume:     {pocket.volume:.1f} Å³")
            print(f"  Spheres:    {pocket.n_spheres}")
            print(f"  Center:     ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})")
            print(f"  Residues:   {n_residues}")


def main():
    """Example usage of PocketExtractor."""
    
    # Initialize
    protein_file = "protein.pdb"  # Replace with your protein file
    extractor = PocketExtractor(protein_file)
    
    # Find pockets (customize parameters as needed)
    pockets = extractor.find_pockets(
        r_min=3.0,
        r_max=6.0,
        min_spheres=30,
        merge_distance=2.5
    )
    
    # Print summary
    extractor.print_summary()
    
    # Export for top pocket (pocket 0)
    print("\n" + "="*70)
    print("EXPORTING TOP POCKET INFORMATION")
    print("="*70 + "\n")
    
    # 1. Get pocket center for cofolding
    center = extractor.get_pocket_center(pocket_id=0)
    print(f"Pocket center: {center}")
    
    # 2. Get docking box
    box = extractor.get_docking_box(pocket_id=0, padding=5.0)
    print(f"\nDocking box:")
    print(f"  Center: ({box['center_x']:.2f}, {box['center_y']:.2f}, {box['center_z']:.2f})")
    print(f"  Size:   ({box['size_x']:.2f}, {box['size_y']:.2f}, {box['size_z']:.2f})")
    
    # 3. Get pocket residues
    residues = extractor.get_pocket_residues(pocket_id=0, distance_cutoff=5.0)
    print(f"\nPocket residues ({len(residues)} total):")
    print("  First 10:", residues[:10])
    
    # 4. Export files
    extractor.export_vina_config(pocket_id=0, output_file='vina_pocket0.txt')
    extractor.export_pymol_selection(pocket_id=0, output_file='pymol_pocket0.pml')
    extractor.export_all_pockets_json(output_file='all_pockets.json')
    
    # Example: Export for multiple pockets (allosteric vs orthosteric)
    print("\n" + "="*70)
    print("EXPORTING MULTIPLE POCKETS (for allosteric/orthosteric)")
    print("="*70 + "\n")
    
    if len(pockets) >= 2:
        # Top pocket (likely orthosteric)
        extractor.export_vina_config(pocket_id=0, output_file='vina_orthosteric.txt')
        
        # Second pocket (potentially allosteric)
        extractor.export_vina_config(pocket_id=1, output_file='vina_allosteric.txt')
        
        print("Exported configs for orthosteric and allosteric sites")


if __name__ == "__main__":
    main()
