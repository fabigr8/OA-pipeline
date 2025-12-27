#!/usr/bin/env python3
"""
Minimalistic PDB file downloader
Downloads PDB files from RCSB PDB database using PDB IDs
"""

import requests
from pathlib import Path
import pandas as pd

def download_pdb(pdb_id, output_dir="."):
    """
    Download a PDB file from RCSB PDB database
    
    Args:
        pdb_id: PDB identifier (e.g., '1ABC')
        output_dir: Directory to save the file (default: current directory)
    """
    pdb_id = pdb_id.lower()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    output_path = Path(output_dir) / f"{pdb_id}.pdb"
    
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        
        output_path.write_text(response.text)
        print(f"✓ Downloaded {pdb_id}.pdb")
        
    except requests.exceptions.RequestException as e:
        print(f"✗ Failed to download {pdb_id}: {e}")


def main():
    # read csv file containing pdb ids
    df = pd.read_csv('data/processed/OSAS_metadata.csv')
    
    # read first coplumn as list
    pdb_ids = df['AS PDB ID'].tolist()
    # split entry that contains multiple pdb ids separated by ;
    pdb_ids = [pdb_id.strip() for entry in pdb_ids for pdb_id in entry.split(';')]
    print(f"Initial PDB IDs from first column: {len(pdb_ids)}")
    
    #read second column as list 
    pdb_ids2 = df['OS PDB ID'].tolist()
    # split entry that contains multiple pdb ids separated by ;
    pdb_ids2 = [pdb_ids2.strip() for entry in pdb_ids2 for pdb_ids2 in entry.split(';')]
    print(f"Initial PDB IDs from first column: {len(pdb_ids2)}")
    
    # Merge both lists
    pdb_ids.extend(pdb_ids2)
    
    # Optional: specify output directory
    output_dir = "data/raw"
    Path(output_dir).mkdir(exist_ok=True)
    
    # Download all PDB files
    for pdb_id in pdb_ids:
        # download with progress indicator
        download_pdb(pdb_id, output_dir)


if __name__ == "__main__":
    main()