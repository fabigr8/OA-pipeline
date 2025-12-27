#!/usr/bin/env python3
"""
Minimalistic UniProt FASTA downloader
Downloads FASTA sequences from UniProt using IDs from a CSV column
"""

import requests
import pandas as pd
from pathlib import Path


def download_fasta(uniprot_id, output_dir="."):
    """
    Download a FASTA file from UniProt
    
    Args:
        uniprot_id: UniProt accession (e.g., 'P12345')
        output_dir: Directory to save the file (default: current directory)
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    output_path = Path(output_dir) / f"{uniprot_id}.fasta"
    
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        
        output_path.write_text(response.text)
        print(f"✓ Downloaded {uniprot_id}.fasta")
        
    except requests.exceptions.RequestException as e:
        print(f"✗ Failed to download {uniprot_id}: {e}")


def main():
    # Configuration
    csv_file = "data/processed/OSAS_metadata.csv"  # Your CSV file
    column_name = "uniprotID"  # Column containing UniProt IDs
    output_dir = "data/raw/fastas"  # Output directory for FASTA files
    
    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)
    
    # Read CSV
    df = pd.read_csv(csv_file)
    
    # Get UniProt IDs from specified column
    uniprot_ids = df[column_name].dropna().unique()
    
    print(f"Found {len(uniprot_ids)} unique UniProt IDs")
    
    # Download all FASTA files for the UniProt IDs
    for uniprot_id in uniprot_ids:
        download_fasta(str(uniprot_id).strip(), output_dir)


if __name__ == "__main__":
    main()