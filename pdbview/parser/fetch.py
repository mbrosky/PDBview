"""Fetch structures from RCSB PDB"""

import requests
from pathlib import Path
import tempfile


RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download"


def fetch_from_rcsb(pdb_id: str, format: str = "pdb") -> str:
    """
    Download a structure from RCSB PDB.

    Args:
        pdb_id: PDB identifier (e.g., '1UBQ', '4HHB')
        format: File format - 'pdb' or 'cif' (default: 'pdb')

    Returns:
        Path to downloaded file

    Raises:
        ValueError: If PDB ID is invalid
        requests.HTTPError: If download fails
    """
    pdb_id = pdb_id.strip().upper()

    if len(pdb_id) != 4:
        raise ValueError(f"Invalid PDB ID: {pdb_id}. Must be 4 characters.")

    # Determine file extension
    if format.lower() == "cif":
        ext = "cif"
    else:
        ext = "pdb"

    # Construct download URL
    url = f"{RCSB_DOWNLOAD_URL}/{pdb_id}.{ext}"

    # Download file
    print(f"Downloading {pdb_id} from RCSB PDB...")
    response = requests.get(url, timeout=30)
    response.raise_for_status()

    # Save to temporary file
    temp_dir = Path(tempfile.gettempdir())
    output_path = temp_dir / f"{pdb_id}.{ext}"

    with open(output_path, 'w') as f:
        f.write(response.text)

    print(f"Downloaded to: {output_path}")
    return str(output_path)


def fetch_ligand_info(ligand_id: str) -> dict:
    """
    Fetch information about a ligand from RCSB PDB.

    Args:
        ligand_id: Three-letter ligand code (e.g., 'HEM', 'ATP')

    Returns:
        Dictionary with ligand information
    """
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_id}"

    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        return {'error': str(e), 'id': ligand_id}
