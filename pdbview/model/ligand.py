"""Ligand detection and classification"""

from enum import Enum
from typing import Set


class ResidueType(Enum):
    """Classification of residue types"""
    STANDARD_AMINO_ACID = "amino_acid"
    LIGAND = "ligand"
    WATER = "water"
    ION = "ion"
    NUCLEOTIDE = "nucleotide"
    MODIFIED_RESIDUE = "modified"


# Standard 20 amino acids
STANDARD_AMINO_ACIDS: Set[str] = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS',
    'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO',
    'SER', 'THR', 'TRP', 'TYR', 'VAL',
}

# Common water molecule names
WATER_NAMES: Set[str] = {'HOH', 'WAT', 'H2O', 'DOD', 'TIP3'}

# Common ions
ION_NAMES: Set[str] = {
    'NA', 'CL', 'CA', 'MG', 'ZN', 'FE', 'K', 'MN', 'CO', 'NI', 'CU',
    'NA+', 'CL-', 'CA2+', 'MG2+', 'ZN2+', 'FE2+', 'FE3+', 'K+',
}

# DNA/RNA nucleotides
NUCLEOTIDES: Set[str] = {
    'A', 'C', 'G', 'T', 'U',  # Single letter
    'DA', 'DC', 'DG', 'DT',    # DNA
    'ADE', 'CYT', 'GUA', 'THY', 'URA',  # Full names
}

# Common modified amino acids
MODIFIED_RESIDUES: Set[str] = {
    'MSE', 'SEP', 'TPO', 'PTR', 'HYP', 'MLY', 'ALY',
}

# Common ligands and cofactors
KNOWN_LIGANDS: Set[str] = {
    'HEM', 'HEME',  # Heme
    'ATP', 'ADP', 'AMP', 'GTP', 'GDP',  # Nucleotides
    'NAD', 'NAP', 'FAD', 'FMN',  # Cofactors
    'COA', 'ACO',  # Coenzyme A
    'PLP', 'B12',  # Vitamins
    'SAM', 'SAH',  # Methylation
}


def classify_residue(residue_name: str) -> ResidueType:
    """
    Classify a residue by its 3-letter code.

    Args:
        residue_name: Three-letter residue code (e.g., 'ALA', 'HEM', 'HOH')

    Returns:
        ResidueType classification
    """
    name = residue_name.strip().upper()

    # Standard amino acids
    if name in STANDARD_AMINO_ACIDS:
        return ResidueType.STANDARD_AMINO_ACID

    # Water molecules
    if name in WATER_NAMES:
        return ResidueType.WATER

    # Ions
    if name in ION_NAMES:
        return ResidueType.ION

    # Nucleotides
    if name in NUCLEOTIDES:
        return ResidueType.NUCLEOTIDE

    # Modified residues
    if name in MODIFIED_RESIDUES:
        return ResidueType.MODIFIED_RESIDUE

    # Everything else is a ligand
    return ResidueType.LIGAND


def is_ligand(residue_name: str) -> bool:
    """
    Check if a residue is a ligand (not protein or water).
    Includes ions, nucleotides, and other heteroatoms.

    Args:
        residue_name: Three-letter residue code

    Returns:
        True if residue is classified as a ligand (including ions)
    """
    res_type = classify_residue(residue_name)
    return res_type in (ResidueType.LIGAND, ResidueType.NUCLEOTIDE, ResidueType.ION)


def is_water(residue_name: str) -> bool:
    """Check if residue is water"""
    return classify_residue(residue_name) == ResidueType.WATER


def is_ion(residue_name: str) -> bool:
    """Check if residue is an ion"""
    return classify_residue(residue_name) == ResidueType.ION


def get_ligand_info(residue_name: str) -> dict:
    """
    Get detailed information about a ligand.

    Args:
        residue_name: Three-letter residue code

    Returns:
        Dictionary with ligand information
    """
    res_type = classify_residue(residue_name)

    info = {
        'name': residue_name,
        'type': res_type.value,
        'is_ligand': is_ligand(residue_name),
        'is_cofactor': residue_name in KNOWN_LIGANDS,
    }

    # Add common names for known ligands
    common_names = {
        'HEM': 'Heme',
        'ATP': 'Adenosine Triphosphate',
        'ADP': 'Adenosine Diphosphate',
        'NAD': 'Nicotinamide Adenine Dinucleotide',
        'FAD': 'Flavin Adenine Dinucleotide',
        'COA': 'Coenzyme A',
    }

    if residue_name in common_names:
        info['common_name'] = common_names[residue_name]

    return info
