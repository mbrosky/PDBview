"""Binding site detection"""

from dataclasses import dataclass, field
from typing import List, Tuple
import numpy as np

from .protein import Protein, Residue, Chain


# Distance cutoff for binding site detection (Angstroms)
BINDING_CUTOFF = 5.0


@dataclass
class BindingResidue:
    """A protein residue in contact with a ligand"""
    chain_id: str
    res_name: str
    res_num: int
    min_distance: float
    contact_atoms: int = 0

    def __str__(self) -> str:
        return f"{self.res_name}{self.res_num} ({self.chain_id}) @ {self.min_distance:.2f}Å"


@dataclass
class BindingSite:
    """A protein-ligand binding site"""
    ligand_name: str
    ligand_chain: str
    ligand_resnum: int
    binding_residues: List[BindingResidue] = field(default_factory=list)

    @property
    def n_contacts(self) -> int:
        """Number of residues in contact"""
        return len(self.binding_residues)

    @property
    def closest_residue(self) -> BindingResidue | None:
        """Get closest residue to ligand"""
        if not self.binding_residues:
            return None
        return min(self.binding_residues, key=lambda r: r.min_distance)

    def get_top_residues(self, n: int = 5) -> List[BindingResidue]:
        """Get N closest residues"""
        sorted_res = sorted(self.binding_residues, key=lambda r: r.min_distance)
        return sorted_res[:n]

    def __str__(self) -> str:
        return f"Binding site: {self.ligand_name} (chain {self.ligand_chain}) - {self.n_contacts} contacts"


def detect_binding_sites(
    protein: Protein,
    cutoff: float = BINDING_CUTOFF,
    exclude_water: bool = True,
) -> List[BindingSite]:
    """
    Detect protein-ligand binding sites.

    Args:
        protein: Protein structure with ligands
        cutoff: Distance cutoff in Angstroms (default 5.0)
        exclude_water: Whether to exclude water molecules (default True)

    Returns:
        List of detected binding sites
    """
    binding_sites = []

    # Get all ligands
    ligands = protein.ligands
    if not ligands:
        return binding_sites

    # For each ligand, find nearby protein residues
    for ligand_res, ligand_chain in ligands:
        # Skip water if requested
        if exclude_water and ligand_res.is_water:
            continue

        binding_residues = []

        # Check all protein residues
        for chain in protein.chains:
            for prot_res in chain.residues:
                # Skip ligands and water
                if prot_res.is_ligand or prot_res.is_water:
                    continue

                # Calculate minimum distance between residues
                min_dist = _min_distance_between_residues(ligand_res, prot_res)

                # If within cutoff, add to binding site
                if min_dist < cutoff:
                    # Count atoms in contact (within cutoff)
                    contact_atoms = _count_contact_atoms(ligand_res, prot_res, cutoff)

                    binding_residues.append(BindingResidue(
                        chain_id=chain.id,
                        res_name=prot_res.name,
                        res_num=prot_res.seq_num,
                        min_distance=min_dist,
                        contact_atoms=contact_atoms,
                    ))

        # Create binding site if contacts found
        if binding_residues:
            site = BindingSite(
                ligand_name=ligand_res.name,
                ligand_chain=ligand_chain.id,
                ligand_resnum=ligand_res.seq_num,
                binding_residues=binding_residues,
            )
            binding_sites.append(site)

    return binding_sites


def _min_distance_between_residues(res1: Residue, res2: Residue) -> float:
    """Calculate minimum atom-atom distance between two residues"""
    min_dist = float('inf')

    for atom1 in res1.atoms:
        for atom2 in res2.atoms:
            dist = atom1.distance_to(atom2)
            if dist < min_dist:
                min_dist = dist

    return min_dist


def _count_contact_atoms(res1: Residue, res2: Residue, cutoff: float) -> int:
    """Count number of atom pairs within cutoff distance"""
    count = 0

    for atom1 in res1.atoms:
        for atom2 in res2.atoms:
            if atom1.distance_to(atom2) < cutoff:
                count += 1

    return count


def analyze_binding_site(site: BindingSite) -> dict:
    """
    Analyze a binding site and return statistics.

    Args:
        site: BindingSite to analyze

    Returns:
        Dictionary with analysis results
    """
    if not site.binding_residues:
        return {
            'n_contacts': 0,
            'avg_distance': 0.0,
            'closest_distance': 0.0,
        }

    distances = [r.min_distance for r in site.binding_residues]

    # Classify residues by properties
    hydrophobic = {'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'PRO'}
    polar = {'SER', 'THR', 'CYS', 'TYR', 'ASN', 'GLN'}
    charged_pos = {'LYS', 'ARG', 'HIS'}
    charged_neg = {'ASP', 'GLU'}

    n_hydrophobic = sum(1 for r in site.binding_residues if r.res_name in hydrophobic)
    n_polar = sum(1 for r in site.binding_residues if r.res_name in polar)
    n_charged_pos = sum(1 for r in site.binding_residues if r.res_name in charged_pos)
    n_charged_neg = sum(1 for r in site.binding_residues if r.res_name in charged_neg)

    return {
        'n_contacts': len(site.binding_residues),
        'avg_distance': np.mean(distances),
        'min_distance': min(distances),
        'max_distance': max(distances),
        'n_hydrophobic': n_hydrophobic,
        'n_polar': n_polar,
        'n_charged_positive': n_charged_pos,
        'n_charged_negative': n_charged_neg,
        'top_residues': [str(r) for r in site.get_top_residues(5)],
    }
