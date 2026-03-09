"""Protein structure data models"""

from dataclasses import dataclass, field
from enum import Enum
from typing import List, Tuple
import numpy as np


class SecondaryStructure(Enum):
    """Secondary structure classification"""
    HELIX = "helix"
    SHEET = "sheet"
    TURN = "turn"
    COIL = "coil"


@dataclass
class Atom:
    """Individual atom in a structure"""
    name: str
    element: str
    x: float
    y: float
    z: float
    b_factor: float = 0.0
    occupancy: float = 1.0
    is_ca: bool = False

    @property
    def coords(self) -> np.ndarray:
        """Return coordinates as numpy array"""
        return np.array([self.x, self.y, self.z])

    def distance_to(self, other: 'Atom') -> float:
        """Calculate distance to another atom"""
        dx = self.x - other.x
        dy = self.y - other.y
        dz = self.z - other.z
        return np.sqrt(dx*dx + dy*dy + dz*dz)


@dataclass
class Residue:
    """Amino acid residue or ligand"""
    name: str
    seq_num: int
    atoms: List[Atom] = field(default_factory=list)
    secondary_structure: SecondaryStructure = SecondaryStructure.COIL
    is_ligand: bool = False
    is_water: bool = False
    chain_id: str = ""

    @property
    def ca_atom(self) -> Atom | None:
        """Get alpha carbon atom if present"""
        for atom in self.atoms:
            if atom.is_ca:
                return atom
        return None

    @property
    def center(self) -> np.ndarray:
        """Calculate geometric center of residue"""
        if not self.atoms:
            return np.array([0.0, 0.0, 0.0])
        coords = np.array([[a.x, a.y, a.z] for a in self.atoms])
        return coords.mean(axis=0)

    def min_distance_to(self, other: 'Residue') -> float:
        """Calculate minimum distance to another residue"""
        min_dist = float('inf')
        for a1 in self.atoms:
            for a2 in other.atoms:
                dist = a1.distance_to(a2)
                if dist < min_dist:
                    min_dist = dist
        return min_dist


@dataclass
class Chain:
    """Polypeptide chain or ligand collection"""
    id: str
    residues: List[Residue] = field(default_factory=list)

    @property
    def sequence(self) -> str:
        """Get single-letter amino acid sequence"""
        AA_MAP = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        }
        return ''.join(AA_MAP.get(r.name, 'X') for r in self.residues if not r.is_ligand)

    @property
    def n_residues(self) -> int:
        """Count protein residues only (excluding ligands)"""
        return sum(1 for r in self.residues if not r.is_ligand)

    @property
    def n_atoms(self) -> int:
        return sum(len(r.atoms) for r in self.residues)


@dataclass
class Protein:
    """Complete protein structure"""
    name: str
    chains: List[Chain] = field(default_factory=list)
    metadata: 'PDBMetadata' = None  # Optional metadata from PDB header

    @property
    def all_atoms(self) -> List[Atom]:
        """Get all atoms in structure"""
        atoms = []
        for chain in self.chains:
            for residue in chain.residues:
                atoms.extend(residue.atoms)
        return atoms

    @property
    def all_residues(self) -> List[Residue]:
        """Get all residues in structure"""
        residues = []
        for chain in self.chains:
            residues.extend(chain.residues)
        return residues

    @property
    def ca_atoms(self) -> List[Tuple[Atom, Residue, Chain]]:
        """Get all C-alpha atoms with context"""
        cas = []
        for chain in self.chains:
            for residue in chain.residues:
                if residue.ca_atom:
                    cas.append((residue.ca_atom, residue, chain))
        return cas

    @property
    def ligands(self) -> List[Tuple[Residue, Chain]]:
        """Get all ligand residues with their chain"""
        ligs = []
        for chain in self.chains:
            for residue in chain.residues:
                if residue.is_ligand:
                    ligs.append((residue, chain))
        return ligs

    @property
    def n_atoms(self) -> int:
        return len(self.all_atoms)

    @property
    def n_residues(self) -> int:
        """Count protein residues only (excluding ligands)"""
        return sum(1 for res in self.all_residues if not res.is_ligand)

    @property
    def n_ligands(self) -> int:
        return len(self.ligands)

    def center(self) -> np.ndarray:
        """Calculate geometric center"""
        atoms = self.all_atoms
        if not atoms:
            return np.array([0.0, 0.0, 0.0])
        coords = np.array([[a.x, a.y, a.z] for a in atoms])
        return coords.mean(axis=0)

    def center_at_origin(self):
        """Translate structure to origin"""
        center = self.center()
        for chain in self.chains:
            for residue in chain.residues:
                for atom in residue.atoms:
                    atom.x -= center[0]
                    atom.y -= center[1]
                    atom.z -= center[2]

    def bounding_radius(self) -> float:
        """Get bounding sphere radius from origin"""
        atoms = self.all_atoms
        if not atoms:
            return 0.0
        max_dist = 0.0
        for atom in atoms:
            dist = np.sqrt(atom.x**2 + atom.y**2 + atom.z**2)
            if dist > max_dist:
                max_dist = dist
        return max_dist
