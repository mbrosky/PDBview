"""PDB/mmCIF parser using BioPython"""

from pathlib import Path
from Bio.PDB import PDBParser, MMCIFParser, DSSP
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import warnings

from ..model.protein import Protein, Chain, Residue, Atom, SecondaryStructure
from ..model.ligand import classify_residue, ResidueType, is_ligand
from .pdb_metadata import parse_pdb_metadata


# Suppress BioPython warnings
warnings.filterwarnings('ignore', category=UserWarning, module='Bio.PDB')


def load_pdb_structure(filepath: str, run_dssp: bool = False, exclude_water: bool = True, exclude_hydrogen: bool = True) -> Protein:
    """
    Load a protein structure from PDB or mmCIF file.

    Args:
        filepath: Path to PDB or mmCIF file
        run_dssp: Whether to run DSSP for secondary structure (requires dssp binary)
        exclude_water: Whether to exclude water molecules (default: True)
        exclude_hydrogen: Whether to exclude hydrogen atoms (default: True)

    Returns:
        Protein object with structure data
    """
    path = Path(filepath)

    if not path.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    # Choose parser based on extension
    ext = path.suffix.lower()
    if ext in ['.cif', '.mmcif']:
        parser = MMCIFParser(QUIET=True)
    elif ext in ['.pdb', '.ent']:
        parser = PDBParser(QUIET=True)
    else:
        # Try PDB parser as default
        parser = PDBParser(QUIET=True)

    # Parse structure
    structure = parser.get_structure(path.stem, str(path))

    # Convert to our data model
    protein = _convert_biopython_structure(structure, path.stem, exclude_water=exclude_water, exclude_hydrogen=exclude_hydrogen)

    # Parse and attach metadata
    protein.metadata = parse_pdb_metadata(filepath)

    # Assign secondary structure from PDB HELIX/SHEET records
    _assign_secondary_structure_from_pdb(protein, filepath)

    # Optionally run DSSP for more accurate secondary structure
    if run_dssp:
        try:
            _assign_secondary_structure_from_dssp(protein, filepath)
        except Exception as e:
            print(f"Warning: DSSP failed: {e}")

    return protein


def _convert_biopython_structure(bio_structure, name: str, exclude_water: bool = True, exclude_hydrogen: bool = True) -> Protein:
    """Convert BioPython Structure to our Protein model"""
    chains = []

    # Usually we just take the first model
    if len(bio_structure) == 0:
        raise ValueError("No models found in structure")

    model = bio_structure[0]

    for bio_chain in model:
        chain_id = bio_chain.id
        residues = []

        for bio_residue in bio_chain:
            # Get residue info
            hetero_flag, seq_num, insertion_code = bio_residue.id
            res_name = bio_residue.resname.strip()

            # Classify residue type
            res_type = classify_residue(res_name)

            # Skip water if requested
            if exclude_water and res_type == ResidueType.WATER:
                continue

            # Convert atoms
            atoms = []
            for bio_atom in bio_residue:
                # Get element
                element = bio_atom.element.strip() if hasattr(bio_atom, 'element') else bio_atom.name[0]

                # Skip hydrogen if requested
                if exclude_hydrogen and element.upper() == 'H':
                    continue

                atom = Atom(
                    name=bio_atom.name,
                    element=element,
                    x=float(bio_atom.coord[0]),
                    y=float(bio_atom.coord[1]),
                    z=float(bio_atom.coord[2]),
                    b_factor=float(bio_atom.bfactor) if hasattr(bio_atom, 'bfactor') else 0.0,
                    occupancy=float(bio_atom.occupancy) if hasattr(bio_atom, 'occupancy') else 1.0,
                    is_ca=(bio_atom.name == 'CA'),
                )
                atoms.append(atom)

            # Create residue
            residue = Residue(
                name=res_name,
                seq_num=seq_num,
                atoms=atoms,
                secondary_structure=SecondaryStructure.COIL,
                is_ligand=is_ligand(res_name),  # Use is_ligand() which includes ions
                is_water=(res_type == ResidueType.WATER),
                chain_id=chain_id,
            )
            residues.append(residue)

        # Create chain (only if it has residues after filtering)
        if residues:
            chain = Chain(id=chain_id, residues=residues)
            chains.append(chain)

    return Protein(name=name, chains=chains)


def _assign_secondary_structure_from_pdb(protein: Protein, filepath: str):
    """Assign secondary structure from PDB HELIX/SHEET records"""
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()

        # Parse HELIX records
        helices = {}
        for line in lines:
            if line.startswith('HELIX'):
                try:
                    chain_id = line[19].strip()
                    start_seq = int(line[21:25].strip())
                    end_seq = int(line[33:37].strip())

                    if chain_id not in helices:
                        helices[chain_id] = []
                    helices[chain_id].append((start_seq, end_seq))
                except (ValueError, IndexError):
                    continue

        # Parse SHEET records
        sheets = {}
        for line in lines:
            if line.startswith('SHEET'):
                try:
                    chain_id = line[21].strip()
                    start_seq = int(line[22:26].strip())
                    end_seq = int(line[33:37].strip())

                    if chain_id not in sheets:
                        sheets[chain_id] = []
                    sheets[chain_id].append((start_seq, end_seq))
                except (ValueError, IndexError):
                    continue

        # Assign to residues
        for chain in protein.chains:
            # Helices
            if chain.id in helices:
                for start, end in helices[chain.id]:
                    for residue in chain.residues:
                        if start <= residue.seq_num <= end and not residue.is_ligand:
                            residue.secondary_structure = SecondaryStructure.HELIX

            # Sheets
            if chain.id in sheets:
                for start, end in sheets[chain.id]:
                    for residue in chain.residues:
                        if start <= residue.seq_num <= end and not residue.is_ligand:
                            residue.secondary_structure = SecondaryStructure.SHEET

    except Exception as e:
        # Silently fail if we can't read secondary structure
        pass


def _assign_secondary_structure_from_dssp(protein: Protein, filepath: str):
    """Assign secondary structure using DSSP (requires dssp binary installed)"""
    try:
        dssp_dict = dssp_dict_from_pdb_file(filepath)

        # Map DSSP codes to our SecondaryStructure enum
        dssp_map = {
            'H': SecondaryStructure.HELIX,  # Alpha helix
            'G': SecondaryStructure.HELIX,  # 3-10 helix
            'I': SecondaryStructure.HELIX,  # Pi helix
            'E': SecondaryStructure.SHEET,  # Beta sheet
            'B': SecondaryStructure.SHEET,  # Beta bridge
            'T': SecondaryStructure.TURN,   # Turn
            'S': SecondaryStructure.COIL,   # Bend
            '-': SecondaryStructure.COIL,   # Coil
        }

        # Assign to residues
        for chain in protein.chains:
            for residue in chain.residues:
                if residue.is_ligand:
                    continue

                key = (chain.id, (' ', residue.seq_num, ' '))
                if key in dssp_dict:
                    dssp_code = dssp_dict[key][2]  # Secondary structure code
                    residue.secondary_structure = dssp_map.get(dssp_code, SecondaryStructure.COIL)

    except Exception:
        # DSSP not available or failed
        pass
