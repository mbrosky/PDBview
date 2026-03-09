"""Command-line interface"""

import argparse
import sys
from pathlib import Path

from .parser.pdb_parser import load_pdb_structure
from .parser.fetch import fetch_from_rcsb
from .ui.app import PDBViewApp


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="PDBView - Terminal protein and ligand structure viewer",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  pdbview protein.pdb              Load local PDB file
  pdbview structure.cif            Load mmCIF file
  pdbview --fetch 4HHB             Download and view from RCSB PDB
  pdbview protein.pdb --dssp       Use DSSP for secondary structure

Keybindings:
  h/l     Rotate left/right
  j/k     Rotate down/up
  u/i     Roll left/right
  +/-     Zoom in/out
  r       Reset view
  c       Cycle color schemes
  v       Toggle ligand visibility
  f       Toggle info sidebar
  space   Toggle auto-rotation
  q       Quit
        """
    )

    parser.add_argument(
        "structure",
        nargs="?",
        help="Path to PDB/mmCIF file or PDB ID (with --fetch)"
    )

    parser.add_argument(
        "--fetch",
        action="store_true",
        help="Download structure from RCSB PDB by ID"
    )

    parser.add_argument(
        "--dssp",
        action="store_true",
        help="Use DSSP for secondary structure assignment (requires dssp binary)"
    )

    parser.add_argument(
        "--include-water",
        action="store_true",
        help="Include water molecules (excluded by default)"
    )

    parser.add_argument(
        "--include-hydrogen",
        action="store_true",
        help="Include hydrogen atoms (excluded by default)"
    )

    parser.add_argument(
        "--version",
        action="version",
        version="PDBView 0.2.0"
    )

    args = parser.parse_args()

    # Check if structure provided
    if not args.structure:
        parser.print_help()
        sys.exit(1)

    try:
        # Fetch from RCSB if requested
        if args.fetch:
            print(f"Fetching {args.structure} from RCSB PDB...")
            filepath = fetch_from_rcsb(args.structure)
        else:
            filepath = args.structure

        # Check file exists
        if not Path(filepath).exists():
            print(f"Error: File not found: {filepath}", file=sys.stderr)
            # Give helpful hint if it looks like a PDB ID
            if len(args.structure) == 4 and args.structure.isalnum():
                print(f"\nDid you mean to use --fetch?", file=sys.stderr)
                print(f"  pdbview --fetch {args.structure}", file=sys.stderr)
            sys.exit(1)

        # Load structure
        print(f"Loading structure from {filepath}...")
        exclude_water = not args.include_water  # Invert: --include-water means exclude_water=False
        exclude_hydrogen = not args.include_hydrogen  # Invert: --include-hydrogen means exclude_hydrogen=False
        protein = load_pdb_structure(filepath, run_dssp=args.dssp, exclude_water=exclude_water, exclude_hydrogen=exclude_hydrogen)

        print(f"Loaded: {protein.name}")
        print(f"  Chains: {len(protein.chains)}")
        print(f"  Residues: {protein.n_residues}")
        print(f"  Atoms: {protein.n_atoms}")
        print(f"  Ligands: {protein.n_ligands}")
        print()

        # Launch TUI
        app = PDBViewApp(protein)
        app.run()

    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error loading structure: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
