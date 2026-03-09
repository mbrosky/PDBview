# PDBView

```
██████╗ ██████╗ ██████╗ ██╗   ██╗██╗███████╗██╗    ██╗
██╔══██╗██╔══██╗██╔══██╗██║   ██║██║██╔════╝██║    ██║
██████╔╝██║  ██║██████╔╝██║   ██║██║█████╗  ██║ █╗ ██║
██╔═══╝ ██║  ██║██╔══██╗╚██╗ ██╔╝██║██╔══╝  ██║███╗██║
██║     ██████╔╝██████╔╝ ╚████╔╝ ██║███████╗╚███╔███╔╝
╚═╝     ╚═════╝ ╚═════╝   ╚═══╝  ╚═╝╚══════╝ ╚══╝╚══╝
```

Visualize protein structures on remote servers via SSH - no GUI needed.


## Features

  - View protein structures in your terminal
  - Interactive 3D rotation and zoom
  - Multiple color schemes
  - Download from RCSB PDB database
  - Ligand visualization

## Installation

```bash
# Clone repository
git clone https://github.com/BrikiMondher/PDBView.git
cd PDBView

# Install with pipx
pipx install .

# Verify installation
pdbview --version
```

## Usage

```bash
# Fetch from RCSB PDB
pdbview --fetch 4HHB

# Load local file
pdbview structure.pdb
```

## Controls

| Key | Action |
|-----|--------|
| `h` / `l` | Rotate left/right |
| `j` / `k` | Rotate up/down |
| `u` / `i` | Roll left/right |
| `+` / `-` | Zoom in/out |
| `r` | Reset view |
| `c` | Cycle color schemes |
| `v` | Toggle ligands |
| `f` | Toggle sidebar |
| `space` | Auto-rotate |
| `q` | Quit |


## Requirements

- Python 3.9+
- Linux terminal with Unicode support
- Dependencies: textual, biopython, numpy, requests

## License

MIT License - see [LICENSE](LICENSE) file


## Email: mondher.briki.isdd@gmail.com
