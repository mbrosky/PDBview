"""Color schemes for visualization"""

from enum import Enum
from typing import Tuple
import colorsys

from ..model.protein import Residue, SecondaryStructure


class ColorScheme(Enum):
    """Available color schemes"""
    SECONDARY_STRUCTURE = "secondary"
    CHAIN = "chain"
    ELEMENT = "element"
    BFACTOR = "bfactor"
    RAINBOW = "rainbow"


def get_color_for_residue(residue: Residue, scheme: ColorScheme, chain_idx: int = 0, total_residues: int = 1) -> Tuple[int, int, int]:
    """
    Get RGB color for a residue based on color scheme.

    Args:
        residue: Residue to color
        scheme: Color scheme to use
        chain_idx: Index of chain (for chain coloring)
        total_residues: Total number of residues (for rainbow)

    Returns:
        RGB tuple (0-255)
    """
    if scheme == ColorScheme.SECONDARY_STRUCTURE:
        return _color_by_secondary_structure(residue)
    elif scheme == ColorScheme.CHAIN:
        return _color_by_chain(chain_idx)
    elif scheme == ColorScheme.ELEMENT:
        # For residue-level, use a default or average atom color
        if residue.atoms:
            return get_element_color(residue.atoms[0].element)
        return (128, 128, 128)
    elif scheme == ColorScheme.BFACTOR:
        return _color_by_bfactor(residue)
    elif scheme == ColorScheme.RAINBOW:
        return _color_rainbow(residue.seq_num, total_residues)
    else:
        return (128, 128, 128)  # Gray default


def _color_by_secondary_structure(residue: Residue) -> Tuple[int, int, int]:
    """Color by secondary structure - VIVID COLORS for better contrast"""
    if residue.is_ligand:
        return (255, 165, 0)  # Bright orange for ligands

    # Much more vivid and saturated colors for clear differentiation
    ss_colors = {
        SecondaryStructure.HELIX: (255, 60, 60),     # Bright vivid red
        SecondaryStructure.SHEET: (255, 235, 60),    # Bright vivid yellow
        SecondaryStructure.COIL: (80, 255, 120),     # Bright vivid green
        SecondaryStructure.TURN: (100, 180, 255),    # Bright vivid blue
    }
    return ss_colors.get(residue.secondary_structure, (150, 150, 150))


def _color_by_chain(chain_idx: int) -> Tuple[int, int, int]:
    """Color by chain using distinct colors"""
    chain_palette = [
        (230, 50, 50),    # Red
        (50, 130, 230),   # Blue
        (80, 200, 120),   # Green
        (230, 180, 50),   # Yellow
        (200, 80, 200),   # Purple
        (255, 128, 0),    # Orange
        (0, 200, 200),    # Cyan
        (255, 100, 180),  # Pink
    ]
    return chain_palette[chain_idx % len(chain_palette)]


def _color_by_bfactor(residue: Residue) -> Tuple[int, int, int]:
    """Color by B-factor (blue=low, red=high)"""
    if not residue.atoms:
        return (128, 128, 128)

    avg_bfactor = sum(a.b_factor for a in residue.atoms) / len(residue.atoms)

    # Normalize to 0-1 (typical B-factors: 10-50)
    normalized = (avg_bfactor - 10) / 40
    normalized = max(0.0, min(1.0, normalized))

    # Blue to red gradient
    r = int(255 * normalized)
    g = 0
    b = int(255 * (1 - normalized))

    return (r, g, b)


def _color_rainbow(residue_num: int, total: int) -> Tuple[int, int, int]:
    """Rainbow color from N-terminus (blue) to C-terminus (red)"""
    if total <= 1:
        hue = 0.5
    else:
        hue = 0.66 - (residue_num / total) * 0.66  # 0.66 (blue) to 0 (red)

    rgb = colorsys.hsv_to_rgb(hue, 0.8, 0.9)
    return tuple(int(c * 255) for c in rgb)


def get_element_color(element: str) -> Tuple[int, int, int]:
    """Get CPK color for an element - VIVID COLORS for better visibility"""
    cpk_colors = {
        'C': (170, 170, 170),   # Lighter gray for contrast
        'N': (70, 130, 255),    # Brighter blue
        'O': (255, 50, 50),     # Brighter red
        'S': (255, 220, 80),    # Brighter yellow
        'P': (255, 150, 50),    # Brighter orange
        'H': (255, 255, 255),   # White
        'FE': (255, 120, 60),   # Brighter orange-brown (iron in heme)
        'ZN': (150, 160, 200),  # Brighter blue-gray
        'CA': (80, 255, 50),    # Brighter green
        'MG': (160, 255, 50),   # Brighter light green
        'CL': (50, 255, 50),    # Brighter green
        'NA': (200, 120, 255),  # Brighter purple
    }
    return cpk_colors.get(element.upper(), (220, 220, 220))


def rgb_to_hex(rgb: Tuple[int, int, int]) -> str:
    """Convert RGB tuple to hex color string"""
    return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"
