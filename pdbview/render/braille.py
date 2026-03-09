"""Braille character rendering for terminal"""

import numpy as np
from typing import List, Tuple
from ..model.protein import Protein, Atom, SecondaryStructure
from .camera import Camera
from .colors import ColorScheme, get_color_for_residue, get_element_color, rgb_to_hex


# Braille Unicode offset
BRAILLE_OFFSET = 0x2800

# Braille dot positions (2x4 grid)
# ⠁ ⠂ ⠄ ⠈
# ⠐ ⠠ ⡀ ⢀
BRAILLE_DOTS = [
    0x01, 0x02, 0x04, 0x40,  # Left column (top to bottom)
    0x08, 0x10, 0x20, 0x80,  # Right column (top to bottom)
]


class BrailleRenderer:
    """Renders 3D structures using Unicode Braille characters"""

    def __init__(self, width: int = 80, height: int = 40):
        """
        Initialize renderer.

        Args:
            width: Width in braille characters
            height: Height in braille characters
        """
        self.width = width
        self.height = height
        # Each braille char is 2x4 dots
        self.pixel_width = width * 2
        self.pixel_height = height * 4

        # Framebuffer for braille dots
        self.buffer = np.zeros((self.pixel_height, self.pixel_width), dtype=bool)
        # Color buffer (stores RGB for each pixel)
        self.color_buffer = np.zeros((self.pixel_height, self.pixel_width, 3), dtype=np.uint8)
        # Depth buffer for proper occlusion
        self.depth_buffer = np.full((self.pixel_height, self.pixel_width), float('inf'))

    def clear(self):
        """Clear the framebuffer"""
        self.buffer.fill(False)
        self.color_buffer.fill(0)
        self.depth_buffer.fill(float('inf'))

    def draw_point(self, x: int, y: int, z: float, color: Tuple[int, int, int]):
        """Draw a single pixel with depth testing"""
        if 0 <= x < self.pixel_width and 0 <= y < self.pixel_height:
            if z < self.depth_buffer[y, x]:
                self.buffer[y, x] = True
                self.color_buffer[y, x] = color
                self.depth_buffer[y, x] = z

    def draw_line(self, x0: int, y0: int, z0: float, x1: int, y1: int, z1: float, color: Tuple[int, int, int]):
        """Draw a line using Bresenham's algorithm"""
        dx = abs(x1 - x0)
        dy = abs(y1 - y0)
        sx = 1 if x0 < x1 else -1
        sy = 1 if y0 < y1 else -1
        err = dx - dy

        x, y = x0, y0
        steps = max(dx, dy)
        if steps == 0:
            self.draw_point(x, y, z0, color)
            return

        for i in range(steps + 1):
            # Interpolate Z
            t = i / steps if steps > 0 else 0
            z = z0 + (z1 - z0) * t
            self.draw_point(x, y, z, color)

            if x == x1 and y == y1:
                break

            e2 = 2 * err
            if e2 > -dy:
                err -= dy
                x += sx
            if e2 < dx:
                err += dx
                y += sy

    def draw_circle(self, cx: int, cy: int, radius: int, z: float, color: Tuple[int, int, int], filled: bool = True):
        """Draw a circle (for atoms in ball-and-stick)"""
        if filled:
            for y in range(-radius, radius + 1):
                for x in range(-radius, radius + 1):
                    if x*x + y*y <= radius*radius:
                        self.draw_point(cx + x, cy + y, z, color)
        else:
            # Draw circle outline
            for angle in range(0, 360, 5):
                rad = np.radians(angle)
                x = int(cx + radius * np.cos(rad))
                y = int(cy + radius * np.sin(rad))
                self.draw_point(x, y, z, color)

    def render_to_string(self) -> str:
        """Convert framebuffer to colored braille string"""
        lines = []

        for row in range(self.height):
            line_chars = []
            for col in range(self.width):
                # Get 2x4 block of pixels
                y_base = row * 4
                x_base = col * 2

                # Calculate braille character
                char_code = BRAILLE_OFFSET
                colors_in_cell = []

                # Map pixels to braille dots
                for dot_idx in range(8):
                    if dot_idx < 4:  # Left column
                        px = x_base
                        py = y_base + dot_idx
                    else:  # Right column
                        px = x_base + 1
                        py = y_base + (dot_idx - 4)

                    if py < self.pixel_height and px < self.pixel_width:
                        if self.buffer[py, px]:
                            char_code |= BRAILLE_DOTS[dot_idx]
                            colors_in_cell.append(tuple(self.color_buffer[py, px]))

                # Choose dominant color
                # Use Rich markup instead of ANSI codes for Textual
                if colors_in_cell:
                    color = colors_in_cell[len(colors_in_cell) // 2]  # Median color
                    # Rich markup format: [rgb(R,G,B)]char[/]
                    line_chars.append(f"[rgb({color[0]},{color[1]},{color[2]})]{chr(char_code)}[/]")
                else:
                    line_chars.append(chr(char_code))

            lines.append(''.join(line_chars))

        return '\n'.join(lines)


def apply_depth_fog(color: Tuple[int, int, int], z_depth: float, max_depth: float) -> Tuple[int, int, int]:
    """
    Apply depth fog to a color based on Z depth.

    Args:
        color: Base RGB color (0-255)
        z_depth: Z coordinate (closer = negative, farther = positive)
        max_depth: Maximum depth for full fog

    Returns:
        Fogged RGB color
    """
    # Normalize depth (0 = close, 1 = far)
    fog_factor = min(1.0, max(0.0, abs(z_depth) / max_depth))

    # Apply exponential fog for better depth perception
    fog_factor = fog_factor ** 1.5

    # Background color (dark gray, not pure black)
    bg = (25, 25, 30)

    # Interpolate between original color and background
    r = int(color[0] * (1 - fog_factor) + bg[0] * fog_factor)
    g = int(color[1] * (1 - fog_factor) + bg[1] * fog_factor)
    b = int(color[2] * (1 - fog_factor) + bg[2] * fog_factor)

    return (r, g, b)


def draw_thick_line(
    renderer: 'BrailleRenderer',
    x0: int, y0: int, z0: float,
    x1: int, y1: int, z1: float,
    color: Tuple[int, int, int],
    thickness: int = 2
):
    """
    Draw a thick line by drawing multiple parallel lines.

    Args:
        renderer: BrailleRenderer instance
        x0, y0, z0: Start point
        x1, y1, z1: End point
        color: RGB color
        thickness: Line thickness in pixels (default: 2)
    """
    # Draw center line
    renderer.draw_line(x0, y0, z0, x1, y1, z1, color)

    # Draw parallel lines for thickness
    for offset in range(1, thickness + 1):
        # Horizontal offsets
        renderer.draw_line(x0 + offset, y0, z0, x1 + offset, y1, z1, color)
        renderer.draw_line(x0 - offset, y0, z0, x1 - offset, y1, z1, color)

        # Vertical offsets
        renderer.draw_line(x0, y0 + offset, z0, x1, y1 + offset, z1, color)
        renderer.draw_line(x0, y0 - offset, z0, x1, y1 - offset, z1, color)

        # Diagonal offsets for smoother look
        renderer.draw_line(x0 + offset, y0 + offset, z0, x1 + offset, y1 + offset, z1, color)
        renderer.draw_line(x0 - offset, y0 - offset, z0, x1 - offset, y1 - offset, z1, color)


def render_protein_backbone(
    protein: Protein,
    camera: Camera,
    renderer: BrailleRenderer,
    color_scheme: ColorScheme = ColorScheme.SECONDARY_STRUCTURE,
):
    """Render protein backbone (CA trace) with depth fog and thick lines"""
    width = renderer.pixel_width
    height = renderer.pixel_height

    # Calculate max depth for fog
    all_z_coords = []

    for chain_idx, chain in enumerate(protein.chains):
        ca_atoms = [r.ca_atom for r in chain.residues if r.ca_atom and not r.is_ligand]

        if len(ca_atoms) < 2:
            continue

        # Project all CA atoms
        coords = np.array([[a.x, a.y, a.z] for a in ca_atoms])
        projected = camera.project(coords, width, height)

        # Collect Z coords for fog calculation
        all_z_coords.extend(projected[:, 2])

    # Calculate max depth for fog effect
    if all_z_coords:
        max_depth = max(abs(z) for z in all_z_coords)
    else:
        max_depth = 100.0

    # Now render with fog
    for chain_idx, chain in enumerate(protein.chains):
        ca_atoms = [r.ca_atom for r in chain.residues if r.ca_atom and not r.is_ligand]

        if len(ca_atoms) < 2:
            continue

        # Project all CA atoms
        coords = np.array([[a.x, a.y, a.z] for a in ca_atoms])
        projected = camera.project(coords, width, height)

        # Draw thick lines between consecutive CA atoms with depth fog
        for i in range(len(ca_atoms) - 1):
            residue = chain.residues[i]
            next_residue = chain.residues[i + 1]

            # Get base color
            base_color = get_color_for_residue(residue, color_scheme, chain_idx, protein.n_residues)

            x0, y0, z0 = projected[i]
            x1, y1, z1 = projected[i + 1]

            # Average Z for fog calculation
            z_avg = (z0 + z1) / 2

            # Apply depth fog
            fogged_color = apply_depth_fog(base_color, z_avg, max_depth)

            # Draw simple lines for better performance
            # All lines same thickness for fluid rotation
            renderer.draw_line(
                int(x0), int(y0), z0,
                int(x1), int(y1), z1,
                fogged_color
            )


def render_ligands_ball_stick(
    protein: Protein,
    camera: Camera,
    renderer: BrailleRenderer,
):
    """Render ligands in ball-and-stick style with enhanced visibility"""
    width = renderer.pixel_width
    height = renderer.pixel_height

    # Calculate max depth for fog
    all_z_coords = []

    for ligand_res, ligand_chain in protein.ligands:
        if ligand_res.is_water:
            continue

        atoms = ligand_res.atoms
        if not atoms:
            continue

        atom_coords = np.array([[a.x, a.y, a.z] for a in atoms])
        projected = camera.project(atom_coords, width, height)
        all_z_coords.extend(projected[:, 2])

    max_depth = max(abs(z) for z in all_z_coords) if all_z_coords else 100.0

    # Render ligands with depth fog
    for ligand_res, ligand_chain in protein.ligands:
        if ligand_res.is_water:
            continue  # Skip water

        atoms = ligand_res.atoms
        if not atoms:
            continue

        # Project atoms
        atom_coords = np.array([[a.x, a.y, a.z] for a in atoms])
        projected = camera.project(atom_coords, width, height)

        # Draw bonds first (so atoms are on top) - THICKER
        for i, atom1 in enumerate(atoms):
            for j, atom2 in enumerate(atoms[i+1:], start=i+1):
                # Simple distance check for bonding
                dist = atom1.distance_to(atom2)
                if dist < 2.0:  # Typical bond length cutoff
                    x0, y0, z0 = projected[i]
                    x1, y1, z1 = projected[j]

                    # Average Z for fog
                    z_avg = (z0 + z1) / 2

                    # Bond color with fog
                    bond_color = apply_depth_fog((140, 140, 140), z_avg, max_depth)

                    # Draw simple bond line (lighter for performance)
                    renderer.draw_line(
                        int(x0), int(y0), z0,
                        int(x1), int(y1), z1,
                        bond_color
                    )

        # Draw atoms as smaller circles for performance
        for i, atom in enumerate(atoms):
            x, y, z = projected[i]
            base_color = get_element_color(atom.element)

            # Apply depth fog to atom color
            fogged_color = apply_depth_fog(base_color, z, max_depth)

            # Smaller atom radius for better performance
            radius = int(2.5 * camera.zoom)  # Reduced for performance
            radius = max(1, min(radius, 6))  # Smaller range

            renderer.draw_circle(int(x), int(y), radius, z, fogged_color, filled=True)


