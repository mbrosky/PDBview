"""Main Textual application"""

from textual.app import App, ComposeResult
from textual.containers import Container, Horizontal, Vertical
from textual.widgets import Header, Footer, Static, Label
from textual.binding import Binding
from textual import events

from ..model.protein import Protein
from ..model.binding_site import detect_binding_sites
from ..render.camera import Camera
from ..render.braille import BrailleRenderer, render_protein_backbone, render_ligands_ball_stick
from ..render.colors import ColorScheme
from .viewport import Viewport
from .ligand_panel import LigandPanel
from .info_panel import InfoPanel


class PDBViewApp(App):
    """PDBView - Terminal Protein and Ligand Viewer"""

    CSS = """
    #viewport {
        width: 100%;
        height: 1fr;
    }

    #sidebar {
        width: 40;
        dock: right;
        background: $panel;
    }

    #info {
        height: auto;
        max-height: 20;
        background: $panel;
    }

    .hidden {
        display: none;
    }
    """

    BINDINGS = [
        Binding("q", "quit", "Quit"),
        Binding("h", "rotate_left", "/l Rotate"),
        Binding("l", "rotate_right", "", show=False),
        Binding("j", "rotate_down", "/k Rotate"),
        Binding("k", "rotate_up", "", show=False),
        Binding("u", "roll_left", "/i Roll"),
        Binding("i", "roll_right", "", show=False),
        Binding("+", "zoom_in", "/- Zoom"),
        Binding("-", "zoom_out", "", show=False),
        Binding("=", "zoom_in", "", show=False),
        Binding("_", "zoom_out", "", show=False),
        Binding("r", "reset_view", "Reset"),
        Binding("c", "cycle_color", "Color"),
        Binding("v", "toggle_ligands", "Ligands"),
        Binding("f", "toggle_sidebar", "Info"),
        Binding("space", "toggle_rotation", "Auto-rotate"),
        Binding("?", "show_help", "Help"),
    ]

    def __init__(self, protein: Protein):
        super().__init__()
        self.protein = protein
        self.camera = Camera(distance=max(50.0, protein.bounding_radius() * 1.5))
        self.color_scheme = ColorScheme.SECONDARY_STRUCTURE
        self.show_ligands = True
        self.show_sidebar = True
        self.auto_rotate = False

        # Detect binding sites
        self.binding_sites = detect_binding_sites(protein)

        # Center protein at origin
        self.protein.center_at_origin()

    def compose(self) -> ComposeResult:
        """Create UI layout"""
        yield Header()

        with Horizontal():
            yield Viewport(
                self.protein,
                self.camera,
                self.color_scheme,
                self.show_ligands,
                id="viewport"
            )

            with Vertical(id="sidebar"):
                yield InfoPanel(self.protein, id="info")
                yield LigandPanel(self.protein, self.binding_sites, id="ligands")

        yield Footer()

    def on_mount(self) -> None:
        """Called when app is mounted"""
        self.title = f"PDBView - {self.protein.name}"
        self.update_viewport()

        # Set up auto-rotation timer if enabled
        if self.auto_rotate:
            self.set_interval(1/30, self.auto_rotate_step)

    def update_viewport(self):
        """Refresh the viewport"""
        viewport = self.query_one("#viewport", Viewport)
        viewport.camera = self.camera
        viewport.color_scheme = self.color_scheme
        viewport.show_ligands = self.show_ligands
        viewport.refresh()

    def action_rotate_left(self):
        """Rotate left"""
        self.camera.rotate(0, -5, 0)
        self.update_viewport()

    def action_rotate_right(self):
        """Rotate right"""
        self.camera.rotate(0, 5, 0)
        self.update_viewport()

    def action_rotate_up(self):
        """Rotate up"""
        self.camera.rotate(-5, 0, 0)
        self.update_viewport()

    def action_rotate_down(self):
        """Rotate down"""
        self.camera.rotate(5, 0, 0)
        self.update_viewport()

    def action_roll_left(self):
        """Roll left"""
        self.camera.rotate(0, 0, -5)
        self.update_viewport()

    def action_roll_right(self):
        """Roll right"""
        self.camera.rotate(0, 0, 5)
        self.update_viewport()

    def action_zoom_in(self):
        """Zoom in"""
        self.camera.zoom_in(1.1)
        self.update_viewport()

    def action_zoom_out(self):
        """Zoom out"""
        self.camera.zoom_out(1.1)
        self.update_viewport()

    def action_reset_view(self):
        """Reset camera"""
        self.camera.reset()
        self.update_viewport()

    def action_cycle_color(self):
        """Cycle through color schemes"""
        schemes = list(ColorScheme)
        current_idx = schemes.index(self.color_scheme)
        self.color_scheme = schemes[(current_idx + 1) % len(schemes)]
        self.update_viewport()

    def action_toggle_ligands(self):
        """Toggle ligand visibility"""
        self.show_ligands = not self.show_ligands
        self.update_viewport()

    def action_toggle_sidebar(self):
        """Toggle sidebar visibility"""
        self.show_sidebar = not self.show_sidebar
        sidebar = self.query_one("#sidebar")
        sidebar.set_class(not self.show_sidebar, "hidden")

    def action_toggle_rotation(self):
        """Toggle auto-rotation"""
        self.auto_rotate = not self.auto_rotate
        if self.auto_rotate:
            self.set_interval(1/30, self.auto_rotate_step)

    def auto_rotate_step(self):
        """Auto-rotation animation step"""
        if self.auto_rotate:
            self.camera.rotate(0, 1, 0)
            self.update_viewport()

    def action_show_help(self):
        """Show help overlay"""
        pass
