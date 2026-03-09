"""3D viewport widget"""

from textual.widget import Widget
from textual.reactive import reactive
from rich.text import Text
from typing import List

from ..model.protein import Protein
from ..render.camera import Camera
from ..render.braille import BrailleRenderer, render_protein_backbone, render_ligands_ball_stick
from ..render.colors import ColorScheme


class Viewport(Widget):
    """3D structure viewport"""

    protein: Protein
    camera: Camera
    color_scheme: ColorScheme = reactive(ColorScheme.SECONDARY_STRUCTURE)
    show_ligands: bool = reactive(True)

    def __init__(
        self,
        protein: Protein,
        camera: Camera,
        color_scheme: ColorScheme,
        show_ligands: bool,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.protein = protein
        self.camera = camera
        self.color_scheme = color_scheme
        self.show_ligands = show_ligands

    def watch_color_scheme(self, new_scheme: ColorScheme) -> None:
        """React to color scheme changes"""
        self.refresh()

    def watch_show_ligands(self, show: bool) -> None:
        """React to ligand visibility changes"""
        self.refresh()

    def render(self) -> Text:
        """Render the 3D structure with Braille"""
        # Get viewport size
        width = self.size.width
        height = self.size.height

        if width == 0 or height == 0:
            return Text("")

        # Create Braille renderer
        renderer = BrailleRenderer(width=width, height=height)
        renderer.clear()

        # Render protein backbone with depth fog and thick lines
        render_protein_backbone(self.protein, self.camera, renderer, self.color_scheme)

        # Render ligands if enabled
        if self.show_ligands:
            render_ligands_ball_stick(self.protein, self.camera, renderer)

        # Convert to Rich Text
        markup_str = renderer.render_to_string()
        return Text.from_markup(markup_str)
