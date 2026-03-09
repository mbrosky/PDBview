"""Information panel widget"""

from textual.widgets import Static
from rich.text import Text
from ..model.protein import Protein


class InfoPanel(Static):
    """Display protein information"""

    def __init__(self, protein: Protein, **kwargs):
        super().__init__(**kwargs)
        self.protein = protein

    def render(self) -> Text:
        """Render protein info with metadata"""
        lines = []

        # PDB ID and title
        if self.protein.metadata and self.protein.metadata.pdb_id:
            lines.append(f"[bold cyan]{self.protein.metadata.pdb_id}[/bold cyan] - {self.protein.metadata.title[:40]}")
        else:
            lines.append(f"[bold cyan]{self.protein.name}[/bold cyan]")

        lines.append("")  # Separator

        # Basic structure info - use official PDB stats if available
        lines.append(f"[dim]Structure:[/dim]")

        # Use official PDB statistics when available
        if self.protein.metadata:
            meta = self.protein.metadata
            if meta.pdb_chain_count is not None:
                lines.append(f"  Chains: {meta.pdb_chain_count}")
            else:
                lines.append(f"  Chains: {len(self.protein.chains)}")

            if meta.pdb_residue_count is not None:
                lines.append(f"  Residues: {meta.pdb_residue_count}")
            else:
                lines.append(f"  Residues: {self.protein.n_residues}")

            if meta.pdb_atom_count is not None:
                lines.append(f"  Atoms: {meta.pdb_atom_count}")
            else:
                lines.append(f"  Atoms: {self.protein.n_atoms}")
        else:
            # Fallback to our own counts
            lines.append(f"  Chains: {len(self.protein.chains)}")
            lines.append(f"  Residues: {self.protein.n_residues}")
            lines.append(f"  Atoms: {self.protein.n_atoms}")

        lines.append(f"  Ligands: {self.protein.n_ligands}")

        # Metadata if available
        if self.protein.metadata:
            meta = self.protein.metadata

            if meta.classification:
                lines.append("")
                lines.append(f"[dim]Classification:[/dim]")
                lines.append(f"  {meta.classification}")

            if meta.organism:
                lines.append("")
                lines.append(f"[dim]Organism:[/dim]")
                org_display = meta.organism
                if meta.organism_common:
                    org_display += f" ({meta.organism_common})"
                lines.append(f"  {org_display}")

            if meta.method:
                lines.append("")
                lines.append(f"[dim]Method:[/dim]")
                method_str = f"  {meta.method}"
                if meta.resolution:
                    method_str += f"\n  Resolution: {meta.resolution:.2f} Å"
                lines.append(method_str)

            if meta.has_mutations:
                lines.append("")
                lines.append(f"[yellow]⚠ Contains mutations[/yellow]")

        return Text.from_markup("\n".join(lines))

