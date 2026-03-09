"""Ligand information panel"""

from textual.widgets import Static
from textual.containers import ScrollableContainer
from typing import List

from ..model.protein import Protein
from ..model.binding_site import BindingSite


class LigandPanel(ScrollableContainer):
    """Display ligand and binding site information"""

    def __init__(self, protein: Protein, binding_sites: List[BindingSite], **kwargs):
        super().__init__(**kwargs)
        self.protein = protein
        self.binding_sites = binding_sites

    def compose(self):
        """Create panel content"""
        yield Static(self._render_content(), id="ligand_content")

    def _render_content(self) -> str:
        """Render ligand information"""
        lines = ["[bold yellow]═══ Ligands & Binding Sites ═══[/bold yellow]", ""]

        ligands = self.protein.ligands
        if not ligands:
            lines.append("[dim]No ligands detected[/dim]")
            return "\n".join(lines)

        lines.append(f"[bold]Total ligands: {len(ligands)}[/bold]")
        lines.append("")

        # Group ligands by type
        ligand_types = {}
        for lig_res, lig_chain in ligands:
            if lig_res.is_water:
                continue
            key = lig_res.name
            if key not in ligand_types:
                ligand_types[key] = []
            ligand_types[key].append((lig_res, lig_chain))

        # Display each ligand type
        for idx, (lig_name, instances) in enumerate(ligand_types.items(), 1):
            lines.append(f"[bold cyan]{idx}. {lig_name}[/bold cyan] ({len(instances)} instance{'s' if len(instances) > 1 else ''})")

            # Find binding site for first instance
            first_res, first_chain = instances[0]
            binding_site = None
            for site in self.binding_sites:
                if site.ligand_name == lig_name and site.ligand_chain == first_chain.id:
                    binding_site = site
                    break

            if binding_site and binding_site.n_contacts > 0:
                lines.append(f"   [green]✓[/green] {binding_site.n_contacts} residue contacts (< 5Å)")

                # Show top 5 closest residues
                top_residues = binding_site.get_top_residues(5)
                for res in top_residues:
                    lines.append(f"     • {res.res_name}{res.res_num} ({res.chain_id}) - {res.min_distance:.2f}Å")
            else:
                lines.append("   [dim]No binding site contacts[/dim]")

            lines.append("")

        # Summary statistics
        total_contacts = sum(site.n_contacts for site in self.binding_sites)
        if total_contacts > 0:
            lines.append("[bold]Summary:[/bold]")
            lines.append(f"Total binding residues: {total_contacts}")

        return "\n".join(lines)
