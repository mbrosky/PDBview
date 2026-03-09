"""PDB metadata parser for extracting header information"""

from dataclasses import dataclass, field
from typing import List, Optional
from pathlib import Path
import requests
import json


@dataclass
class PDBMetadata:
    """Metadata extracted from PDB file header"""
    pdb_id: str = ""
    title: str = ""
    classification: str = ""
    organism: str = ""
    organism_common: str = ""
    expression_system: str = ""
    method: str = ""
    resolution: Optional[float] = None
    keywords: str = ""
    ec_number: str = ""
    has_mutations: bool = False
    deposition_date: str = ""
    function_description: str = ""  # Functional description from API
    # Official PDB statistics (from RCSB)
    pdb_atom_count: Optional[int] = None
    pdb_residue_count: Optional[int] = None
    pdb_chain_count: Optional[int] = None

    def summary(self) -> str:
        """Get formatted summary of metadata"""
        lines = []
        if self.title:
            lines.append(self.title)
        if self.classification:
            lines.append(f"Classification: {self.classification}")
        if self.organism:
            lines.append(f"Organism: {self.organism}")
        if self.method:
            method_str = self.method
            if self.resolution:
                method_str += f" ({self.resolution:.2f} Å)"
            lines.append(f"Method: {method_str}")
        return "\n".join(lines)


def parse_pdb_metadata(filepath: str) -> PDBMetadata:
    """
    Parse metadata from PDB file header records.

    Args:
        filepath: Path to PDB file

    Returns:
        PDBMetadata object with extracted information
    """
    path = Path(filepath)
    if not path.exists():
        return PDBMetadata()

    metadata = PDBMetadata()

    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()

        # Parse header records
        title_lines = []
        compnd_lines = []
        source_lines = []

        for line in lines:
            # Stop at ATOM records (header is before them)
            if line.startswith('ATOM') or line.startswith('HETATM'):
                break

            # HEADER record
            if line.startswith('HEADER'):
                # Format: HEADER    CLASSIFICATION            DD-MMM-YY   IDCODE
                metadata.classification = line[10:50].strip()
                metadata.deposition_date = line[50:59].strip()
                metadata.pdb_id = line[62:66].strip()

            # TITLE records
            elif line.startswith('TITLE'):
                title_text = line[10:].strip()
                # Remove continuation number (e.g., "TITLE   2")
                if title_text and not title_text[0].isdigit():
                    title_lines.append(title_text)
                else:
                    # Parse "TITLE   2 text" format
                    parts = title_text.split(None, 1)
                    if len(parts) > 1:
                        title_lines.append(parts[1])

            # COMPND records
            elif line.startswith('COMPND'):
                compnd_lines.append(line[10:].strip())

            # SOURCE records
            elif line.startswith('SOURCE'):
                source_lines.append(line[10:].strip())

            # KEYWDS records
            elif line.startswith('KEYWDS'):
                if metadata.keywords:
                    metadata.keywords += " " + line[10:].strip()
                else:
                    metadata.keywords = line[10:].strip()

            # EXPDTA record
            elif line.startswith('EXPDTA'):
                metadata.method = line[10:].strip()

            # REMARK 2 for resolution
            elif line.startswith('REMARK   2 RESOLUTION'):
                # Format: REMARK   2 RESOLUTION.    1.80 ANGSTROMS.
                try:
                    res_str = line.split()[3]
                    metadata.resolution = float(res_str)
                except (ValueError, IndexError):
                    pass

        # Join title lines
        if title_lines:
            metadata.title = " ".join(title_lines)

        # Parse COMPND records
        compnd_text = " ".join(compnd_lines)
        if "MOLECULE:" in compnd_text:
            mol_start = compnd_text.find("MOLECULE:") + 9
            mol_end = compnd_text.find(";", mol_start)
            if mol_end > mol_start:
                metadata.title = compnd_text[mol_start:mol_end].strip()

        if "EC:" in compnd_text:
            ec_start = compnd_text.find("EC:") + 3
            ec_end = compnd_text.find(";", ec_start)
            if ec_end > ec_start:
                metadata.ec_number = compnd_text[ec_start:ec_end].strip()

        if "MUTATION: YES" in compnd_text:
            metadata.has_mutations = True

        # Parse SOURCE records
        source_text = " ".join(source_lines)
        if "ORGANISM_SCIENTIFIC:" in source_text:
            org_start = source_text.find("ORGANISM_SCIENTIFIC:") + 20
            org_end = source_text.find(";", org_start)
            if org_end > org_start:
                metadata.organism = source_text[org_start:org_end].strip()

        if "ORGANISM_COMMON:" in source_text:
            common_start = source_text.find("ORGANISM_COMMON:") + 16
            common_end = source_text.find(";", common_start)
            if common_end > common_start:
                metadata.organism_common = source_text[common_start:common_end].strip()

        if "EXPRESSION_SYSTEM:" in source_text:
            exp_start = source_text.find("EXPRESSION_SYSTEM:") + 18
            exp_end = source_text.find(";", exp_start)
            if exp_end > exp_start:
                metadata.expression_system = source_text[exp_start:exp_end].strip()

    except Exception as e:
        # Return empty metadata if parsing fails
        pass

    # Fetch additional metadata from RCSB API if we have a PDB ID
    if metadata.pdb_id:
        api_metadata = fetch_rcsb_metadata(metadata.pdb_id)
        if api_metadata:
            # Enrich with API data - function description
            if api_metadata.get('function_description') and not metadata.function_description:
                metadata.function_description = api_metadata['function_description']

            # Enrich with official PDB statistics
            if api_metadata.get('pdb_atom_count') is not None:
                metadata.pdb_atom_count = api_metadata['pdb_atom_count']
            if api_metadata.get('pdb_residue_count') is not None:
                metadata.pdb_residue_count = api_metadata['pdb_residue_count']
            if api_metadata.get('pdb_chain_count') is not None:
                metadata.pdb_chain_count = api_metadata['pdb_chain_count']

    return metadata


def fetch_rcsb_metadata(pdb_id: str) -> Optional[dict]:
    """
    Fetch additional metadata from RCSB PDB REST API, UniProt, and PubMed.

    Args:
        pdb_id: 4-letter PDB identifier

    Returns:
        Dictionary with additional metadata, or None if fetch fails
    """
    result = {}

    # Fetch official PDB statistics from RCSB
    try:
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.upper()}"
        response = requests.get(url, timeout=5)

        if response.status_code == 200:
            data = response.json()

            # Get official atom/residue/chain counts
            if 'rcsb_entry_info' in data:
                info = data['rcsb_entry_info']
                # Deposited atom count (non-hydrogen)
                if 'deposited_atom_count' in info:
                    result['pdb_atom_count'] = info['deposited_atom_count']
                # Deposited residue count
                if 'deposited_modeled_polymer_monomer_count' in info:
                    result['pdb_residue_count'] = info['deposited_modeled_polymer_monomer_count']
                # Chain count
                if 'polymer_entity_count_protein' in info:
                    result['pdb_chain_count'] = info['polymer_entity_count_protein']
    except Exception:
        pass

    # Try UniProt first (most detailed biological function)
    uniprot_desc = fetch_uniprot_description(pdb_id)
    if uniprot_desc:
        result['function_description'] = uniprot_desc
    else:
        # Fallback to PubMed abstract
        pubmed_desc = fetch_pubmed_description(pdb_id)
        if pubmed_desc:
            result['function_description'] = pubmed_desc

    return result if result else None


def fetch_uniprot_description(pdb_id: str) -> Optional[str]:
    """
    Fetch functional description from UniProt via RCSB mapping.

    Args:
        pdb_id: 4-letter PDB identifier

    Returns:
        Functional description string or None
    """
    try:
        # Get UniProt accession from RCSB (try entity 1 first)
        url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id.upper()}/1"
        response = requests.get(url, timeout=5)

        if response.status_code != 200:
            return None

        data = response.json()

        # Extract UniProt accession
        uniprot_id = None
        if 'rcsb_polymer_entity_container_identifiers' in data:
            identifiers = data['rcsb_polymer_entity_container_identifiers']
            # Try direct uniprot_ids first
            if 'uniprot_ids' in identifiers and identifiers['uniprot_ids']:
                uniprot_id = identifiers['uniprot_ids'][0]
            # Fallback to reference_sequence_identifiers
            elif 'reference_sequence_identifiers' in identifiers:
                for ref in identifiers['reference_sequence_identifiers']:
                    if ref.get('database_name') == 'UniProt':
                        uniprot_id = ref.get('database_accession')
                        break

        # Fetch function from UniProt
        if uniprot_id:
            uniprot_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
            uniprot_response = requests.get(uniprot_url, timeout=5)

            if uniprot_response.status_code == 200:
                uniprot_data = uniprot_response.json()

                # Get function comment (primary source of biological function)
                if 'comments' in uniprot_data:
                    for comment in uniprot_data['comments']:
                        if comment.get('commentType') == 'FUNCTION':
                            texts = comment.get('texts', [])
                            if texts and 'value' in texts[0]:
                                func_text = texts[0]['value']
                                # Return full description, cleaned up
                                # Remove PubMed references like (PubMed:12345678)
                                import re
                                func_text = re.sub(r'\(PubMed:\d+\)', '', func_text)
                                func_text = re.sub(r'\(Ref\.\d+\)', '', func_text)
                                # Clean up extra spaces
                                func_text = ' '.join(func_text.split())
                                return func_text

                # Fallback to protein name/description
                if 'proteinDescription' in uniprot_data:
                    rec_name = uniprot_data['proteinDescription'].get('recommendedName', {})
                    if rec_name:
                        full_name = rec_name.get('fullName', {}).get('value', '')
                        if full_name and full_name.lower() not in ['uncharacterized protein', 'unknown']:
                            return f"Functions as {full_name.lower()}"

    except Exception as e:
        # Silently fail if API unavailable
        pass

    return None


def fetch_pubmed_description(pdb_id: str) -> Optional[str]:
    """
    Fetch description from PubMed primary citation.

    Args:
        pdb_id: 4-letter PDB identifier

    Returns:
        Description from PubMed abstract or None
    """
    try:
        # Get primary citation PMID from RCSB
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.upper()}"
        response = requests.get(url, timeout=5)

        if response.status_code != 200:
            return None

        data = response.json()

        # Extract PMID from primary citation
        pmid = None
        if 'rcsb_primary_citation' in data:
            pmid = data['rcsb_primary_citation'].get('pdbx_database_id_PubMed')

        if not pmid:
            return None

        # Fetch abstract from PubMed
        pubmed_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            'db': 'pubmed',
            'id': pmid,
            'retmode': 'xml',
            'rettype': 'abstract'
        }

        pubmed_response = requests.get(pubmed_url, params=params, timeout=5)

        if pubmed_response.status_code == 200:
            # Parse XML to extract abstract
            import xml.etree.ElementTree as ET
            root = ET.fromstring(pubmed_response.content)

            # Find AbstractText
            abstract_texts = root.findall('.//AbstractText')
            if abstract_texts:
                # Combine all abstract sections
                abstract_parts = []
                for abstract in abstract_texts:
                    if abstract.text:
                        # Skip if it's just a label
                        label = abstract.get('Label', '')
                        if label and label.upper() in ['BACKGROUND', 'OBJECTIVE', 'METHODS', 'RESULTS', 'CONCLUSIONS']:
                            # For structured abstracts, only use BACKGROUND and CONCLUSIONS
                            if label.upper() in ['BACKGROUND', 'CONCLUSIONS']:
                                abstract_parts.append(abstract.text.strip())
                        else:
                            # For unstructured abstracts, take everything
                            abstract_parts.append(abstract.text.strip())

                if abstract_parts:
                    full_abstract = ' '.join(abstract_parts)
                    # Limit to reasonable length (first 3-4 sentences)
                    sentences = full_abstract.split('. ')
                    if len(sentences) > 4:
                        return '. '.join(sentences[:4]) + '.'
                    return full_abstract

    except Exception:
        pass

    return None
