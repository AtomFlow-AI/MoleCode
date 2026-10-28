"""MoleCode — small-molecule representation.

Deterministic, RDKit-backed, lossless conversion between RDKit molecules /
SMILES and the MoleCode Mermaid graph language for small molecules.

    >>> from molecode.molecule import mol_to_mermaid, mermaid_to_mol, mol_to_smiles
    >>> from rdkit import Chem
    >>> mermaid = mol_to_mermaid(Chem.MolFromSmiles("CCO"), name="Ethanol")
    >>> mol_to_smiles(mermaid_to_mol(mermaid))
    'CCO'
"""

from .rdkit_to_mermaid import MolToMermaidConverter, mol_to_mermaid
from . import stereo
from .stereo import (
    CIP_IMPL,
    PROFILE_ID,
    StereoGroupRecord,
    StereoInfo,
    AxialRecord,
    ParsedStereoLayer,
    apply_stereo_layer,
    extract_stereo_info,
    parse_stereo_layer,
    prepare_mol,
    render_stereo_layer,
    stereo_coverage,
    stereo_fingerprint,
    validate_document,
)
from .mermaid_to_rdkit import (
    MermaidMolParser,
    mermaid_to_mol,
    mol_to_smiles,
    mol_to_inchi,
    has_invalid_atoms,
    get_invalid_atom_labels,
    visualize_mol,
    visualize_mols_grid,
)

__all__ = [
    "MolToMermaidConverter",
    "mol_to_mermaid",
    "stereo",
    "CIP_IMPL",
    "PROFILE_ID",
    "StereoInfo",
    "StereoGroupRecord",
    "AxialRecord",
    "ParsedStereoLayer",
    "prepare_mol",
    "extract_stereo_info",
    "render_stereo_layer",
    "parse_stereo_layer",
    "apply_stereo_layer",
    "stereo_fingerprint",
    "stereo_coverage",
    "validate_document",
    "MermaidMolParser",
    "mermaid_to_mol",
    "mol_to_smiles",
    "mol_to_inchi",
    "has_invalid_atoms",
    "get_invalid_atom_labels",
    "visualize_mol",
    "visualize_mols_grid",
]
