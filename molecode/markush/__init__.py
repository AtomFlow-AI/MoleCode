"""MoleCode — Markush / abbreviated-group representation.

Extends the small-molecule grammar with **abbreviation nodes** written with
curly braces (e.g. ``Mol_X_1{Boc}``, ``{R1}``, ``{Ar}``) so that variable
R-groups and named substituents can be represented without expanding them into
atoms. Includes an RDKit-free graph-isomorphism comparator (``egl_isomorphic``)
that scores predictions up to abbreviation expansion and Kekulé ambiguity.

    >>> from molecode.markush import mol_to_mermaid, mermaid_to_mol
    >>> from molecode.markush import EGLGraph, egl_isomorphic, EXPAND_MAP
"""

from .rdkit_to_mermaid import MolToMermaidConverter, mol_to_mermaid
from .mermaid_to_rdkit import (
    MermaidMolParser,
    mermaid_to_mol,
    mol_to_smiles,
    mol_to_inchi,
    has_invalid_atoms,
    get_invalid_atom_labels,
)
from .egl_graph import (
    EGLGraph,
    NodeInfo,
    EdgeInfo,
    normalize_abbrev_name,
    egl_isomorphic,
)
from .abbreviation_map import (
    SINGLE_ATOM_MAP,
    SUBGRAPH_MAP,
    NON_EXPANDABLE,
    build_expand_map,
    EXPAND_MAP,
)

__all__ = [
    "MolToMermaidConverter",
    "mol_to_mermaid",
    "MermaidMolParser",
    "mermaid_to_mol",
    "mol_to_smiles",
    "mol_to_inchi",
    "has_invalid_atoms",
    "get_invalid_atom_labels",
    "EGLGraph",
    "NodeInfo",
    "EdgeInfo",
    "normalize_abbrev_name",
    "egl_isomorphic",
    "SINGLE_ATOM_MAP",
    "SUBGRAPH_MAP",
    "NON_EXPANDABLE",
    "build_expand_map",
    "EXPAND_MAP",
]
