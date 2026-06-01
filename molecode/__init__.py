"""MoleCode — an LLM-native, graph-explicit molecular language.

MoleCode serializes molecules as Mermaid graphs in which every atom and bond is
a typed declaration with a persistent identifier, so molecular topology is
directly readable, editable, and auditable inside an LLM context window. It is
deterministically and losslessly inter-convertible with SMILES / MOL via RDKit
(no learned model). The same Subgraph–Node–Edge grammar spans three domains:

    molecode.molecule  — small molecules
    molecode.polymer   — polymers (explicit repeat unit + ``×n``)
    molecode.markush   — Markush structures (``{}`` abbreviation nodes)
    molecode.prompts   — LLM system prompts for each grammar

Quick start::

    from molecode import mol_to_mermaid, mermaid_to_mol, mol_to_smiles
    from rdkit import Chem

    mermaid = mol_to_mermaid(Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O"),
                             name="Aspirin")
    print(mermaid)
    assert mol_to_smiles(mermaid_to_mol(mermaid))  # round-trips
"""

from . import molecule, polymer, markush, prompts

# Small molecules are the default top-level converters.
from .molecule import (
    mol_to_mermaid,
    mermaid_to_mol,
    mol_to_smiles,
    mol_to_inchi,
)

# Optional, dependency-free OpenAI-compatible client for driving LLM tasks.
from .llm import LLMClient

__version__ = "0.1.0"

__all__ = [
    "molecule",
    "polymer",
    "markush",
    "prompts",
    "LLMClient",
    "mol_to_mermaid",
    "mermaid_to_mol",
    "mol_to_smiles",
    "mol_to_inchi",
    "__version__",
]
