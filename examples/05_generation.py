#!/usr/bin/env python3
"""Example 5 — molecular GENERATION with an LLM.

Task family: de novo design under constraints. The model emits a MoleCode
Mermaid graph; because the output IS the structure, you can parse and validate
it deterministically with RDKit — no fragile SMILES post-processing.

    python examples/05_generation.py
"""

from molecode.molecule import mermaid_to_mol, mol_to_smiles, has_invalid_atoms
from molecode.prompts import MOLECULE_SYSTEM_PROMPT
from _llm import call_llm

USER = """Task: Design a molecule and DRAW it as a Mermaid graph.

Constraints:
- Exactly 8 carbon atoms
- Contains one aromatic ring
- Contains exactly one carboxylic acid group (-C(=O)OH)

Output ONLY a fenced Mermaid code block:
```mermaid
...graph...
```
Use the MoleCode atom/bond grammar (Kekulé form for aromatic rings)."""


def extract_mermaid(text: str) -> str:
    """Pull the fenced ```mermaid ... ``` block out of an LLM reply."""
    if "```" not in text:
        return text
    block = text.split("```", 2)[1]
    if block.startswith("mermaid"):
        block = block[len("mermaid"):]
    return block.strip()


def main() -> None:
    reply = call_llm(MOLECULE_SYSTEM_PROMPT, USER)
    if reply is None:
        print("\n(With an API key set, this validates the generated graph below.)")
        return

    mermaid = extract_mermaid(reply)
    mol = mermaid_to_mol(mermaid, strict=False)
    if mol is None or has_invalid_atoms(mol):
        print("Generated graph is INVALID.")
        return
    print("Generated molecule SMILES:", mol_to_smiles(mol))


if __name__ == "__main__":
    main()
