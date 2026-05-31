#!/usr/bin/env python3
"""Example 4 — molecular UNDERSTANDING with an LLM.

Task family: read a molecule and answer a structural question (here: count the
carbon atoms; the same pattern covers molecular-formula prediction, ring
counting, functional-group identification, IUPAC naming, etc.).

The point of MoleCode: because every atom is an explicit, typed node, the model
reads structure directly instead of reconstructing it from SMILES syntax.

Runs offline as a "dry run" (prints the exact prompt). Set MOLECODE_API_KEY to
actually call an LLM — see examples/_llm.py.

    python examples/04_understanding.py
"""

from rdkit import Chem

from molecode.molecule import mol_to_mermaid
from molecode.prompts import MOLECULE_SYSTEM_PROMPT
from _llm import call_llm

SMILES = "CC(=O)Oc1ccccc1C(=O)O"   # Aspirin — 9 carbons

USER_TEMPLATE = """Task: Count the number of carbon atoms in the given molecule \
represented as a Mermaid graph.

You must analyze ONLY the Mermaid graph below. Do not use any SMILES information.

Molecule (Mermaid):
```mermaid
{mermaid_graph}
```

Instructions:
- Analyze the Mermaid graph structure carefully
- Count all carbon atoms (nodes labeled with 'C')
- Give your answer as an integer only

Answer:"""


def main() -> None:
    mermaid = mol_to_mermaid(Chem.MolFromSmiles(SMILES), name="Aspirin")
    user = USER_TEMPLATE.format(mermaid_graph=mermaid)

    reply = call_llm(MOLECULE_SYSTEM_PROMPT, user)
    if reply is not None:
        print("\nLLM answer:", reply.strip(), "  (ground truth: 9)")


if __name__ == "__main__":
    main()
