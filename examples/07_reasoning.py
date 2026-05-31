#!/usr/bin/env python3
"""Example 7 — molecular REASONING with an LLM (reaction prediction).

Task family: multi-step chemical reasoning. The model receives the reactants as
MoleCode graphs and must draw the product as a MoleCode graph. Because atom
identities are persistent nodes, the model can track atoms across the
transformation and you can validate the predicted product deterministically.

    python examples/07_reasoning.py
"""

from rdkit import Chem

from molecode.molecule import mol_to_mermaid, mermaid_to_mol, mol_to_smiles
from molecode.prompts import MOLECULE_SYSTEM_PROMPT
from _llm import call_llm

# Fischer esterification: acetic acid + ethanol -> ethyl acetate (+ water)
REACTANTS = [("AceticAcid", "CC(=O)O"), ("Ethanol", "CCO")]


def main() -> None:
    blocks = []
    for name, smi in REACTANTS:
        blocks.append(f"{name}:\n```mermaid\n"
                      f"{mol_to_mermaid(Chem.MolFromSmiles(smi), name=name)}\n```")
    reactants_text = "\n\n".join(blocks)

    user = f"""Task: Predict the main organic product of the reaction and draw it \
as a MoleCode Mermaid graph.

Reactant Molecules:
{reactants_text}

Instructions:
- Reaction: acid-catalysed Fischer esterification.
- Give ONLY the main organic product (ignore eliminated water).
- You MAY include %% comments inside the graph to show your reasoning.
- Output ONLY a fenced ```mermaid ... ``` block."""

    reply = call_llm(MOLECULE_SYSTEM_PROMPT, user)
    if reply is None:
        print("\nExpected product: ethyl acetate (CCOC(C)=O)")
        return

    block = reply.split("```", 2)[1] if "```" in reply else reply
    if block.startswith("mermaid"):
        block = block[len("mermaid"):]
    mol = mermaid_to_mol(block.strip(), strict=False)
    if mol is not None:
        print("Predicted product SMILES:", mol_to_smiles(mol),
              "  (expected ethyl acetate: CCOC(C)=O)")


if __name__ == "__main__":
    main()
