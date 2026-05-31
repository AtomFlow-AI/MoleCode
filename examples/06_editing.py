#!/usr/bin/env python3
"""Example 6 — molecular EDITING with an LLM.

Task family: structure-aware editing (add / delete / substitute atoms). In
MoleCode an edit is a *local graph operation* — e.g. adding a methyl group is
literally one new node plus one new edge — so edits are localized, auditable,
and diff-able instead of a whole-string rewrite.

This example shows the edit two ways:
  1. mechanically, by editing the Mermaid text directly (offline, runnable);
  2. via an LLM, by sending the graph + a natural-language instruction.

    python examples/06_editing.py
"""

from rdkit import Chem

from molecode.molecule import mol_to_mermaid, mermaid_to_mol, mol_to_smiles
from molecode.prompts import MOLECULE_SYSTEM_PROMPT
from _llm import call_llm

SMILES = "CCO"   # ethanol -> we will add a methyl to make propan-1-ol


def mechanical_edit() -> None:
    """Add one carbon node + one bond directly in the Mermaid text."""
    mermaid = mol_to_mermaid(Chem.MolFromSmiles(SMILES), name="Ethanol")
    print("Original graph:\n" + mermaid)

    # Graft a new carbon onto the terminal CH3 (Ethanol_C_1). Because MoleCode
    # encodes hydrogen counts explicitly, a correct local edit is two changes:
    #   1. the new methyl node + the new bond, and
    #   2. relabel Ethanol_C_1 from [CH3] to [CH2] (it lost one H to the bond).
    lines = mermaid.splitlines()
    out = []
    for line in lines:
        if line.strip() == "Ethanol_C_1[CH3]":
            line = line.replace("[CH3]", "[CH2]")  # H-count bookkeeping
        if line.strip() == "end":
            out.append("        Ethanol_C_3[CH3]            %% new methyl node")
            out.append("        Ethanol_C_1 --- Ethanol_C_3 %% new bond (the edit)")
        out.append(line)
    edited = "\n".join(out)
    print("\nEdited graph (C_1: CH3->CH2, added Ethanol_C_3 + one bond):\n" + edited)

    new_smiles = mol_to_smiles(mermaid_to_mol(edited))
    print("\nResult SMILES:", new_smiles, "(propan-1-ol = CCCO)")


def llm_edit() -> None:
    mermaid = mol_to_mermaid(Chem.MolFromSmiles(SMILES), name="Ethanol")
    user = f"""Here is a molecule as a MoleCode Mermaid graph:
```mermaid
{mermaid}
```
Edit instruction: add a methyl group (-CH3) to the terminal CH3 carbon, turning
ethanol into propan-1-ol. Output ONLY the edited Mermaid graph in a fenced block."""
    call_llm(MOLECULE_SYSTEM_PROMPT, user)


def main() -> None:
    print("### Mechanical edit (offline) ###\n")
    mechanical_edit()
    print("\n\n### LLM edit ###\n")
    llm_edit()


if __name__ == "__main__":
    main()
