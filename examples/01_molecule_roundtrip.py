#!/usr/bin/env python3
"""Example 1 — small-molecule round-trip (offline, RDKit only).

SMILES -> MoleCode Mermaid graph -> SMILES, showing that the conversion is
deterministic and lossless. Run:

    python examples/01_molecule_roundtrip.py
"""

from rdkit import Chem

from molecode.molecule import mol_to_mermaid, mermaid_to_mol, mol_to_smiles

EXAMPLES = [
    ("Ethanol", "CCO"),
    ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("Caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("para-Chlorophenol", "Oc1ccc(Cl)cc1"),
]


def main() -> None:
    for name, smiles in EXAMPLES:
        mol = Chem.MolFromSmiles(smiles)
        mermaid = mol_to_mermaid(mol, name=name)

        recovered = mermaid_to_mol(mermaid)
        recovered_smiles = mol_to_smiles(recovered)
        match = Chem.CanonSmiles(recovered_smiles) == Chem.CanonSmiles(smiles)

        print(f"\n{'=' * 64}\n{name}: {smiles}\n{'=' * 64}")
        print(mermaid)
        print(f"-> back to SMILES: {recovered_smiles}   "
              f"[{'round-trip OK' if match else 'MISMATCH'}]")


if __name__ == "__main__":
    main()
