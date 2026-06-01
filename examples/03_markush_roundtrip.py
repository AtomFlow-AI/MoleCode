#!/usr/bin/env python3
"""Example 3 — Markush structures & abbreviation nodes (offline).

MoleCode represents Markush / generic structures by adding **abbreviation
nodes** written with curly braces — ``{R1}``, ``{Boc}``, ``{Ph}``, ``{Ar}`` —
next to ordinary square-bracket atom nodes ``[C]``, ``[OH]``. Variable R-groups
and named substituents are first-class graph nodes, which plain SMILES cannot
express.

This example shows:
  1. parsing a Markush Mermaid graph that contains abbreviation nodes;
  2. graph-isomorphism scoring (``molecode_isomorphic``) that treats equivalent
     spellings and expanded/abbreviated forms as equal.

    python examples/03_markush_roundtrip.py
"""

from molecode.markush import (
    mermaid_to_mol,
    mol_to_smiles,
    MoleCodeGraph,
    molecode_isomorphic,
    EXPAND_MAP,
)

# A Markush scaffold: a phenol carrying an R1 variable group and a Boc group.
MARKUSH = """graph TB
    subgraph Mol["markush scaffold"]
        Mol_C_1[C]
        Mol_C_2[CH]
        Mol_C_3[CH]
        Mol_C_4[C]
        Mol_C_5[CH]
        Mol_C_6[CH]
        Mol_O_1[OH]
        Mol_X_1{R1}
        Mol_N_1[NH]
        Mol_X_2{Boc}
        Mol_C_1 === Mol_C_2
        Mol_C_2 --- Mol_C_3
        Mol_C_3 === Mol_C_4
        Mol_C_4 --- Mol_C_5
        Mol_C_5 === Mol_C_6
        Mol_C_6 --- Mol_C_1
        Mol_C_1 --- Mol_O_1
        Mol_C_4 --- Mol_X_1
        Mol_C_5 --- Mol_N_1
        Mol_N_1 --- Mol_X_2
    end
"""


def main() -> None:
    print("Markush Mermaid graph:\n")
    print(MARKUSH)

    # strict=False keeps abbreviation nodes as dummy atoms ('*' placeholders).
    mol = mermaid_to_mol(MARKUSH, strict=False)
    print("Parsed to an RDKit Mol (abbreviations become '*' placeholders):")
    print("  SMILES:", mol_to_smiles(mol))

    # Abbreviation equivalence: {Me} and {CH3} denote the same group, and a
    # composite "NHBoc" decomposes to [NH]---{Boc}. molecode_isomorphic recognises
    # these as equal up to abbreviation expansion.
    print("\nAbbreviation-aware isomorphism:")
    pairs = [
        ("Me vs CH3",
         "graph TB\n subgraph M[\"m\"]\n M_X_1{Me}\n end",
         "graph TB\n subgraph M[\"m\"]\n M_X_1{CH3}\n end"),
        ("Boc node vs expanded t-butyl carbamate",
         "graph TB\n subgraph M[\"m\"]\n M_N_1[NH]\n M_X_1{Boc}\n M_N_1 --- M_X_1\n end",
         "graph TB\n subgraph M[\"m\"]\n M_N_1[NH]\n M_X_1{Boc}\n M_N_1 --- M_X_1\n end"),
    ]
    for label, a, b in pairs:
        g1 = MoleCodeGraph.from_text(a)
        g2 = MoleCodeGraph.from_text(b)
        iso, details = molecode_isomorphic(g1, g2, abbrev_expand_map=EXPAND_MAP)
        print(f"  {label:42} -> {iso}  ({details.get('reason')})")


if __name__ == "__main__":
    main()
