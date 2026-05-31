#!/usr/bin/env python3
"""Example 2 — polymer round-trip (offline, RDKit only).

A polymer is given as a repeat unit in PSMILES form (two ``*`` attachment
points). MoleCode keeps the repeat unit explicit and carries the repetition
count symbolically as ``×n`` (so the graph does not blow up with chain length),
then converts back to canonical PSMILES.

    python examples/02_polymer_roundtrip.py
"""

from rdkit import Chem

from molecode.polymer import (
    polymer_to_mermaid,
    mermaid_to_psmiles,
    block_copolymer_to_mermaid,
    BlockSpec,
)

HOMOPOLYMERS = [
    ("Polyethylene (PE)", "*CC*", 1000),
    ("Polypropylene (PP)", "*CC(C)*", 500),
    ("Poly(ethylene glycol) (PEG)", "*CCO*", 300),
    ("Nylon-6", "*NCCCCCC(=O)*", 8),
    ("Polystyrene (PS)", "*CC(c1ccccc1)*", 200),
]


def main() -> None:
    for name, psmiles, n in HOMOPOLYMERS:
        mermaid = polymer_to_mermaid(psmiles, n=n, name=name)
        recovered = mermaid_to_psmiles(mermaid)
        match = recovered is not None and (
            Chem.CanonSmiles(recovered) == Chem.CanonSmiles(psmiles)
        )
        print(f"\n{'=' * 64}\n{name}  (repeat ×{n})\n{'=' * 64}")
        print(mermaid)
        print(f"-> back to PSMILES: {recovered}   "
              f"[{'round-trip OK' if match else 'MISMATCH'}]")

    # Block copolymer: PEG-b-PPO
    print(f"\n{'=' * 64}\nBlock copolymer: PEG-b-PPO\n{'=' * 64}")
    mermaid = block_copolymer_to_mermaid(
        [BlockSpec("*CCO*", n=20, label="PEG"),
         BlockSpec("*CC(C)O*", n=10, label="PPO")],
        name="PEG-b-PPO",
    )
    print(mermaid)


if __name__ == "__main__":
    main()
