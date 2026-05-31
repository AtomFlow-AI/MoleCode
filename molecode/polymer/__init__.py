"""MoleCode — polymer representation.

Converts polymer repeat units (PSMILES with two ``*`` attachment points) and
block copolymers into a MoleCode Mermaid graph that keeps the repeat unit
explicit and carries the repetition count symbolically as ``×n``, and parses
that graph back to a canonical PSMILES.

    >>> from molecode.polymer import polymer_to_mermaid, mermaid_to_psmiles
    >>> mermaid = polymer_to_mermaid("*CC*", n=1000, name="Polyethylene")
    >>> mermaid_to_psmiles(mermaid)
    '*CC*'
"""

from .polymer_to_mermaid import (
    BlockSpec,
    PolymerSpec,
    RepeatUnitConverter,
    PolymerToMermaidConverter,
    polymer_to_mermaid,
    block_copolymer_to_mermaid,
)
from .mermaid_to_psmiles import mermaid_to_psmiles, parse_element

__all__ = [
    "BlockSpec",
    "PolymerSpec",
    "RepeatUnitConverter",
    "PolymerToMermaidConverter",
    "polymer_to_mermaid",
    "block_copolymer_to_mermaid",
    "mermaid_to_psmiles",
    "parse_element",
]
