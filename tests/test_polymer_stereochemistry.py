"""Polymer round-trips preserve R/S chirality and E/Z double-bond stereo.

Mirrors tests/test_molecule_stereochemistry.py but for the polymer converters,
which encode a repeat unit (PSMILES with two ``*`` attachment points) as a
Mermaid graph and back. Before the CIP fix the polymer path dropped both
chirality and E/Z; these tests pin the corrected behaviour.
"""

from rdkit import Chem

from molecode.polymer import polymer_to_mermaid, mermaid_to_psmiles


def _round_trip(psmiles: str, n: int = 10) -> str:
    graph = polymer_to_mermaid(psmiles, n=n, name="Test")
    back = mermaid_to_psmiles(graph)
    assert back is not None
    return Chem.CanonSmiles(back)


def test_r_chirality_round_trip_uses_absolute_cip():
    ps = "*OC(=O)[C@@H](C)*"          # PLA-like, one stereocentre
    graph = polymer_to_mermaid(ps, n=10, name="Test")
    assert "_R" in graph or "_S" in graph
    assert _round_trip(ps) == Chem.CanonSmiles(ps)


def test_s_chirality_round_trip_uses_absolute_cip():
    ps = "*OC(=O)[C@H](C)*"
    assert _round_trip(ps) == Chem.CanonSmiles(ps)


def test_trans_double_bond_round_trip():
    ps = "*C/C=C/C*"
    graph = polymer_to_mermaid(ps, n=5, name="Test")
    assert "===|E|" in graph
    assert _round_trip(ps, n=5) == Chem.CanonSmiles(ps)


def test_cis_double_bond_round_trip():
    ps = "*C/C=C\\C*"
    graph = polymer_to_mermaid(ps, n=5, name="Test")
    assert "===|Z|" in graph
    assert _round_trip(ps, n=5) == Chem.CanonSmiles(ps)


def test_cis_and_trans_are_distinguished():
    assert _round_trip("*C/C=C/C*", n=5) != _round_trip("*C/C=C\\C*", n=5)


def test_chirality_and_ez_combined():
    ps = "*O[C@@H](C)/C=C/C*"
    assert _round_trip(ps) == Chem.CanonSmiles(ps)


def test_plain_polymers_still_round_trip():
    for ps in ("*CC*", "*CC(C)*", "*CCO*", "*NCCCCCC(=O)*", "*CC(c1ccccc1)*"):
        assert _round_trip(ps) == Chem.CanonSmiles(ps)
