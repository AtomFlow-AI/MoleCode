from rdkit import Chem

from molecode.molecule import mermaid_to_mol, mol_to_mermaid


def _canonical_isomeric_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def _round_trip_smiles(smiles: str, *, kekulize: bool = True) -> tuple[str, str]:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    graph = mol_to_mermaid(mol, name="Test", kekulize=kekulize)
    recovered = mermaid_to_mol(graph)
    assert recovered is not None
    return graph, Chem.MolToSmiles(recovered, isomericSmiles=True)


def test_aromatic_bond_round_trip_when_not_kekulized():
    graph, recovered = _round_trip_smiles("c1ccccc1", kekulize=False)

    assert "<-->" in graph
    assert recovered == _canonical_isomeric_smiles("c1ccccc1")


def test_double_bond_e_stereo_round_trip():
    graph, recovered = _round_trip_smiles("F/C=C/F")

    assert "===|E|" in graph
    assert recovered == _canonical_isomeric_smiles("F/C=C/F")


def test_double_bond_z_stereo_round_trip():
    graph, recovered = _round_trip_smiles("F/C=C\\F")

    assert "===|Z|" in graph
    assert recovered == _canonical_isomeric_smiles("F/C=C\\F")


def test_tetrahedral_r_stereo_round_trip_uses_absolute_cip_label():
    graph, recovered = _round_trip_smiles("C[C@H](O)Cl")

    assert "_R[CH]" in graph
    assert recovered == _canonical_isomeric_smiles("C[C@H](O)Cl")


def test_tetrahedral_s_stereo_round_trip_uses_absolute_cip_label():
    graph, recovered = _round_trip_smiles("C[C@@H](O)Cl")

    assert "_S[CH]" in graph
    assert recovered == _canonical_isomeric_smiles("C[C@@H](O)Cl")
