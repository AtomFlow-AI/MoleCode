"""Advanced stereo/1 features: planar/helical chirality, explicit E/Z reference
atoms and bond-level enhanced stereo groups.

Run from the repository root:

    python examples/10_advanced_stereo.py
"""

from rdkit import Chem
from rdkit.Chem import AllChem

from molecode.molecule import mermaid_to_mol, mol_to_mermaid, mol_to_smiles
from molecode.molecule import stereo


def show(title: str, doc: str) -> None:
    print("\n" + "-" * 78)
    print(title)
    print("-" * 78)
    for line in doc.splitlines():
        if any(key in line for key in ("_TET_", "_AX_", "_EZ_", "_PL_",
                                       "_HL_", "_SG_", "_ISO_", "===|")):
            print(line.strip())


def ez_explicit_references() -> None:
    """The base E/Z token has no reference atoms; add them only when needed."""
    mol = Chem.MolFromSmiles("C/C=C(/F)Br")  # first-neighbour rule flips this
    doc = mol_to_mermaid(mol, name="EZ", parity_mode="all")
    show("E/Z with explicit reference atoms (ambiguous substitution)", doc)
    back = mermaid_to_mol(doc)
    print("round-trip InChIKey equal:",
          stereo.inchi_key(mol) == stereo.inchi_key(back))


def bond_level_group() -> None:
    """Enhanced stereo can group axis members (atropisomer AND group)."""
    mol = Chem.MolFromSmiles("Cc1ccccc1-c1ccccc1C |wU:6.5|")
    rw = Chem.RWMol(mol)
    axis = [b.GetIdx() for b in mol.GetBonds() if b.GetStereo() not in
            (Chem.BondStereo.STEREONONE,)][0]
    group = Chem.CreateStereoGroup(Chem.StereoGroupType.STEREO_AND, rw, [], [axis])
    rw.SetStereoGroups([group])
    with_group = rw.GetMol()
    doc = mol_to_mermaid(with_group, name="Bip", parity_mode="all")
    show("Atropisomer with a bond-level AND group", doc)
    back = mermaid_to_mol(doc)
    print("group bonds preserved:",
          [b.GetIdx() for b in back.GetStereoGroups()[0].GetBonds()] == [axis])
    print("CXSMILES equal:",
          stereo.cx_smiles(stereo.prepare_mol(with_group))
          == stereo.cx_smiles(stereo.prepare_mol(back)))


def helical_chirality() -> None:
    """M/P helicity derived from a 3D conformer and carried by the document."""
    mol = Chem.MolFromSmiles("c1ccc2c(c1)ccc1c2ccc2ccccc21")  # hexahelicene
    AllChem.EmbedMolecule(mol, randomSeed=11)
    path = list(range(0, 8))
    sign = stereo.set_helical_chirality(mol, path)
    doc = mol_to_mermaid(mol, name="Helicene", parity_mode="all")
    show("Helical chirality (hl=M/P derived from 3D)", doc)
    back = mermaid_to_mol(doc)
    _planar, helical, _warnings = stereo.extra_stereo_records(back)
    print("derived sign:", sign, "| round-trip units:", len(helical))


def planar_chirality() -> None:
    """A planar unit (plane + two face references) derived from 3D coordinates."""
    mol = Chem.MolFromSmiles("CCCCC")
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, point in enumerate([(0, 0, 0), (1, 0, 0), (0, 1, 0),
                               (0.2, 0.2, 1.0), (0.8, 0.1, -1.0)]):
        conf.SetAtomPosition(i, point)
    mol.RemoveAllConformers()
    mol.AddConformer(conf)
    sign = stereo.set_planar_chirality(mol, [0, 1, 2], 3, 4)
    doc = mol_to_mermaid(mol, name="Planar", parity_mode="all")
    show("Planar chirality (toy 3D geometry)", doc)
    back = mermaid_to_mol(doc)
    planar, _helical, _warnings = stereo.extra_stereo_records(back)
    print("derived sign:", sign, "| round-trip units:", len(planar))


if __name__ == "__main__":
    print("=" * 78)
    print("MoleCode stereo/1 - advanced features")
    print("=" * 78)
    ez_explicit_references()
    bond_level_group()
    helical_chirality()
    planar_chirality()
    print("\nAll advanced stereo features exercised.")
