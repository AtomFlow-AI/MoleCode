"""Stereo round-trip demo for the MoleCode stereo/1 extension.

Shows, for a small set of stereochemically demanding molecules, that
RDKit -> MoleCode -> RDKit preserves the full stereo fingerprint:

* absolute tetrahedral centres (R/S) and pseudo-asymmetric centres (r/s),
* E/Z double bonds,
* enhanced stereo groups (AND / OR / ABS),
* axial chirality (atropisomers),
* isotope-labelled stereocentres.

Run from the repository root:

    python examples/09_stereo_roundtrip.py
"""

from rdkit import Chem

from molecode.molecule import mol_to_mermaid, mermaid_to_mol, mol_to_smiles
from molecode.molecule import stereo


CASES = [
    ("(S)-lactic acid", "C[C@H](O)C(=O)O"),
    ("pseudo-asymmetric pentane-2,3,4-triol",
     "C[C@H](O)[C@H](O)[C@@H](C)O"),
    ("13C-labelled lactic acid", "[13CH3][C@H](O)C(=O)O"),
    ("(E)-2-butene", "C/C=C/C"),
    ("(Z)-2-butene", "C/C=C\\C"),
    ("threonine, relative configuration (AND group)",
     "N[C@@H](C)[C@H](O)C(=O)O |&1:1,3|"),
    ("racemic mixture (OR group)",
     "N[C@@H](C)[C@H](O)C(=O)O |o1:1,3|"),
    ("2,2'-dimethylbiphenyl atropisomer",
     "Cc1ccccc1-c1ccccc1C |wU:6.5|"),
]


def main() -> None:
    print("=" * 78)
    print("MoleCode stereo/1 round-trip demo (parity_mode 'all', opt-in source truth)")
    print("=" * 78)
    all_ok = True
    for name, smiles in CASES:
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, smiles
        doc = mol_to_mermaid(mol, name=name.split(",")[0][:14], parity_mode="all")
        back = mermaid_to_mol(doc)
        fp_in = stereo.stereo_fingerprint(mol)
        fp_out = stereo.stereo_fingerprint(back)
        keys = ("inchi_key", "cxsmiles", "cip", "groups", "axial", "isotopes")
        ok = all(fp_in[k] == fp_out[k] for k in keys)
        all_ok = all_ok and ok
        print(f"\n[{ 'OK ' if ok else 'FAIL' }] {name}")
        print(f"        in : {mol_to_smiles(mol)}")
        print(f"        out: {mol_to_smiles(back)}")
        if not ok:
            for k in keys:
                if fp_in[k] != fp_out[k]:
                    print(f"        mismatch {k}: {fp_in[k]!r} != {fp_out[k]!r}")

    print("\n" + "=" * 78)
    print("A relative-configuration document (extension layer at the bottom):")
    print("=" * 78)
    demo = Chem.MolFromSmiles("N[C@@H](C)[C@H](O)C(=O)O |&1:1,3|")
    print(mol_to_mermaid(demo, name="Thr"))
    print("\nAll cases round-trip:", all_ok)


if __name__ == "__main__":
    main()
