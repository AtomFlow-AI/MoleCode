"""Unit tests for the MoleCode stereo/1 extension.

Run with:
    set PYTHONPATH=<repo-root>
    python -m unittest discover -s tests -v

The tests use only RDKit + the local molecode package; they do not require
network access or an LLM.
"""

import unittest

from rdkit import Chem

from molecode.molecule import mol_to_mermaid, mermaid_to_mol
from molecode.molecule import stereo as S


def _mol(smiles):
    m = Chem.MolFromSmiles(smiles)
    assert m is not None, smiles
    return m


def _roundtrip(smiles, parity_mode="all", name="Mol"):
    mol = _mol(smiles)
    doc = mol_to_mermaid(mol, name=name, parity_mode=parity_mode)
    out = mermaid_to_mol(doc)
    assert out is not None, "reader returned None"
    return mol, out, doc


def _group_sig(mol):
    return sorted(
        (S.group_type_name(g.GetGroupType()), tuple(sorted(a.GetIdx() for a in g.GetAtoms())))
        for g in mol.GetStereoGroups()
    )


def _axial_sig(mol):
    return sorted(
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx(), str(b.GetStereo()))
        for b in mol.GetBonds()
        if b.GetStereo() in (Chem.BondStereo.STEREOATROPCW, Chem.BondStereo.STEREOATROPCCW)
    )


def _inchi_key(mol):
    return S.inchi_key(mol)


class TestExtraction(unittest.TestCase):
    def test_absolute_and_local_parity(self):
        mol = _mol("C[C@H](O)C(=O)O")
        info = S.extract_stereo_info(S.prepare_mol(mol), cip=S.modern_cip_map(mol))
        self.assertEqual(info.abs_suffix, {1: "S"})
        self.assertIn(1, info.local_par)
        self.assertIn(info.local_par[1], ("cw", "ccw"))

    def test_pseudo_asymmetric_lowercase(self):
        mol = _mol("C[C@H](O)[C@H](O)[C@@H](C)O")
        info = S.extract_stereo_info(S.prepare_mol(mol), cip=S.modern_cip_map(mol))
        self.assertIn("s", info.abs_suffix.values())

    def test_isotope_descriptor_conflict(self):
        # legacy CIP says R, modern CIP (and InChI) say S -- the cache must use
        # the modern descriptor.
        mol = _mol("[13CH3][C@H](O)C(=O)O")
        legacy = S.legacy_cip_map(mol)
        modern = S.modern_cip_map(mol)
        self.assertNotEqual(legacy, modern)
        info = S.extract_stereo_info(S.prepare_mol(mol), cip=modern)
        self.assertEqual(info.abs_suffix[1], "S")
        self.assertEqual(info.isotopes[0], 13)

    def test_group_members_have_no_absolute_suffix(self):
        mol = _mol("N[C@@H](C)[C@H](O)C(=O)O |&1:1,3|")
        info = S.extract_stereo_info(S.prepare_mol(mol), cip=S.modern_cip_map(mol))
        self.assertEqual(info.abs_suffix, {})
        self.assertEqual(sorted(info.local_par), [1, 3])
        self.assertEqual(len(info.groups), 1)
        self.assertEqual(info.groups[0].group_type, "AND")

    def test_axial_extraction_resolves_cxsmiles_frame(self):
        mol = _mol("Cc1ccccc1-c1ccccc1C |wU:6.5|")
        info = S.extract_stereo_info(S.prepare_mol(mol), cip=S.modern_cip_map(mol))
        self.assertEqual(len(info.axial), 1)
        rec = info.axial[0]
        rebuilt = S._mol_with_axial(Chem.Mol(mol), info.axial)
        self.assertEqual(S._cx_axial_tokens(rebuilt), S._cx_axial_tokens(mol))


class TestRoundTrip(unittest.TestCase):
    CASES = [
        "C[C@H](O)C(=O)O",
        "N[C@@H](C)C(=O)O",
        "C[C@H](O)[C@H](O)[C@@H](C)O",
        "N[C@@H]1CC[C@H](O)CC1",
        "C[C@H](O)[C@@H](N)C(=O)O",
    ]

    def test_tetrahedral_both_modes(self):
        for smi in self.CASES:
            for mode in ("all", "group_only"):
                with self.subTest(smiles=smi, mode=mode):
                    _in, out, _doc = _roundtrip(smi, parity_mode=mode)
                    self.assertEqual(_inchi_key(_in), _inchi_key(out))

    def test_ez(self):
        for smi in ["C/C=C/C", "C/C=C\\C", "F/C=C/Cl"]:
            with self.subTest(smiles=smi):
                _in, out, _doc = _roundtrip(smi)
                self.assertEqual(_inchi_key(_in), _inchi_key(out))

    def test_enhanced_stereo_groups(self):
        for smi in [
            "N[C@@H](C)[C@H](O)C(=O)O |&1:1,3|",
            "N[C@@H](C)[C@H](O)C(=O)O |o1:1,3|",
            "N[C@@H](C)[C@H](O)C(=O)O |a:1|",
        ]:
            with self.subTest(smiles=smi):
                _in, out, doc = _roundtrip(smi)
                self.assertEqual(_group_sig(_in), _group_sig(out))
                self.assertEqual(S.cx_smiles(S.prepare_mol(_in)),
                                 S.cx_smiles(S.prepare_mol(out)))
                # invariant I1: no absolute suffix on AND/OR members
                self.assertNotIn("_R[", doc.split("stereo layer")[-1])
                self.assertNotIn("_S[", doc.split("stereo layer")[-1])

    def test_isotopes(self):
        for smi in ["[13CH3][C@H](O)C(=O)O", "[2H][C@H](F)Cl",
                    "[2H]C([2H])([2H])O"]:
            with self.subTest(smiles=smi):
                _in, out, _doc = _roundtrip(smi)
                self.assertEqual(_inchi_key(_in), _inchi_key(out))

    def test_axial(self):
        for smi in ["Cc1ccccc1-c1ccccc1C |wU:6.5|",
                    "Cc1ccccc1-c1ccccc1C |wD:6.5|"]:
            with self.subTest(smiles=smi):
                _in, out, doc = _roundtrip(smi)
                self.assertEqual(_axial_sig(_in), _axial_sig(out))
                self.assertIn("_AX_", doc)


class TestBaseProjection(unittest.TestCase):
    def test_and_members_are_unspecified_in_base_projection(self):
        smi = "N[C@@H](C)[C@H](O)C(=O)O |&1:1,3|"
        doc = mol_to_mermaid(_mol(smi), name="Thr")
        parser = S.__dict__.get("MermaidMolParser")
        from molecode.molecule.mermaid_to_rdkit import MermaidMolParser
        base = MermaidMolParser(parse_stereo_layer=False).parse_mermaid_graph(doc)
        self.assertIsNotNone(base)
        self.assertNotIn("@", Chem.MolToSmiles(base))
        stripped = Chem.MolFromSmiles(Chem.MolToSmiles(_mol(smi)))
        self.assertEqual(Chem.MolToSmiles(base, isomericSmiles=False),
                         Chem.MolToSmiles(stripped, isomericSmiles=False))

    def test_extension_lines_are_ignored_by_base_parser(self):
        from molecode.molecule.mermaid_to_rdkit import MermaidMolParser
        smi = "Cc1ccccc1-c1ccccc1C |wU:6.5|"
        doc = mol_to_mermaid(_mol(smi), name="Biphenyl")
        base = MermaidMolParser(parse_stereo_layer=False).parse_mermaid_graph(doc)
        self.assertIsNotNone(base)
        self.assertEqual(Chem.MolToSmiles(base, isomericSmiles=False),
                         "Cc1ccccc1-c1ccccc1C")

    def test_profile_header_present_only_with_extension(self):
        # default (group_only) writes no TET record for a plain absolute centre
        # -> no extension object -> no profile header (pure base document)
        simple = mol_to_mermaid(_mol("C[C@H](O)C(=O)O"), name="Lac")
        self.assertEqual(simple, mol_to_mermaid(_mol("C[C@H](O)C(=O)O"),
                                                name="Lac", parity_mode="group_only"))
        self.assertNotIn("molecode-profile", simple)
        # parity=all stores the local parity source record -> header present
        v2 = mol_to_mermaid(_mol("C[C@H](O)C(=O)O"), name="Lac", parity_mode="all")
        self.assertIn("molecode-profile: stereo/1", v2)
        self.assertIn("cip-impl: hanson2018", v2)
        group = mol_to_mermaid(_mol("N[C@@H](C)[C@H](O)C(=O)O |&1:1,3|"), name="Thr")
        self.assertIn("molecode-profile: stereo/1", group)


class TestValidation(unittest.TestCase):
    def test_over_claim_detected(self):
        text = "\n".join([
            "Thr_SG_1{sg=AND:1; mem=Thr_C_1,Thr_C_2}",
        ])
        report = S.validate_document(text, suffix_map={"Thr_C_1": "S"})
        self.assertFalse(report["ok"])
        self.assertTrue(any("I1" in i for i in report["issues"]))

    def test_judge_warning_for_groups_and_axial(self):
        text = "M_SG_1{sg=OR:1; mem=M_C_1,M_C_2}"
        report = S.validate_document(text)
        self.assertTrue(any("I6" in i for i in report["issues"]))

    def test_malformed_axial_ref_reported(self):
        layer = S.ParsedStereoLayer()
        layer.axial.append(("cw", "A_C_1", "A_C_2", "A_C_9", "A_C_10"))
        base = Chem.MolFromSmiles("Cc1ccccc1-c1ccccc1C")
        idx = {f"A_C_{i + 1}": i for i in range(base.GetNumAtoms())}
        _mol_out, issues = S.apply_stereo_layer(base, layer, idx)
        self.assertTrue(issues)


class TestIdempotence(unittest.TestCase):
    def test_doc_serialization_is_stable(self):
        # Non-aromatic molecules: byte-identical re-serialization.
        for smi in ["C[C@H](O)[C@@H](N)C(=O)O",
                    "N[C@@H](C)[C@H](O)C(=O)O |&1:1,3|",
                    "[13CH3][C@H](O)C(=O)O"]:
            with self.subTest(smiles=smi):
                doc1 = mol_to_mermaid(_mol(smi), name="M", parity_mode="all")
                mol2 = mermaid_to_mol(doc1)
                doc2 = mol_to_mermaid(mol2, name="M", parity_mode="all")
                self.assertEqual(doc1, doc2)

    def test_aromatic_docs_stable_up_to_kekule(self):
        # Aromatic rings admit several valid Kekule forms; the writer is not
        # required to reproduce the same resonance structure, only the same
        # molecule (documented limitation of the base projection).
        smi = "Cc1ccccc1-c1ccccc1C |wU:6.5|"
        doc1 = mol_to_mermaid(_mol(smi), name="M", parity_mode="all")
        mol2 = mermaid_to_mol(doc1)
        doc2 = mol_to_mermaid(mol2, name="M", parity_mode="all")
        self.assertEqual(_inchi_key(_mol(smi)), _inchi_key(mol2))
        self.assertEqual(S.cx_smiles(S.prepare_mol(_mol(smi))),
                         S.cx_smiles(S.prepare_mol(mol2)))
        self.assertEqual(_axial_sig(_mol(smi)), _axial_sig(mol2))


class TestPublicAPI(unittest.TestCase):
    def test_fingerprint_keys(self):
        fp = S.stereo_fingerprint(_mol("C[C@H](O)C(=O)O"))
        for key in ("smiles", "inchi_key", "cxsmiles", "cip", "groups",
                    "axial", "isotopes"):
            self.assertIn(key, fp)

    def test_coverage(self):
        cov = S.stereo_coverage(_mol("N[C@@H](C)[C@H](O)C(=O)O |&1:1,3|"))
        self.assertEqual(cov["groups"], 1)
        self.assertGreaterEqual(cov["centres"], 2)


class TestValidatorRegressions(unittest.TestCase):
    def test_reader_enforces_i1_soundness(self):
        smi = "N[C@@H](C)[C@H](O)C(=O)O |&1:1,3|"
        doc = mol_to_mermaid(_mol(smi), name="Thr")
        # tamper: give an AND member an absolute suffix (rename it everywhere,
        # like an editor would do when regenerating the graph)
        tampered = doc.replace("Thr_C_1", "Thr_C_1_S")
        from molecode.molecule.mermaid_to_rdkit import MermaidMolParser
        parser = MermaidMolParser()
        mol = parser.parse_mermaid_graph(tampered)
        self.assertIsNotNone(mol)
        self.assertTrue(any("I1" in w for w in parser.warnings))
        # the extension-aware reader ignores the tampered suffix: the group is
        # still relative and the molecule is unchanged
        self.assertTrue(mol.GetStereoGroups())
        self.assertEqual(_inchi_key(_mol(smi)), _inchi_key(mol))

    def test_diene_ez_cycle_stable(self):
        smi = "C/C=C\\C=C/C"
        doc1 = mol_to_mermaid(_mol(smi), name="M", parity_mode="all")
        mol2 = mermaid_to_mol(doc1)
        doc2 = mol_to_mermaid(mol2, name="M", parity_mode="all")
        self.assertEqual(doc1, doc2)
        self.assertEqual(_inchi_key(_mol(smi)), _inchi_key(mol2))

    def test_duplicate_tet_flagged(self):
        text = "\n".join([
            "M_TET_1{tet=M_C_1; par=cw}",
            "M_TET_2{tet=M_C_1; par=ccw}",
        ])
        report = S.validate_document(text)
        self.assertTrue(any("duplicate TET" in i for i in report["issues"]))

    def test_group_overlap_flagged(self):
        text = "\n".join([
            "M_SG_1{sg=AND:1; mem=M_C_1,M_C_2}",
            "M_SG_2{sg=OR:1; mem=M_C_1}",
        ])
        report = S.validate_document(text)
        self.assertTrue(any("group conflict" in i for i in report["issues"]))


if __name__ == "__main__":
    unittest.main()
