# MoleCode stereo extension (profile `stereo/1`)

> Design, grammar, invariants and API of the MoleCode stereo extension.  This page is the
> authoritative reference; the implementation lives in
> [`molecode/molecule/stereo.py`](../molecode/molecule/stereo.py).

## 1. Why a separate stereo layer?

MoleCode's base grammar (named atoms, typed bonds, explicit H counts, charge,
`===|E|`) is a **constitution** language.  Chiral chemistry does not fit
one fixed granularity:

| Stereo fact | Unary? | Relational? | Globally derived? |
| --- | --- | --- | --- |
| tetrahedral R/S | yes (one centre) | no | **yes** (CIP depends on the whole ligand graph) |
| pseudo-asymmetric r/s | yes | **yes** (only exists relative to the other centre) | yes |
| enhanced stereo AND / OR / ABS | no | **yes** (n-ary group) | no |
| axial chirality (atropisomer) | no (axis bond + 2 references) | **yes** | yes (aR/aS needs CIP priorities) |
| isotope | yes | no | no (but changes CIP) |
| E/Z | no (one bond, two references) | no | partly |

Two failure modes follow when everything is shoved into an atom-id suffix:

1. **Descriptor as truth.**  Absolute CIP is a global function.  A local edit
   changes remote labels (measured: a single terminal CH3 -> F edit rewrites
   ~30% of `_R/_S` suffixes in multi-centre molecules), so a
   representation that stores CIP as truth forces global re-labelling.
2. **Silent over-claim.**  Enhanced stereo AND/OR means "relative configuration
   known, absolute unknown".  Writing `_R/_S` for such atoms turns a
   *relative* claim into an *absolute* one.  The stock writer does exactly this,
   and full InChIKey cannot see it (InChI drops enhanced stereo groups).

The extension therefore separates **source truth** (local parity + relations)
from **derived descriptors** (CIP labels, E/Z), and keeps the base projection
**sound** (it never asserts more than the extension entails).

## 2. Layers

| Layer | Content | Stock parser sees it? |
| --- | --- | --- |
| **S0 base projection** | atoms, bonds, H counts, charge, `===|E|/Z|` | yes |
| **S1 local parity (source truth)** | `Pfx_TET_1{tet=Pfx_C_2; par=cw}` | no (ignored) |
| **S2 relations / higher-order objects** | stereo groups, axial units, isotopes | no (ignored) |
| **S3 derived caches** | atom-id suffix `_R/_S/_r/_s` (modern CIP), E/Z | yes (suffix) / yes (bond) |

S0 is what every existing MoleCode consumer reads.  The stock parser simply
skips `%%` profile comments, `{...}` object lines and
`-.->` member edges, so a stereo document is still a valid MoleCode
document (with the declared loss listed in §7).

## 3. Grammar

### 3.1 Header (only when S1/S2 is present)

    %% molecode-profile: stereo/1
    %% cip-impl: hanson2018 (rdkit.Chem.rdCIPLabeler)

`cip-impl` declares which CIP implementation produced the `_R/_S`
cache; a reader must validate the cache with that implementation.

### 3.2 Tetrahedral centres

* Absolute configuration: atom-id suffix `_R` / `_S` (unchanged base grammar).
* Pseudo-asymmetric: `_r` / `_s` (lowercase; invisible to the stock parser).
* Local parity record (source truth, V2 mode):

      Pfx_TET_1{tet=Pfx_C_2; par=ccw}

* Enhanced stereo AND/OR members get **only** the `TET` record and **no**
  absolute suffix in the base projection (invariant I1).

### 3.3 Enhanced stereo groups

      Pfx_SG_1{sg=AND:1; mem=Pfx_C_2,Pfx_C_4}
      Pfx_SG_2{sg=OR:1;  mem=Pfx_C_5,Pfx_C_6}
      Pfx_SG_3{sg=ABS:1; mem=Pfx_C_7}

Group types: `AND` (relative, absolute unknown), `OR` (mixture / racemic),
`ABS` (explicit absolute set; members also carry `_R/_S` suffixes).

### 3.4 Axial chirality (atropisomer)

      Pfx_AX_1{ax=cw; axis=Pfx_C_7,Pfx_C_8; ref1=Pfx_C_6; ref2=Pfx_C_9}

`axis` is the rotatable single bond, `ref1` a neighbour of its first atom
and `ref2` a neighbour of its second atom.  `ax=cw|ccw` is the RDKit
atropisomer sign in that reference frame; the writer resolves the frame so that
the molecule's canonical CXSMILES (`wU:` / `wD:`) is reproduced (RDKit does
not persist reference atoms, so a floating sign alone is ambiguous).

### 3.5 Isotopes

      Pfx_ISO_1{iso=13; atom=Pfx_C_1}

Isotopes are applied **before** CIP validation (invariant I4): isotopic
substitution can create a stereocentre and always changes CIP priorities.

## 4. The local-parity reference frame

`par=cw` / `par=ccw` is the RDKit tetrahedral chiral tag of the atom with
respect to **the order in which its incident edges appear in the document**
(RDKit / SMILES convention: first neighbour in bond order, implicit H last).
Because every bond is one line and the reader adds bonds in document order, the
frame is explicit, reproducible and round-trip stable:

* the writer emits bonds in its molecule's bond order;
* the reader inserts them in the same order;
* hence writer and reader agree on "first neighbour" without any hidden atom
  numbering.

**Editing rule:** append new bonds for an existing centre rather than inserting
them before existing incident edges; the CIP cache (S3) is re-derived on read,
so a broken frame is *detected* (suffix vs parity mismatch) instead of silently
flipping the stereocentre.

This is intentionally different from SMILES `@/@@`, whose reference order is
implicit in the string, and from ConfSeq's pseudo-chirality, whose CW/CCW also
depends on the (randomised) SMILES neighbour order.

## 5. Examples

### 5.1 Relative configuration (threonine, AND group + ABS centre)

    graph TB
        %% molecode-profile: stereo/1
        %% cip-impl: hanson2018 (rdkit.Chem.rdCIPLabeler)
        subgraph Thr["2-amino-3-hydroxybutanoic acid"]
            Thr_N_1[NH2]
            Thr_C_1[CH]
            Thr_C_2[CH3]
            Thr_C_3_S[CH]
            Thr_O_1[OH]
            Thr_C_4[C]
            Thr_O_2[O]
            Thr_O_3[OH]

            Thr_N_1 --- Thr_C_1
            Thr_C_1 --- Thr_C_2
            Thr_C_1 --- Thr_C_3_S
            Thr_C_3_S --- Thr_O_1
            Thr_C_3_S --- Thr_C_4
            Thr_C_4 === Thr_O_2
            Thr_C_4 --- Thr_O_3

            %% ---- stereo layer (stereo/1) ----
            Thr_TET_1{tet=Thr_C_1; par=cw}
            Thr_TET_2{tet=Thr_C_3_S; par=ccw}
            Thr_SG_1{sg=AND:1; mem=Thr_C_1,Thr_C_3_S}
        end

Note that `Thr_C_1` has no absolute suffix: the base projection is sound (it
does not claim that the relative pair is a particular enantiomer).

### 5.2 Atropisomer

    Biaryl_AX_1{ax=cw; axis=Biaryl_C_7,Biaryl_C_8; ref1=Biaryl_C_6; ref2=Biaryl_C_9}

### 5.3 Isotopically labelled centre

    Lac_ISO_1{iso=13; atom=Lac_C_1}
    Lac_TET_1{tet=Lac_C_2_S; par=ccw}

## 6. Invariants

| Id | Invariant | Checked by |
| --- | --- | --- |
| I1 | Soundness: no absolute `_R/_S` on AND/OR members in the base projection | `stereo.validate_document` |
| I2 | Every suffix equals the modern CIP descriptor of the final molecule | reader (re-derived) |
| I3 | Each stereo element belongs to at most one group; ABS and AND/OR are disjoint | apply / validate |
| I4 | Isotopes applied before CIP computation/validation | reader order |
| I5 | Axial `ref1/ref2` are neighbours of the two axis atoms | apply / validate |
| I6 | Full InChIKey is **not** a valid judge when AND/OR or axial is present | validate_document |
| I7 | `-.->` member edges must reference declared nodes | base parser contract |
| I8 | Re-serialization is stable up to Kekule resonance | unit tests |

## 7. Round-trip guarantees and judges

| Information | Base projection (stock parser) | Extension-aware reader |
| --- | --- | --- |
| tetrahedral R/S | preserved (suffix) | preserved + validated |
| tetrahedral r/s | lost (suffix invisible) | preserved |
| E/Z | preserved (`===|E|`) | preserved |
| enhanced AND/OR/ABS | **lost** (sound: unspecified) | preserved |
| axial | **lost** | preserved |
| isotope | **lost** | preserved |
| local parity | lost | preserved |

Because InChI drops enhanced stereo groups and is blind to axial chirality, the
judge must be *declared*: use full InChIKey + CXSMILES + StereoGroup signatures
together.  `stereo.validate_document` emits the I6 warning when a document
contains AND/OR or axial objects.

## 8. API

    from molecode.molecule import (
        mol_to_mermaid, mermaid_to_mol,
        prepare_mol, modern_cip_map, extract_stereo_info,
        render_stereo_layer, parse_stereo_layer, apply_stereo_layer,
        stereo_fingerprint, stereo_coverage, validate_document,
        PROFILE_ID, CIP_IMPL,
    )

    doc = mol_to_mermaid(mol, name="Thr", stereo_layer=True, parity_mode="all")  # opt-in source-truth mode
    mol2 = mermaid_to_mol(doc)                     # extension-aware
    base = mermaid_to_mol(doc, parse_stereo_layer=False)   # stock semantics

    info = extract_stereo_info(prepare_mol(mol), cip=modern_cip_map(mol))
    report = validate_document(doc)                # static invariant checks
    fp = stereo_fingerprint(mol)                   # CIP / groups / axial / isotopes

`parity_mode` controls how much source truth is written:

| Mode | Stored | Use |
| --- | --- | --- |
| `"group_only"` (default) | local parity only where CIP cannot stand alone (AND/OR members, pseudo-asymmetric) | compact documents; plain R/S uses the CIP cache |
| `"all"` (opt-in) | local parity for every tetrahedral centre + relations | edit-locality research mode; a remote local edit never rewrites source records |
| `"off"` | no S1 records (relations only) | ablations |

## 9. Measured behaviour (144-molecule stereo benchmark)

Corpus: 144 molecules / 14 categories (tetrahedral, pseudo-asymmetric, isotope
chirality, E/Z, axial, enhanced AND/OR/ABS/mixed, unspecified, ring/bridged,
charged hetero) with an independent RDKit oracle.

| Variant | full InChIKey | CXSMILES | StereoGroups | axial | isotopes | over-claim |
| --- | --- | --- | --- | --- | --- | --- |
| V0 stock MoleCode 0.1.0 | 122/144 | 85/144 | 0/17 | 0/20 | 0/10 | 9.5% |
| V1 (this extension, group_only) | 144/144 | 143/144 | 17/17 | 20/20 | 10/10 | 0 |
| V2 (this extension, all) | 144/144 | 143/144 | 17/17 | 20/20 | 10/10 | 0 |

The one CXSMILES difference is a meso compound (myo-inositol): two SMILES strings
that are mirror images of the same achiral molecule; the InChIKey is identical.
The token cost of the stereo layer and the edit-locality experiment are
reported in the accompanying study; `tests/test_stereo.py` covers the
round-trip and invariant checks in this repository.

## 10. Limitations

* **Planar / helical chirality** are not represented yet (no object type).
* **Bond-level enhanced stereo** (AND/OR of atropisomers) is rejected with a
  warning; only atom-level groups are round-tripped.
* **Kekule resonance**: aromatic documents are not byte-identical after a
  read/write cycle (the base projection already has this property).
* **Conformation**: this profile encodes configuration, not conformation.  A
  separate geometry layer (bond lengths / angles / dihedrals, cf. ConfSeq) is
  orthogonal and does not replace configuration.
* The local-parity frame is document-order based: inserting a new incident edge
  *before* existing ones for a centre changes its reference frame.  The CIP
  cache validation catches the resulting mismatch (I2) instead of failing
  silently; append-only edits are safe.

## 11. References

* Daylight SMILES: neighbour-order dependence of `@/@@`.
* CXSMILES / OpenSMILES: enhanced stereo `|a:|`, `|&1:|`, `|o1:|`.
* MDL MOL V3000: `CFG` and `STERAC / STEREL / STEABS` collections.
* InChI Technical FAQ: stereochemistry is treated as absolute; enhanced stereo is
  dropped; axial chirality is not covered.
* Hanson et al., *J. Chem. Inf. Model.* 2018, 58, 1755-1765 (CIP implementation
  and its revision; the `hanson2018` cache refers to this).
* ConfSeq (bioRxiv 2025 / Nat. Mach. Intell. 2026): 3D internal-coordinate
  language; its pseudo-chirality CW/CCW is neighbour-order dependent and only
  applies to non-stereogenic centres, so it is complementary, not a substitute.
