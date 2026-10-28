# MoleCode stereo extension (profile `stereo/1`)

> Design, grammar, invariants and API of the MoleCode stereo extension.
> Implementation: [`molecode/molecule/stereo.py`](../molecode/molecule/stereo.py).

## 1. Why a separate stereo layer?

MoleCode's base grammar (named atoms, typed bonds, explicit H counts, charge,
`===|E|`) is a **constitution** language.  Chiral chemistry does not fit
one fixed granularity:

| Stereo fact | Unary? | Relational? | Globally derived? |
| --- | --- | --- | --- |
| tetrahedral R/S | yes (one centre) | no | **yes** (CIP depends on the whole ligand graph) |
| pseudo-asymmetric r/s | yes | **yes** (only exists relative to the other centre) | yes |
| enhanced stereo AND / OR / ABS (atoms) | no | **yes** (n-ary group) | no |
| enhanced stereo on an axis (atoms + bonds) | no | **yes** | no |
| axial chirality (atropisomer) | no (axis bond + 2 references) | **yes** | yes (aR/aS needs CIP priorities) |
| planar chirality | no (plane + 2 face references) | **yes** | yes |
| helical chirality | no (helix path) | **yes** | yes |
| E/Z reference atoms | no (one bond, two references) | no | partly |
| isotope | yes | no | no (but changes CIP) |

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
| **S2 relations / higher-order objects** | stereo groups (atom and/or axis members), axial units, E/Z reference atoms, planar units, helical units, isotopes | no (ignored) |
| **S3 derived caches** | atom-id suffix `_R/_S/_r/_s` (modern CIP), E/Z token | yes (suffix) / yes (bond) |

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

### 3.3 Enhanced stereo groups (atom and axis members)

      Pfx_SG_1{sg=AND:1; mem=Pfx_C_2,Pfx_C_4}
      Pfx_SG_2{sg=OR:1;  mem=Pfx_C_5,Pfx_C_6}
      Pfx_SG_3{sg=ABS:1; mem=Pfx_C_7}
      Pfx_SG_4{sg=AND:2; axes=Pfx_AX_1,Pfx_AX_2}          # atropisomers
      Pfx_SG_5{sg=OR:1;  mem=Pfx_C_8; axes=Pfx_AX_3}      # mixed members

Group types: `AND` (relative, absolute unknown), `OR` (mixture / racemic),
`ABS` (explicit absolute set; atom members also carry `_R/_S` suffixes).
`axes` references AX object ids (§3.4); an axis can belong to at most one group.

### 3.4 Axial chirality (atropisomer)

      Pfx_AX_1{ax=cw; axis=Pfx_C_7,Pfx_C_8; ref1=Pfx_C_6; ref2=Pfx_C_9}

`axis` is the rotatable single bond, `ref1` a neighbour of its first atom
and `ref2` a neighbour of its second atom.  `ax=cw|ccw` is the RDKit
atropisomer sign in that reference frame; the writer resolves the frame so that
the molecule's canonical CXSMILES (`wU:` / `wD:`) is reproduced (RDKit does
not persist reference atoms, so a floating sign alone is ambiguous).

### 3.5 Explicit E/Z reference atoms (only when needed)

      Pfx_EZ_1{ez=Z; bond=Pfx_C_2,Pfx_C_3; ref1=Pfx_C_1; ref2=Pfx_Br_1}

The base token (`===|E|` / `===|Z|`) remains unchanged; the reader's
documented default is "first neighbour of each end".  For some substitution
patterns that default would flip the geometry (e.g. `C/C=C(/F)Br`).  The
writer therefore checks the round-trip with the reader's exact algorithm and
emits an E/Z object only when the default is ambiguous.  Reader: the object's
`ref1/ref2` override the default reference atoms.

### 3.6 Planar chirality

      Pfx_PL_1{pl=+; plane=Pfx_C_1,Pfx_C_2,Pfx_C_3,Pfx_C_4; ref1=Pfx_C_5; ref2=Pfx_C_6}

`plane` contains at least three atoms (the first three define the normal
`n = (p1-p0) x (p2-p0)`); `ref1` / `ref2` are substituents on the two faces.
`pl` is a local sign (`+`/`-`; `pR`/`pS` accepted as a label) computed from 3D
coordinates as `sign((r1-p0).n - (r2-p0).n)`.  Mirroring the molecule or swapping
the two references flips the sign.

RDKit has no native planar-chirality stereo, so the record is stored in the mol
property `_molecode_stereo_extra` and re-emitted on write; doc -> mol -> doc
is lossless for it.  Use `stereo.set_planar_chirality(mol, plane, ref1, ref2)` to
attach a unit (the sign is derived from the conformer when omitted).

### 3.7 Helical chirality

      Pfx_HL_1{hl=M; path=Pfx_C_1,Pfx_C_2,Pfx_C_3,Pfx_C_4,Pfx_C_5,Pfx_C_6}

`path` lists at least four consecutive atoms of the helix backbone.  `hl=M|P`
is computed from 3D coordinates as the sign of the sum of consecutive triple
products; mirroring flips M <-> P.  Like planar chirality it is stored in the
mol property `_molecode_stereo_extra`.  Use
`stereo.set_helical_chirality(mol, path)` to attach a unit.

### 3.8 Isotopes

      Pfx_ISO_1{iso=13; atom=Pfx_C_1}

Isotopes are applied **before** CIP validation (invariant I4): isotopic
substitution can create a stereocentre and always changes CIP priorities.

## 4. The local-parity reference frame

`par=cw` / `par=ccw` is the RDKit tetrahedral chiral tag of the atom with
respect to the **persistent atom-id order of its explicit neighbours**: the
neighbours are sorted by (prefix, element, counter) as written in the document,
with the implicit hydrogen treated as a fixed ligand.  The writer converts the
RDKit tag into this frame; the reader converts back.

Compared with the earlier "order of incident edges in the document" convention:

* **edit stability** — reordering bond lines, inserting new bonds elsewhere, or
  appending atoms does not change the meaning of an existing `par` record;
* **no hidden numbering** — the frame is exactly the ids the document already
  shows (no dependence on line order);
* the CIP cache (S3) is still re-derived on read, so a *chemically* invalid edit
  is detected instead of silently flipping a centre.

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

### 5.2 Atropisomer with a bond-level AND group

    Biaryl_AX_1{ax=cw; axis=Biaryl_C_7,Biaryl_C_8; ref1=Biaryl_C_6; ref2=Biaryl_C_9}
    Biaryl_SG_1{sg=AND:1; axes=Biaryl_AX_1}

### 5.3 E/Z with explicit references

    EZ_EZ_1{ez=Z; bond=EZ_C_2,EZ_C_3; ref1=EZ_C_1; ref2=EZ_Br_1}

### 5.4 Planar and helical units

    Pl_PL_1{pl=+; plane=Pl_C_1,Pl_C_2,Pl_C_3; ref1=Pl_C_4; ref2=Pl_C_5}
    Hl_HL_1{hl=P; path=Hl_C_1,Hl_C_2,Hl_C_3,Hl_C_4,Hl_C_5,Hl_C_6}

## 6. Invariants

| Id | Invariant | Checked by |
| --- | --- | --- |
| I1 | Soundness: no absolute `_R/_S` on AND/OR members in the base projection | `validate_document` (+ reader enforcement) |
| I2 | Every suffix equals the modern CIP descriptor; no duplicate/contradictory local records | reader re-derivation / validate |
| I3 | Each stereo element (atom or axis) belongs to at most one group; ABS and AND/OR are disjoint | apply / validate |
| I4 | Isotopes applied before CIP computation/validation | reader order |
| I5 | Axial `ref1/ref2` are neighbours of the two axis atoms | apply / validate |
| I6 | Full InChIKey is **not** a valid judge when AND/OR, axial, planar or helical chirality is present | validate_document |
| I7 | `-.->` member edges must reference declared nodes | base parser contract |
| I8 | Re-serialization is byte-stable (including aromatic Kekule forms) | unit tests |
| I9 | Planar record: >=3 plane atoms, both references present, sign in {+,-,pR,pS} | validate / apply |
| I10 | Helical record: >=4 path atoms, sign in {M,P} | validate / apply |
| I11 | A group member (atom or axis) appears in exactly one group | validate / apply |
| I12 | A group's `axes` reference exists and resolves to an AX object | validate / apply |

## 7. Round-trip guarantees and judges

| Information | Base projection (stock parser) | Extension-aware reader |
| --- | --- | --- |
| tetrahedral R/S | preserved (suffix) | preserved + validated |
| tetrahedral r/s | lost (suffix invisible) | preserved |
| E/Z | preserved (`===|E|`) | preserved + explicit refs when ambiguous |
| enhanced AND/OR/ABS (atoms) | **lost** (sound: unspecified) | preserved |
| enhanced AND/OR/ABS (axes) | **lost** (sound: unspecified) | preserved |
| axial | **lost** | preserved |
| planar / helical | **lost** | preserved (mol property; RDKit has no native representation) |
| isotope | **lost** | preserved |
| local parity | lost | preserved (edit-stable frame) |

Because InChI drops enhanced stereo groups and is blind to axial/planar/helical
chirality, the judge must be *declared*: use full InChIKey + CXSMILES +
StereoGroup signatures together.  `stereo.validate_document` emits the I6
warning when a document contains any of those elements.

## 8. API

    from molecode.molecule import (
        mol_to_mermaid, mermaid_to_mol,
        prepare_mol, modern_cip_map, extract_stereo_info,
        render_stereo_layer, parse_stereo_layer, apply_stereo_layer,
        stereo_fingerprint, stereo_coverage, validate_document,
        canonical_kekulize, set_double_bond_stereo,
        planar_chirality_sign, helical_chirality_sign,
        set_planar_chirality, set_helical_chirality,
        read_extra_stereo, write_extra_stereo, extra_stereo_records,
        PROFILE_ID, CIP_IMPL,
    )

    doc = mol_to_mermaid(mol, name="Thr", parity_mode="all")
    mol2 = mermaid_to_mol(doc)                     # extension-aware
    base = mermaid_to_mol(doc, parse_stereo_layer=False)   # stock semantics

    info = extract_stereo_info(prepare_mol(mol), cip=modern_cip_map(mol))
    report = validate_document(doc)                # static invariant checks
    fp = stereo_fingerprint(mol)                   # CIP / groups / axial / isotopes / planar / helical

    # attach geometry-derived higher-order chirality
    set_planar_chirality(mol, plane_atoms=[...], ref1=..., ref2=...)
    set_helical_chirality(mol, path_atoms=[...])

`parity_mode` controls how much source truth is written:

| Mode | Stored | Use |
| --- | --- | --- |
| `"group_only"` (default) | local parity only where CIP cannot stand alone (AND/OR members, pseudo-asymmetric) | compact documents; plain R/S uses the CIP cache |
| `"all"` (opt-in) | local parity for every tetrahedral centre + relations | edit-locality research mode; a remote local edit never rewrites source records |
| `"off"` | no S1 records (relations only) | ablations |

## 9. Measured behaviour

### 9.1 Chiral core corpus (144 molecules / 14 categories)

Independent RDKit oracle; official MoleCode 0.1.0 as the stock base parser.

| Variant | full InChIKey | CXSMILES | StereoGroups | axial | isotopes | over-claim |
| --- | --- | --- | --- | --- | --- | --- |
| V0 stock MoleCode 0.1.0 | 122/144 | 85/144 | 0/17 | 0/20 | 0/10 | 9.5% |
| V1 (this extension, group_only) | 144/144 | 143/144 | 17/17 | 20/20 | 10/10 | 0 |
| V2 (this extension, all) | 144/144 | 143/144 | 17/17 | 20/20 | 10/10 | 0 |

The one CXSMILES difference is a meso compound (myo-inositol): two SMILES strings
that are mirror images of the same achiral molecule; the InChIKey is identical.

### 9.2 Complex natural products (16 molecules, up to 21 centres / 101 heavy atoms)

Includes paclitaxel, docetaxel, erythromycin, azithromycin, digoxin, reserpine,
morphine, colchicine, artemisinin, sirolimus, cyclosporin A, vancomycin,
ivermectin, bryostatin 1, strychnine and quinine (PubChem isomeric SMILES).

| Variant | full InChIKey | CXSMILES |
| --- | --- | --- |
| V0 stock MoleCode 0.1.0 | 14/16 | 14/16 |
| V1 / V2 (this extension) | **16/16** | **16/16** |

The two V0 failures (ivermectin, bryostatin 1) are macrocyclic E/Z losses; the
explicit E/Z reference atoms in §3.5 fix them.  A 101-heavy-atom vancomycin
round-trips as well.

Document idempotence: doc -> mol -> doc is byte-stable for the 144 core molecules
in both parity modes (canonical Kekule form, §6 I8).

## 10. Limitations

* **RDKit has no native planar/helical stereo** — the extension stores those
  units on the mol property `_molecode_stereo_extra`; the document round-trip and
  the local 3D sign derivation are supported, but RDKit itself cannot judge them
  (use CXSMILES/InChIKey as appropriate and treat planar/helical as
  representation-level records).
* **Conformation** is out of scope: this profile encodes configuration, not
  geometry.  A separate geometry/conformation layer (bond lengths, angles,
  dihedrals) is orthogonal and does not replace configuration.
* Absolute CIP descriptors (`_R/_S`) remain a derived cache; they are
  recomputed on read and validated, not treated as source truth.
* No LLM-task experiments are bundled with this profile (it is a deterministic
  representation layer; downstream model evaluation is separate).

## 11. References

* Daylight SMILES: neighbour-order dependence of `@/@@`.
* CXSMILES / OpenSMILES: enhanced stereo `|a:|`, `|&1:|`, `|o1:|` (atom and bond members).
* MDL MOL V3000: `CFG` and `STERAC / STEREL / STEABS` collections.
* InChI Technical FAQ: stereochemistry is treated as absolute; enhanced stereo is
  dropped; axial chirality is not covered.
* Hanson et al., *J. Chem. Inf. Model.* 2018, 58, 1755-1765 (CIP implementation
  and its revision; the `hanson2018` cache refers to this).
* ConfSeq (bioRxiv 2025 / Nat. Mach. Intell. 2026): 3D internal-coordinate
  language; its pseudo-chirality CW/CCW is neighbour-order dependent and only
  applies to non-stereogenic centres, so it is complementary, not a substitute.
