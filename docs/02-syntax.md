# MoleCode — Syntax reference (small molecules)

> The Mermaid grammar for atoms, bonds, stereochemistry, and subgraphs.

A MoleCode document is a Mermaid graph. The authoritative, LLM-readable
specification lives in [`molecode.prompts.MOLECULE_SYSTEM_PROMPT`](../molecode/prompts);
this page is the human summary.

## Document structure

```mermaid
graph TB
    subgraph Ethanol["Ethanol"]
        Ethanol_C_1[C]
        Ethanol_C_2[C]
        Ethanol_O_3[OH]
        Ethanol_C_1 --- Ethanol_C_2
        Ethanol_C_2 --- Ethanol_O_3
    end
```

- `graph TB` (top-bottom) or `graph LR` (left-right) opens the document.
- A `subgraph ID["display name"]` … `end` block scopes one structural object. A
  molecule may be a single subgraph or split into several chemically meaningful
  subgraphs; bonds can cross subgraph boundaries.
- `%%` begins a comment — ignored by the parser, used to attach reasoning and
  validation notes inline.

## Atom nodes

Format: **`prefix_Element_Number[DisplayLabel]`**

```
chlorophenol_C_1[C]      ethanol_O_2[OH]      Mol_N_3[NH2]
```

- `prefix` — a namespace (usually the molecule/subgraph name) that keeps IDs
  unique and stable across prompt → reasoning → output.
- `Element_Number` — element symbol plus a per-element counter, giving each atom
  a **persistent identifier** (`C_1`, `C_2`, …).
- `[DisplayLabel]` — element + explicit hydrogen count + formal charge:

| Label | Meaning | Label | Meaning |
| --- | --- | --- | --- |
| `[CH3]` | carbon, 3 H | `[N(+)]` | nitrogen, +1 |
| `[CH2]` | carbon, 2 H | `[O(-)]` | oxygen, −1 |
| `[CH]` | carbon, 1 H | `[NH3(+)]` | ammonium |
| `[C]` | carbon, 0 H | `[O(2-)]` | oxide, −2 |
| `[OH]`, `[NH2]` | … | … | … |

> Because hydrogen counts are **explicit**, a local edit must keep them
> consistent: grafting a bond onto a `[CH3]` carbon turns it into `[CH2]`
> (see [`examples/06_editing.py`](../examples/06_editing.py)).

## Bond (edge) operators

| Operator | Bond | RDKit |
| --- | --- | --- |
| `---` | single | `SINGLE` |
| `===` | double | `DOUBLE` |
| `-.-` | triple | `TRIPLE` |
| `-->` | dative / coordinate | `DATIVE` (metal complexes) |
| `<-->` | aromatic | `AROMATIC` |

```
Acetylene_C_1 -.- Acetylene_C_2        %% triple bond
```

## Aromatic rings — Kekulé form

Aromatic rings are written in **Kekulé form** (explicit alternating single/double
bonds) for small molecules, so the topology stays unambiguous:

```
chlorophenol_C_1 === chlorophenol_C_2
chlorophenol_C_2 --- chlorophenol_C_3
chlorophenol_C_3 === chlorophenol_C_4
...
```

## Stereochemistry

- **Double-bond E/Z** — annotate the double bond:

  ```
  EB_C2 ===|E| EB_C3        %% also ===|Z|
  ```

- **Chirality R/S** — suffix the atom ID:

  ```
  R2butanol_C_2_R[CH]       %% _R or _S
  ```

- **Pseudo-asymmetric r/s** — lowercase suffix:

  ```
  Pentitol_C_3_s[CH]        %% modern CIP pseudo-asymmetric descriptor
  ```

### Stereo extension (profile `stereo/1`)

The base grammar above is unchanged.  Documents that need relative (AND/OR)
configuration, axial chirality, isotopes or an explicit local-parity source of
truth add an **extension layer** that the stock parser ignores:

```
graph TB
    %% molecode-profile: stereo/1
    %% cip-impl: hanson2018 (rdkit.Chem.rdCIPLabeler)
    subgraph Thr["relative configuration defined"]
        Thr_C_1[CH]
        Thr_C_3_S[CH]
        ... atoms and bonds as usual ...
        %% ---- stereo layer (stereo/1) ----
        Thr_TET_1{tet=Thr_C_1; par=cw}
        Thr_SG_1{sg=AND:1; mem=Thr_C_1,Thr_C_3_S}
    end
```

- `{...}` object lines and `-.->` member edges are skipped by the base parser;
  `%%` comments are ignored as usual.
- AND/OR members carry **no** absolute `_R/_S` suffix (sound base projection).
- Axial: `Pfx_AX_1{ax=cw; axis=Pfx_C_7,Pfx_C_8; ref1=Pfx_C_6; ref2=Pfx_C_9}`
- Groups may also contain axis members: `Pfx_SG_1{sg=AND:1; axes=Pfx_AX_1}`
- E/Z reference atoms (only when the default is ambiguous): `Pfx_EZ_1{ez=Z; bond=Pfx_C_2,Pfx_C_3; ref1=Pfx_C_1; ref2=Pfx_Br_1}`
- Planar: `Pfx_PL_1{pl=+; plane=Pfx_C_1,Pfx_C_2,Pfx_C_3; ref1=Pfx_C_4; ref2=Pfx_C_5}`
- Helical: `Pfx_HL_1{hl=M; path=Pfx_C_1,Pfx_C_2,Pfx_C_3,Pfx_C_4}`
- Isotope: `Pfx_ISO_1{iso=13; atom=Pfx_C_1}`

See [07-stereo.md](07-stereo.md) for the full grammar, reference-frame definition,
invariants and API.  Code:

```python
from molecode.molecule import mol_to_mermaid, mermaid_to_mol
doc = mol_to_mermaid(Chem.MolFromSmiles("..."), parity_mode="all")   # opt-in; default is "group_only"
mol = mermaid_to_mol(doc)
```

## Round-tripping in code

```python
from rdkit import Chem
from molecode.molecule import mol_to_mermaid, mermaid_to_mol, mol_to_smiles

graph = mol_to_mermaid(Chem.MolFromSmiles("CCO"), name="Ethanol")
mol   = mermaid_to_mol(graph)          # strict=True rejects unknown atoms
mol_to_smiles(mol)                     # -> 'CCO'
```

Next: [03-polymers.md](03-polymers.md) · [04-markush.md](04-markush.md)
