---
name: molecode
description: >-
  Use for deterministic molecule understanding, graph-level editing, generation,
  and validation with MoleCode — an explicit Mermaid graph in which every atom
  and bond is a typed, named node/edge. Converting SMILES to MoleCode exposes the
  molecule as explicit atoms, hydrogen counts, bonds, bond orders, stereochemistry
  (R/S, E/Z), aromatic/ring systems, functional groups, and linkers, so an agent
  can locate substructures and plan edits far more reliably than from linear
  SMILES. Trigger when the task is to understand or edit a molecular structure,
  replace a functional group, add/delete/substitute atoms, close or fuse rings,
  change a linker or scaffold, count atoms/elements, check formula/connectivity/
  valence, design a molecule under constraints, or convert between SMILES and
  MoleCode. Also covers polymers (PSMILES repeat units) and Markush structures
  (variable R-groups). Prefer MoleCode graph inspection and editing before
  hand-writing SMILES; the bundled script provides the conversion and validation API.
---

# MoleCode

Use the bundled script first. Do not call an external service for conversion.

MoleCode is the preferred working representation when an LLM needs to understand
or edit a molecule. Converting SMILES to MoleCode exposes the molecule as an
explicit [Mermaid](https://mermaid.js.org/) graph: every atom is a named node
`prefix_Element_N[Label]` carrying its element, hydrogen count and charge, and
every bond is an explicit edge (`---` single, `===` double, `-.-` triple). This
graph form is usually strong enough for the model to locate substructures, reason
about aromatic systems, identify fused rings, follow linkers, and plan edits
directly. Reconstruct the whole molecular topology from the explicit nodes and
edges before making any edit.

Treat MoleCode as a **code-like molecular data structure**: inspect it, add
concise `%%` comments when useful, patch it with normal file-editing tools, and
validate it with the bundled script. Directly editing the graph is usually
simpler and safer than generating a new SMILES by hand — adding a methyl is one
new node plus one new edge, not a whole-string rewrite.

## MoleCode-first decision rule

For any molecule understanding or editing request, work at the explicit graph
level first:

1. Convert the input SMILES (or PSMILES) to MoleCode with the bundled script.
2. Read the nodes and edges to understand the structure and the requested change.
3. Edit the graph directly (add/remove nodes and edges; relabel hydrogen counts).
4. Validate and convert back to the requested representation with the script.

Use direct SMILES editing only when the change is trivial and the structure is
already unambiguous. Reach for specialized cheminformatics code only when graph
inspection plus the script are not enough (unresolved aromaticity choices,
ambiguous canonical ranking, subtle stereochemical perception, persistent valence
debugging, conformer/geometric reasoning).

## The six conversion forms

MoleCode spans three structural domains, each convertible in both directions:

| Domain | To MoleCode | From MoleCode |
| --- | --- | --- |
| Small molecule | `smiles-to-molecode` | `molecode-to-smiles` |
| Polymer (PSMILES, `×n`) | `psmiles-to-molecode` | `molecode-to-psmiles` |
| Markush (`{}` R-groups) | *(author by hand)* | `markush-to-smiles` |

Plus utilities: `validate` (formula / atom & ring counts / round-trip),
`compare` (markush-aware graph isomorphism), `doctor` (environment check), and
`image-to-molecode` (OCSR: read a molecule **image** with a vision model and
emit a MoleCode graph — needs a vision-capable model + API key).

## File-based workflow (preferred for nontrivial edits)

Multiline graphs are easier to inspect and patch on disk. When the molecule is
large, the change is broad, or you are not confident emitting a complete correct
graph in one pass, use the file workflow:

1. Write source MoleCode to disk: `smiles-to-molecode "<smiles>" --output source.mmd`
2. Edit `source.mmd` into `edited.mmd` with normal text edits. Keep node ids
   stable unless a rename is necessary.
3. Convert back: `molecode-to-smiles --input edited.mmd --output product.smi`
4. Round-trip to confirm: `smiles-to-molecode --input product.smi`

Use task-local names like `source.mmd`, `edited.mmd`, `product.smi`, or a scratch
dir such as `.molecode-work/`. Delete scratch files when the workflow is complete
and the user did not ask to keep them — but first record the final SMILES/graph
and the validation evidence in your answer. If subagents are available, a
supervising agent may delegate the file-based edit to a subagent with a concrete
output contract (edited graph path, final SMILES, validation output, list of
atom/bond edits) and review the result rather than redo the whole pass.

## Authoring guardrails

When writing or editing a MoleCode graph by hand, read
`references/molecode-syntax.md` first. For polymers, Markush, stereochemistry,
multi-subgraph molecules, or full format details, read
`references/molecode-syntax-full.md`. Minimum contract:

- Each atom node is `prefix_Element_N[Label]`; the `Label` states the element plus
  explicit hydrogen count and charge (`[CH3]`, `[OH]`, `[N(+)]`, `[O(2-)]`).
- **Hydrogen counts are explicit.** When you add a bond to an atom, decrement its
  H count: bonding onto a `[CH3]` makes it `[CH2]`.
- Every bond is an explicit edge between two existing node ids.
- Aromatic rings are written in Kekulé form (alternating `===` / `---`) for small
  molecules.
- Stereochemistry: `===|E|` / `===|Z|` on double bonds; `_R` / `_S` (absolute CIP)
  suffix on chiral atom ids.
- Markush abbreviations use curly braces `{Boc}`, `{R1}`, `{Ar}`; each `{}` is one
  chemically meaningful group (split `NHBoc` into `[NH] --- {Boc}`).
- Validate every edited graph with the script before reporting a result.

## Quick start

```bash
python scripts/molecode_convert.py doctor
python scripts/molecode_convert.py smiles-to-molecode "CCO" --name Ethanol
python scripts/molecode_convert.py molecode-to-smiles --input edited.mmd
python scripts/molecode_convert.py validate --input edited.mmd
```

The script auto-locates the repository's `molecode` package, so it runs from
anywhere inside the repo. It needs `rdkit` (and `networkx` for `compare`); run
`doctor` first if the environment is unknown, and see
`references/dependencies.md` for install options.

## Workflow

1. Run `doctor` before SMILES work if the environment is unknown.
2. For atom/element counts, formula checks, connectivity checks, or any edit that
   starts from SMILES, convert to MoleCode first instead of inspecting the string.
3. For structural edits (functional-group replacement, atom add/delete/substitute,
   ring closure/fusion, linker/scaffold change), edit the graph directly.
4. For larger molecules or broad changes, use the file-based workflow.
5. Validate edited graphs with `validate` (formula, counts, round-trip), then
   export the requested SMILES/PSMILES with the matching `*-to-*` command.
6. For polymers, work on the repeat unit (two `*`); for Markush, author the graph
   with `{}` abbreviation nodes and check candidates with `compare`.
7. Clean up scratch files when they are not requested outputs.

## Resources

- `scripts/molecode_convert.py` — stable CLI over the repo's `molecode` package.
- `references/molecode-syntax.md` — compact syntax for manual molecule edits.
- `references/molecode-syntax-full.md` — full syntax: polymers, Markush, stereo,
  multi-subgraph molecules and reactions, edge cases.
- `references/operations.md` — command examples and supported conversion paths.
- `references/dependencies.md` — dependency policy and install options.
