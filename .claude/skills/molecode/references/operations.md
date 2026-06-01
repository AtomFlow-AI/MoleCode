# Operations

Run `python scripts/molecode_convert.py <command> [args]` from anywhere inside
the MoleCode repository (the script finds the package automatically), or pass the
absolute script path. Every conversion command accepts `--input FILE` /
`--output FILE` for the file-based workflow.

## Small molecules

```bash
# SMILES -> MoleCode graph
python scripts/molecode_convert.py smiles-to-molecode "CCO" --name Ethanol
python scripts/molecode_convert.py smiles-to-molecode "c1ccccc1" --no-kekulize   # aromatic <-->

# MoleCode graph -> SMILES
python scripts/molecode_convert.py molecode-to-smiles --input edited.mmd
python scripts/molecode_convert.py molecode-to-smiles --input markushy.mmd --allow-invalid
```

`--no-kekulize` emits aromatic bonds as `<-->` instead of explicit alternating
single/double bonds. `--allow-invalid` keeps unrecognized/dummy atoms as `*`
rather than failing (useful for partial or Markush graphs).

## Polymers

A repeat unit is a PSMILES string with exactly two `*` attachment points.

```bash
# PSMILES -> MoleCode graph (repeat unit kept explicit, count carried as x n)
python scripts/molecode_convert.py psmiles-to-molecode "*NCCCCCC(=O)*" --n 8 --name Nylon-6

# MoleCode graph -> PSMILES
python scripts/molecode_convert.py molecode-to-psmiles --input polymer.mmd
```

R/S chirality and E/Z double-bond geometry round-trip:
`*C/C=C/C*` (trans) and `*C/C=C\C*` (cis) stay distinct.

## Markush structures

Markush graphs are authored by hand with `{}` abbreviation nodes (see
`molecode-syntax-full.md`). The script parses them and compares candidates.

```bash
# Markush MoleCode graph -> SMILES ({abbrev} become '*')
python scripts/molecode_convert.py markush-to-smiles --input scaffold.mmd

# Graph-isomorphism compare (markush-aware: {Me} == {CH3}, {Boc} == expanded form)
python scripts/molecode_convert.py compare --input candidate.mmd --ref reference.mmd
```

`compare` exits 0 when isomorphic, 1 otherwise, and prints the reason.

## OCSR (molecule image → MoleCode)

Recover a structure from a picture with a **vision-capable** model. Requires
`MOLECODE_API_KEY` and a vision model (`MOLECODE_MODEL=gpt-4o-mini`, `gpt-4o`,
a Gemini/Claude vision model, …). Uses the Markush MoleCode prompt, so it works
for ordinary molecules and for structures with abbreviated/R-group labels.

```bash
python scripts/molecode_convert.py image-to-molecode molecule.png
python scripts/molecode_convert.py image-to-molecode diagram.png --model gpt-4o --output out.mmd
```

It prints the predicted MoleCode graph (and the parsed SMILES to stderr). Pass an
image file path or an image URL.

## Validate

```bash
python scripts/molecode_convert.py validate --input edited.mmd
```

Reports `valid`, `smiles`, `formula`, `heavy_atoms`, `carbons`, `rings`, and
`round_trip_ok`. Use it after every manual edit to check formula and counts
against the requested constraints.

## File-based edit workflow

```bash
python scripts/molecode_convert.py smiles-to-molecode "OC12CC(Cl)C(F)(C1)C(Br)C2" --output source.mmd
# ... edit source.mmd -> edited.mmd with normal text edits (keep node ids stable) ...
python scripts/molecode_convert.py validate         --input edited.mmd
python scripts/molecode_convert.py molecode-to-smiles --input edited.mmd --output product.smi
python scripts/molecode_convert.py smiles-to-molecode --input product.smi --output roundtrip.mmd
```

Remove scratch `.mmd` / `.smi` files when they are not requested outputs, after
recording the final result and validation evidence.

## Environment check

```bash
python scripts/molecode_convert.py doctor
```

Reports Python version, located repo root, and `rdkit` / `networkx` / `molecode`
availability.
