# AGENTS.md — MoleCode

Guidance for coding agents (Claude Code, Codex, and compatible tools) working in
this repository.

## What this repo is

**MoleCode** is an LLM-native molecular representation: it serializes a molecule
as an explicit [Mermaid](https://mermaid.js.org/) graph where every atom and bond
is a typed, named node/edge, losslessly inter-convertible with SMILES / MOL via
RDKit. It covers three domains with one grammar — small molecules, polymers
(PSMILES repeat units), and Markush structures (`{}` R-groups).

The installable library lives in [`molecode/`](molecode/); user-facing docs are
in [`docs/`](docs/); runnable examples in [`examples/`](examples/).

## Use the MoleCode skill

A ready-to-use **Agent Skill** ships in this repo at
[`.claude/skills/molecode/`](.claude/skills/molecode/):

- **Claude Code** auto-discovers it (standard `.claude/skills/` location). Just
  ask about understanding or editing a molecule and the `molecode` skill applies.
- **Codex / other agents**: read
  [`.claude/skills/molecode/SKILL.md`](.claude/skills/molecode/SKILL.md) and use
  the bundled CLI below. Interface metadata is in
  [`.claude/skills/molecode/agents/openai.yaml`](.claude/skills/molecode/agents/openai.yaml).

### MoleCode-first rule

For any molecule understanding or editing task, work at the explicit graph level:
convert SMILES → MoleCode, inspect/edit the named nodes and edges, then validate
and convert back. This is more reliable than reasoning over linear SMILES and
safer than hand-writing a new SMILES for structural edits.

### CLI quick start

```bash
# from the repo root (the script finds the molecode package automatically)
python .claude/skills/molecode/scripts/molecode_convert.py doctor
python .claude/skills/molecode/scripts/molecode_convert.py smiles-to-molecode "CCO" --name Ethanol
python .claude/skills/molecode/scripts/molecode_convert.py validate --input edited.mmd
python .claude/skills/molecode/scripts/molecode_convert.py molecode-to-smiles --input edited.mmd
```

The six conversion forms: `smiles-to-molecode` / `molecode-to-smiles`
(molecules), `psmiles-to-molecode` / `molecode-to-psmiles` (polymers),
`markush-to-smiles` and `compare` (Markush). Plus `validate` and `doctor`.

## Setup

```bash
pip install -e .          # installs rdkit + networkx; makes `molecode` importable
python examples/01_molecule_roundtrip.py                       # smoke check
python .claude/skills/molecode/scripts/molecode_convert.py doctor
```

If dependencies are missing, run the CLI `doctor` command and see
[`.claude/skills/molecode/references/dependencies.md`](.claude/skills/molecode/references/dependencies.md).
Do not install packages into the user's environment without asking.

## Conventions

- Hydrogen counts in MoleCode labels are **explicit**: when adding a bond to an
  atom, decrement its H count (`[CH3]` → `[CH2]`).
- Keep node ids stable when editing a graph unless a rename is necessary.
- Prefer the file-based workflow (`--input` / `--output`, scratch `.mmd` files)
  for larger molecules; clean up scratch files when done.
- Full grammar: [`.claude/skills/molecode/references/molecode-syntax-full.md`](.claude/skills/molecode/references/molecode-syntax-full.md).
