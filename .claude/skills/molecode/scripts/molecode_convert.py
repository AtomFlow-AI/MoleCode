#!/usr/bin/env python3
"""Stable CLI for the MoleCode skill.

A thin, dependency-light wrapper around the repository's ``molecode`` package.
It exposes the six MoleCode conversion forms (molecule / polymer / markush, each
direction) plus ``validate`` and ``doctor`` utilities, with ``--input`` /
``--output`` file support for the file-based editing workflow that coding agents
should prefer on larger molecules.

Run from anywhere; the script locates the repository root automatically:

    python molecode_convert.py doctor
    python molecode_convert.py smiles-to-molecode "CCO" --name Ethanol
    python molecode_convert.py molecode-to-smiles --input edited.mmd
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


# ── Locate and import the repository's molecode package ───────────────────────

def _find_molecode_root() -> Path | None:
    """Walk up from this script to find the dir containing the molecode package."""
    here = Path(__file__).resolve()
    for parent in [here.parent, *here.parents]:
        if (parent / "molecode" / "__init__.py").is_file():
            return parent
    return None


_ROOT = _find_molecode_root()
if _ROOT is not None and str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))


RDKIT_HELP = """MoleCode needs RDKit (and networkx for markush comparison).

Least-invasive options (ask the user before installing):
  1. Use a Python env that already has rdkit + networkx.
  2. pip:    python -m pip install rdkit networkx
  3. conda:  conda install -c conda-forge rdkit networkx
  4. Install the repo itself:  python -m pip install -e .  (from the repo root)
"""


def _import_molecode():
    try:
        import molecode  # noqa: F401
        return molecode
    except ModuleNotFoundError as exc:
        if "molecode" in str(exc) and _ROOT is None:
            sys.stderr.write(
                "Could not locate the molecode package. Run this script from "
                "inside the MoleCode repository, or `pip install -e .` first.\n"
            )
        else:
            sys.stderr.write(f"{exc}\n\n{RDKIT_HELP}")
        raise SystemExit(2)


# ── I/O helpers ───────────────────────────────────────────────────────────────

def _read_payload(positional: str | None, input_path: str | None) -> str:
    if input_path:
        return Path(input_path).read_text(encoding="utf-8").strip()
    if positional is not None:
        return positional.strip()
    sys.stderr.write("Provide input as a positional argument or via --input FILE.\n")
    raise SystemExit(2)


def _emit(text: str, output_path: str | None) -> None:
    if output_path:
        Path(output_path).write_text(text + "\n", encoding="utf-8")
        print(output_path)
    else:
        print(text)


# ── Commands ──────────────────────────────────────────────────────────────────

def cmd_smiles_to_molecode(a) -> int:
    from rdkit import Chem
    from molecode.molecule import mol_to_mermaid
    smiles = _read_payload(a.smiles, a.input)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        sys.stderr.write(f"RDKit could not parse SMILES: {smiles!r}\n")
        return 1
    graph = mol_to_mermaid(mol, name=a.name, kekulize=not a.no_kekulize)
    _emit(graph, a.output)
    return 0


def cmd_molecode_to_smiles(a) -> int:
    from molecode.molecule import mermaid_to_mol, mol_to_smiles, has_invalid_atoms
    graph = _read_payload(None, a.input) if a.input else _read_payload(a.graph, None)
    mol = mermaid_to_mol(graph, strict=not a.allow_invalid)
    if mol is None:
        sys.stderr.write("Failed to parse MoleCode graph (invalid atoms or syntax).\n")
        return 1
    if has_invalid_atoms(mol) and not a.allow_invalid:
        sys.stderr.write("Graph contains unrecognized/dummy atoms; use --allow-invalid "
                         "to keep them as '*'.\n")
        return 1
    _emit(mol_to_smiles(mol), a.output)
    return 0


def cmd_psmiles_to_molecode(a) -> int:
    from molecode.polymer import polymer_to_mermaid
    psmiles = _read_payload(a.psmiles, a.input)
    graph = polymer_to_mermaid(psmiles, n=a.n, name=a.name,
                               kekulize=not a.no_kekulize)
    _emit(graph, a.output)
    return 0


def cmd_molecode_to_psmiles(a) -> int:
    from molecode.polymer import mermaid_to_psmiles
    graph = _read_payload(a.graph, a.input)
    psmiles = mermaid_to_psmiles(graph)
    if psmiles is None:
        sys.stderr.write("Failed to parse polymer MoleCode graph back to PSMILES.\n")
        return 1
    _emit(psmiles, a.output)
    return 0


def cmd_markush_to_smiles(a) -> int:
    from molecode.markush import mermaid_to_mol, mol_to_smiles
    graph = _read_payload(a.graph, a.input)
    mol = mermaid_to_mol(graph, strict=False)  # keep {abbrev} as '*' placeholders
    if mol is None:
        sys.stderr.write("Failed to parse markush MoleCode graph.\n")
        return 1
    _emit(mol_to_smiles(mol), a.output)
    return 0


def cmd_compare(a) -> int:
    """Graph-isomorphism comparison (markush-aware) of two MoleCode graphs."""
    from molecode.markush import MoleCodeGraph, molecode_isomorphic, EXPAND_MAP
    g1 = MoleCodeGraph.from_text(_read_payload(a.graph, a.input))
    g2 = MoleCodeGraph.from_text(_read_payload(None, a.ref))
    same, details = molecode_isomorphic(g1, g2, abbrev_expand_map=EXPAND_MAP)
    print(f"isomorphic: {same}")
    print(f"reason:     {details.get('reason')}")
    if details.get("unmatched_abbrevs"):
        print(f"unmatched_abbrevs: {details['unmatched_abbrevs']}")
    return 0 if same else 1


def cmd_image_to_molecode(a) -> int:
    """OCSR: read a molecule image with a vision LLM, output a MoleCode graph."""
    from molecode.prompts import MARKUSH_SYSTEM_PROMPT
    from molecode.markush import mermaid_to_mol, mol_to_smiles, has_invalid_atoms
    from molecode.llm import LLMClient

    try:
        client = LLMClient(model=a.model) if a.model else LLMClient()
    except ValueError as exc:
        sys.stderr.write(f"{exc}\n\nOCSR needs a vision-capable model. Set "
                         "MOLECODE_API_KEY (and MOLECODE_MODEL, e.g. gpt-4o-mini).\n")
        return 2

    instruction = ("Read the molecular structure in this image and output it as "
                   "a MoleCode graph. Respond with ONLY a single fenced ```mermaid "
                   "code block.")
    reply = client.chat(instruction, system=MARKUSH_SYSTEM_PROMPT, images=[a.image])

    graph = reply.strip()
    if "```" in graph:
        graph = graph.split("```", 2)[1]
        if graph.startswith("mermaid"):
            graph = graph[len("mermaid"):]
        graph = graph.strip()

    _emit(graph, a.output)
    if not a.output:
        mol = mermaid_to_mol(graph, strict=False)
        if mol is not None:
            note = " (with abbreviation placeholders)" if has_invalid_atoms(mol) else ""
            sys.stderr.write(f"# parsed SMILES: {mol_to_smiles(mol)}{note}\n")
    return 0


def cmd_validate(a) -> int:
    """Parse a molecule MoleCode graph and report structure facts + round-trip."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    from molecode.molecule import mermaid_to_mol, mol_to_smiles, has_invalid_atoms
    graph = _read_payload(a.graph, a.input)
    mol = mermaid_to_mol(graph, strict=False)
    if mol is None:
        print("valid: False")
        print("error: could not parse graph")
        return 1
    invalid = has_invalid_atoms(mol)
    smiles = mol_to_smiles(mol)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    n_atoms = mol.GetNumAtoms()
    n_c = sum(1 for at in mol.GetAtoms() if at.GetSymbol() == "C")
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    print(f"valid: {not invalid}")
    print(f"smiles: {smiles}")
    print(f"formula: {formula}")
    print(f"heavy_atoms: {n_atoms}")
    print(f"carbons: {n_c}")
    print(f"rings: {n_rings}")
    if invalid:
        print("note: graph contains unrecognized/dummy atoms (markush '*' or typo)")
    # round-trip check
    rt = mermaid_to_mol(graph, strict=False)
    print(f"round_trip_ok: {rt is not None and mol_to_smiles(rt) == smiles}")
    return 0 if not invalid else 1


def cmd_doctor(a) -> int:
    import importlib.util
    print(f"python: {sys.version.split()[0]}")
    print(f"script: {Path(__file__).resolve()}")
    print(f"repo_root: {_ROOT if _ROOT else 'NOT FOUND'}")
    have_rdkit = importlib.util.find_spec("rdkit") is not None
    have_nx = importlib.util.find_spec("networkx") is not None
    print(f"rdkit: {'available' if have_rdkit else 'missing'}")
    print(f"networkx: {'available' if have_nx else 'missing (markush compare needs it)'}")
    molecode_ok = False
    if _ROOT is not None:
        try:
            import molecode
            molecode_ok = True
            print(f"molecode: {molecode.__version__}")
        except Exception as exc:  # noqa: BLE001
            print(f"molecode: import failed ({exc})")
    else:
        print("molecode: package not found on path")
    if not (have_rdkit and have_nx and molecode_ok):
        print()
        print(RDKIT_HELP.strip())
    return 0


# ── Argument parsing ──────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="molecode_convert.py",
        description="Convert and validate molecules in the MoleCode (Mermaid) graph language.",
    )
    sub = p.add_subparsers(dest="command", required=True)

    def add_io(sp):
        sp.add_argument("--input", help="Read input from FILE instead of the argument.")
        sp.add_argument("--output", help="Write output to FILE instead of stdout.")

    s = sub.add_parser("smiles-to-molecode", help="Small-molecule SMILES -> MoleCode graph.")
    s.add_argument("smiles", nargs="?", help="SMILES string.")
    s.add_argument("--name", default="Molecule", help="Subgraph/molecule name.")
    s.add_argument("--no-kekulize", action="store_true",
                   help="Emit aromatic bonds as <--> instead of Kekule single/double.")
    add_io(s); s.set_defaults(func=cmd_smiles_to_molecode)

    s = sub.add_parser("molecode-to-smiles", help="MoleCode graph -> SMILES.")
    s.add_argument("graph", nargs="?", help="MoleCode graph text.")
    s.add_argument("--allow-invalid", action="store_true",
                   help="Keep unrecognized atoms as '*' instead of failing.")
    add_io(s); s.set_defaults(func=cmd_molecode_to_smiles)

    s = sub.add_parser("psmiles-to-molecode", help="Polymer PSMILES (two *) -> MoleCode graph.")
    s.add_argument("psmiles", nargs="?", help="Repeat-unit PSMILES with two * attachment points.")
    s.add_argument("--n", type=int, default=10, help="Repeat count for the x n label.")
    s.add_argument("--name", default="Polymer", help="Polymer name.")
    s.add_argument("--no-kekulize", action="store_true")
    add_io(s); s.set_defaults(func=cmd_psmiles_to_molecode)

    s = sub.add_parser("molecode-to-psmiles", help="Polymer MoleCode graph -> PSMILES.")
    s.add_argument("graph", nargs="?", help="Polymer MoleCode graph text.")
    add_io(s); s.set_defaults(func=cmd_molecode_to_psmiles)

    s = sub.add_parser("markush-to-smiles",
                       help="Markush MoleCode graph -> SMILES ({abbrev} become '*').")
    s.add_argument("graph", nargs="?", help="Markush MoleCode graph text.")
    add_io(s); s.set_defaults(func=cmd_markush_to_smiles)

    s = sub.add_parser("compare",
                       help="Graph-isomorphism compare two MoleCode graphs (markush-aware).")
    s.add_argument("graph", nargs="?", help="First graph (or use --input).")
    s.add_argument("--ref", required=True, help="Second graph FILE to compare against.")
    s.add_argument("--input", help="Read the first graph from FILE.")
    s.set_defaults(func=cmd_compare)

    s = sub.add_parser("image-to-molecode",
                       help="OCSR: molecule image -> MoleCode graph (needs a vision model).")
    s.add_argument("image", help="Path to a molecule image (PNG/JPG), or an image URL.")
    s.add_argument("--model", help="Override the (vision-capable) model, e.g. gpt-4o.")
    s.add_argument("--output", help="Write the MoleCode graph to FILE instead of stdout.")
    s.set_defaults(func=cmd_image_to_molecode)

    s = sub.add_parser("validate",
                       help="Parse a molecule graph; report SMILES/formula/counts/round-trip.")
    s.add_argument("graph", nargs="?", help="MoleCode graph text.")
    s.add_argument("--input", help="Read the graph from FILE.")
    s.set_defaults(func=cmd_validate)

    s = sub.add_parser("doctor", help="Report environment and dependency status.")
    s.set_defaults(func=cmd_doctor)

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    # doctor must run even if deps are missing; other commands import lazily.
    if args.command != "doctor":
        _import_molecode()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
