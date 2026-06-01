#!/usr/bin/env python3
"""Example 8 — OCSR: molecule image -> MoleCode graph.

Optical Chemical Structure Recognition (OCSR) with a vision LLM: given a picture
of a molecule, recover its structure as a MoleCode graph. The Markush MoleCode
system prompt is used (the most general grammar — it also handles plain atoms
and aromatic rings), so this works for ordinary molecules and for structures
carrying abbreviated / R-group labels.

Pipeline:
  1. obtain a molecule image (either a path you pass on the command line, or a
     built-in demo molecule rendered to PNG with RDKit so the example is
     self-contained);
  2. send `image + MARKUSH_SYSTEM_PROMPT` to a vision-capable model;
  3. extract the ```mermaid block and parse it back to a structure with
     `molecode.markush.mermaid_to_mol` (validation);
  4. if the image came from a known molecule, compare the recovered SMILES.

Usage:
    # self-contained demo (renders aspirin, needs a vision model to recognise it)
    python examples/08_image_to_molecode.py
    # your own image
    python examples/08_image_to_molecode.py path/to/molecule.png

Needs a vision-capable model. Configure credentials (any OpenAI-compatible,
vision-capable endpoint):
    export MOLECODE_API_KEY="sk-..."
    export MOLECODE_MODEL="gpt-4o-mini"     # or gpt-4o, gemini-*, claude vision, ...
With no API key set it runs as a "dry run": it still renders/locates the image
and prints exactly what would be sent.
"""

import sys

from molecode.prompts import MARKUSH_SYSTEM_PROMPT
from molecode.markush import mermaid_to_mol, mol_to_smiles, has_invalid_atoms
from molecode.llm import LLMClient

DEMO_NAME = "Aspirin"
DEMO_SMILES = "CC(=O)Oc1ccccc1C(=O)O"

USER_INSTRUCTION = (
    "Read the molecular structure in this image and output it as a MoleCode "
    "graph. Respond with ONLY a single fenced ```mermaid code block."
)


def render_demo_image(path: str) -> None:
    """Render the demo molecule to a PNG so the example is self-contained."""
    from rdkit import Chem
    from rdkit.Chem import Draw
    mol = Chem.MolFromSmiles(DEMO_SMILES)
    Draw.MolToFile(mol, path, size=(450, 450))


def extract_mermaid(text: str) -> str:
    """Pull the fenced ```mermaid ... ``` block out of an LLM reply."""
    if "```" not in text:
        return text.strip()
    block = text.split("```", 2)[1]
    if block.startswith("mermaid"):
        block = block[len("mermaid"):]
    return block.strip()


def main() -> None:
    # 1. image: CLI arg, or render the built-in demo molecule.
    expected_smiles = None
    if len(sys.argv) > 1:
        image_path = sys.argv[1]
        print(f"Using image: {image_path}")
    else:
        image_path = "molecode_ocsr_demo.png"
        render_demo_image(image_path)
        expected_smiles = DEMO_SMILES
        print(f"Rendered demo molecule {DEMO_NAME} ({DEMO_SMILES}) -> {image_path}")

    # 2. send image + Markush MoleCode prompt to a vision model.
    try:
        client = LLMClient()
    except ValueError:
        print("\n" + "=" * 70)
        print("DRY RUN — no MOLECODE_API_KEY/OPENAI_API_KEY set.")
        print(f"Would send the image '{image_path}' plus this system prompt to a")
        print("vision-capable model and parse the returned MoleCode graph:\n")
        print(f"  system = MARKUSH_SYSTEM_PROMPT  ({len(MARKUSH_SYSTEM_PROMPT)} chars)")
        print(f"  user   = {USER_INSTRUCTION}")
        print(f"  images = ['{image_path}']")
        print("=" * 70)
        return

    print(f"Calling vision model '{client.model}' for OCSR ...")
    reply = client.chat(USER_INSTRUCTION, system=MARKUSH_SYSTEM_PROMPT,
                        images=[image_path])

    # 3. extract + validate the MoleCode graph.
    graph = extract_mermaid(reply)
    print("\nPredicted MoleCode graph:\n")
    print(graph)

    mol = mermaid_to_mol(graph, strict=False)
    if mol is None:
        print("\nCould not parse the predicted graph.")
        return
    recovered = mol_to_smiles(mol)
    print("\nRecovered SMILES:", recovered,
          "(contains abbreviation placeholders)" if has_invalid_atoms(mol) else "")

    # 4. compare against ground truth when we rendered a known molecule.
    if expected_smiles is not None:
        from rdkit import Chem
        try:
            match = Chem.CanonSmiles(recovered) == Chem.CanonSmiles(expected_smiles)
        except Exception:
            match = False
        print(f"Matches ground truth ({expected_smiles}): {match}")


if __name__ == "__main__":
    main()
