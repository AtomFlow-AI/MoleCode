"""
polymer_to_mermaid.py

Converts polymers (repeat units) to Mermaid Graph format, using subgraph to represent
Repeat Units, supporting homopolymers and block copolymers.

Consistent with the style of rdkit_to_mermaid.py; can be integrated directly into
the ChemIQ evaluation framework.

Usage examples (run from the mgl_code/ directory):
    python3 polymer_to_mermaid.py                    # run the built-in demo
    python3 polymer_to_mermaid.py --smiles "*CC*" --n 100   # custom input

Repeat unit SMILES convention:
    - Use * to mark the two attachment points, e.g. "*CC*" (polyethylene repeat unit)
    - The first * corresponds to the left connection, the second * to the right connection
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


# ── Bond type map (consistent with rdkit_to_mermaid.py) ──────────────────────

BOND_TYPE_MAP = {
    Chem.BondType.SINGLE:   "---",
    Chem.BondType.DOUBLE:   "===",
    Chem.BondType.TRIPLE:   "-.-",
    Chem.BondType.AROMATIC: "<-->",
}


# ── Data structures ────────────────────────────────────────────────────────────

@dataclass
class BlockSpec:
    """Describes one block in a copolymer."""
    smiles: str               # SMILES with two * attachment points, e.g. "*CC*"
    n: int                    # repeat count
    label: Optional[str] = None  # display label (None = auto-generated)


@dataclass
class PolymerSpec:
    """Complete polymer description."""
    blocks: List[BlockSpec]
    terminus_left:  Optional[str] = None  # left end-group SMILES (e.g. "C" for methyl)
    terminus_right: Optional[str] = None  # right end-group SMILES
    name: str = "Polymer"


# ── Atom/bond label utilities (counterpart to rdkit_to_mermaid._generate_atom_label) ──

def _sanitize_id(s: str) -> str:
    """Generate a valid Mermaid node ID (keep only alphanumeric characters)."""
    cleaned = re.sub(r"[^a-zA-Z0-9]", "", s)
    return ("M" + cleaned) if cleaned and cleaned[0].isdigit() else (cleaned or "M")


def _atom_label(atom: Chem.Atom) -> str:
    """Generate the atom display label, consistent with rdkit_to_mermaid._generate_atom_label logic."""
    symbol = atom.GetSymbol()
    total_h = atom.GetTotalNumHs()
    charge = atom.GetFormalCharge()

    label = symbol
    if total_h == 1:
        label += "H"
    elif total_h > 1:
        label += f"H{total_h}"
    if charge == 1:
        label += "(+)"
    elif charge == -1:
        label += "(-)"
    elif charge > 1:
        label += f"({charge}+)"
    elif charge < -1:
        label += f"({-charge}-)"
    return label


def _bond_symbol(bond: Chem.Bond) -> str:
    """Get the Mermaid symbol for a bond (including stereochemistry)."""
    bt = bond.GetBondType()
    stereo = bond.GetStereo()
    base = BOND_TYPE_MAP.get(bt, "---")
    if bt == Chem.BondType.DOUBLE:
        if stereo == Chem.BondStereo.STEREOE:
            return "===|E|"
        elif stereo == Chem.BondStereo.STEREOZ:
            return "===|Z|"
    return base


# ── Core converter class ───────────────────────────────────────────────────────

class RepeatUnitConverter:
    """
    Converts a single repeat-unit SMILES into a Mermaid subgraph fragment.

    Returns:
        - list of subgraph text lines
        - entry_node_id: node ID connecting to the left side
        - exit_node_id: node ID connecting to the right side
    """

    def __init__(self, block: BlockSpec, prefix: str, kekulize: bool = True):
        self.block = block
        self.prefix = prefix
        self.kekulize = kekulize
        self._element_counter: Dict[str, int] = defaultdict(int)
        self._node_map: Dict[int, str] = {}  # atom_idx -> node_id

    def _new_node_id(self, atom: Chem.Atom) -> str:
        symbol = atom.GetSymbol()
        self._element_counter[symbol] += 1
        cnt = self._element_counter[symbol]
        # Use the absolute CIP configuration computed by RDKit (_CIPCode) rather than
        # treating CHI_TETRAHEDRAL_CW/CCW directly as R/S -- CW/CCW depends on atom
        # ordering, while only _CIPCode is a serializable absolute R/S label.
        # AssignStereochemistry is called in convert() to ensure _CIPCode is present.
        suffix = ""
        if atom.HasProp("_CIPCode"):
            cip = atom.GetProp("_CIPCode")
            if cip in ("R", "S"):
                suffix = f"_{cip}"
        return f"{self.prefix}_{symbol}{cnt}{suffix}"

    def convert(self) -> Tuple[List[str], str, str]:
        """
        Returns:
            (lines, entry_node_id, exit_node_id)
        """
        mol = Chem.MolFromSmiles(self.block.smiles)
        if mol is None:
            raise ValueError(f"Cannot parse SMILES: {self.block.smiles!r}")

        # Kekulize (explicit single/double bonds) -- optional; skip for aromatic <--> notation
        if self.kekulize:
            try:
                Chem.Kekulize(mol, clearAromaticFlags=True)
            except Exception:
                pass

        try:
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        except Exception:
            pass

        # Find attachment points (* = atomic number 0)
        ap_idxs = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        if len(ap_idxs) != 2:
            raise ValueError(
                f"Repeat unit SMILES must contain exactly 2 * attachment points, "
                f"but found {len(ap_idxs)}: {self.block.smiles!r}"
            )

        ap_set = set(ap_idxs)

        # Find the real atom connected to each * (entry / exit)
        def _real_neighbor(ap_idx: int) -> int:
            for nbr in mol.GetAtomWithIdx(ap_idx).GetNeighbors():
                if nbr.GetAtomicNum() != 0:
                    return nbr.GetIdx()
            raise ValueError(f"Attachment point * (idx={ap_idx}) has no non-* neighbor atom")

        entry_real = _real_neighbor(ap_idxs[0])
        exit_real  = _real_neighbor(ap_idxs[1])

        # Assign node IDs to real atoms
        for atom in mol.GetAtoms():
            if atom.GetIdx() in ap_set:
                continue
            self._node_map[atom.GetIdx()] = self._new_node_id(atom)

        # Generate subgraph text
        sg_label = self.block.label or f"×{self.block.n}"
        lines: List[str] = [f'    subgraph {self.prefix}["{sg_label}"]']

        # Atom node definitions
        for atom in mol.GetAtoms():
            if atom.GetIdx() in ap_set:
                continue
            nid = self._node_map[atom.GetIdx()]
            lbl = _atom_label(atom)
            lines.append(f"        {nid}[{lbl}]")

        lines.append("")

        # Bond connections (skip bonds involving *)

        for bond in mol.GetBonds():
            bi, bj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if bi in ap_set or bj in ap_set:
                continue
            ni = self._node_map[bi]
            nj = self._node_map[bj]
            bs = _bond_symbol(bond)
            lines.append(f"        {ni} {bs} {nj}")

        lines.append("    end")

        return lines, self._node_map[entry_real], self._node_map[exit_real]


class PolymerToMermaidConverter:
    """Converts a PolymerSpec into a complete Mermaid Graph."""

    def __init__(self, kekulize: bool = True):
        self.kekulize = kekulize

    def convert(self, spec: PolymerSpec) -> str:
        lines: List[str] = ["graph LR"]
        lines.append(f"    %% Polymer: {spec.name}")
        lines.append("")

        # 1. Convert each block, collecting subgraph text and entry/exit node IDs
        block_data: List[Tuple[List[str], str, str]] = []
        for i, block in enumerate(spec.blocks):
            prefix = f"B{i}"
            conv = RepeatUnitConverter(block, prefix=prefix, kekulize=self.kekulize)
            sg_lines, entry_nid, exit_nid = conv.convert()
            block_data.append((sg_lines, entry_nid, exit_nid))

        # 2. Left end-group node
        tl_nid = "TL"
        tl_label = self._terminus_label(spec.terminus_left, default="H")
        lines.append(f"    {tl_nid}[{tl_label}]")

        # 3. Insert all subgraphs
        for sg_lines, _, _ in block_data:
            lines.extend(sg_lines)
            lines.append("")

        # 4. Right end-group node
        tr_nid = "TR"
        tr_label = self._terminus_label(spec.terminus_right, default="H")
        lines.append(f"    {tr_nid}[{tr_label}]")
        lines.append("")

        # 5. Connections: TL --- entry0, exit0 --- entry1, ..., exitN --- TR
        _, first_entry, _ = block_data[0]
        lines.append(f"    {tl_nid} --- {first_entry}")

        for i in range(len(block_data) - 1):
            _, _, exit_i      = block_data[i]
            _, entry_next, _  = block_data[i + 1]
            lines.append(f"    {exit_i} --- {entry_next}")

        _, _, last_exit = block_data[-1]
        lines.append(f"    {last_exit} --- {tr_nid}")

        return "\n".join(lines)

    @staticmethod
    def _terminus_label(smiles: Optional[str], default: str = "H") -> str:
        """Generate a display label from an end-group SMILES."""
        if smiles is None:
            return default
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() > 0:
                return _atom_label(atom)
        return default


# ── Convenience functions ─────────────────────────────────────────────────────

def polymer_to_mermaid(
    repeat_smiles: str,
    n: int,
    label: Optional[str] = None,
    terminus_left:  Optional[str] = None,
    terminus_right: Optional[str] = None,
    name: str = "Polymer",
    kekulize: bool = True,
) -> str:
    """
    Convenience interface for converting a homopolymer.

    Args:
        repeat_smiles:  Repeat unit SMILES with two *, e.g. "*CC*"
        n:              Repeat count
        label:          subgraph display label (default "xn")
        terminus_left:  Left end-group SMILES (e.g. "C" for methyl; None displays H)
        terminus_right: Right end-group SMILES
        name:           Graph annotation name
        kekulize:       True  -> Kekulize aromatic bonds to alternating single/double (default, for counting tasks)
                        False -> Keep aromatic bonds as <--> (for editing tasks)

    Returns:
        Mermaid format string
    """
    spec = PolymerSpec(
        blocks=[BlockSpec(smiles=repeat_smiles, n=n, label=label)],
        terminus_left=terminus_left,
        terminus_right=terminus_right,
        name=name,
    )
    return PolymerToMermaidConverter(kekulize=kekulize).convert(spec)


def block_copolymer_to_mermaid(
    blocks: List[BlockSpec],
    terminus_left:  Optional[str] = None,
    terminus_right: Optional[str] = None,
    name: str = "BlockCopolymer",
) -> str:
    """
    Convenience interface for converting a block copolymer.

    Args:
        blocks:         Ordered list of BlockSpec, arranged along the chain
        terminus_left:  Left end-group SMILES
        terminus_right: Right end-group SMILES
        name:           Graph annotation name

    Returns:
        Mermaid format string
    """
    spec = PolymerSpec(
        blocks=blocks,
        terminus_left=terminus_left,
        terminus_right=terminus_right,
        name=name,
    )
    return PolymerToMermaidConverter().convert(spec)


# ── Demo / CLI ─────────────────────────────────────────────────────────────────

DEMO_CASES = [
    {
        "name": "Polyethylene (PE)",
        "type": "homopolymer",
        "smiles": "*CC*",
        "n": 1000,
        "terminus_left": "C",
        "terminus_right": "C",
        "description": "Simplest linear homopolymer; equivalent to a straight-chain alkane for n<=100",
    },
    {
        "name": "Polypropylene (PP)",
        "type": "homopolymer",
        "smiles": "*CC(C)*",
        "n": 500,
        "terminus_left": "C",
        "terminus_right": "C",
        "description": "Homopolymer with methyl side chain",
    },
    {
        "name": "Polyethylene_glycol (PEG)",
        "type": "homopolymer",
        "smiles": "*CCO*",
        "n": 200,
        "terminus_left": "CO",
        "terminus_right": "O",
        "description": "Biomedical polymer; repeat unit contains oxygen",
    },
    {
        "name": "Polystyrene (PS)",
        "type": "homopolymer",
        "smiles": "*CC(c1ccccc1)*",
        "n": 300,
        "terminus_left": "C",
        "terminus_right": "C",
        "description": "Homopolymer with phenyl side chain (aromatic system)",
    },
    {
        "name": "Nylon-6,6",
        "type": "homopolymer",
        "smiles": "*C(=O)CCCCCC(=O)NCCCCCCN*",
        "n": 100,
        "terminus_left": None,
        "terminus_right": None,
        "description": "Amide-bond polymer with a relatively complex repeat unit",
    },
    {
        "name": "PEG-b-PPO block copolymer",
        "type": "block_copolymer",
        "blocks": [
            BlockSpec("*CCO*",     n=50,  label="PEG x50"),
            BlockSpec("*CC(C)O*",  n=30,  label="PPO x30"),
        ],
        "terminus_left": "CO",
        "terminus_right": "O",
        "description": "Diblock copolymer: hydrophilic PEG + hydrophobic PPO",
    },
]


def run_demo():
    print("=" * 70)
    print("polymer_to_mermaid  --  Polymer -> Mermaid Graph demo")
    print("=" * 70)

    for case in DEMO_CASES:
        print(f"\n{'─'*60}")
        print(f"[{case['name']}]")
        print(f"  Description: {case['description']}")

        if case["type"] == "homopolymer":
            mol = Chem.MolFromSmiles(case["smiles"].replace("*", "[*]"))
            formula_note = ""
            if mol:
                # For homopolymers, count the heavy atoms in the repeat unit
                ap_cnt = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 0)
                heavy_cnt = mol.GetNumAtoms() - ap_cnt
                formula_note = f"  Heavy atoms per unit: {heavy_cnt}  |  Total chain atoms ~= {heavy_cnt * case['n']}"

            print(f"  SMILES (repeat): {case['smiles']}  n={case['n']}")
            if formula_note:
                print(formula_note)

            graph = polymer_to_mermaid(
                case["smiles"],
                n=case["n"],
                terminus_left=case.get("terminus_left"),
                terminus_right=case.get("terminus_right"),
                name=case["name"],
            )
        else:
            blocks = case["blocks"]
            total_atoms = 0
            for b in blocks:
                mol = Chem.MolFromSmiles(b.smiles.replace("*", "[*]"))
                if mol:
                    total_atoms += (mol.GetNumAtoms() - 2) * b.n  # subtract 2 attachment points
            print(f"  Blocks: " + " | ".join(f"{b.label or b.smiles} n={b.n}" for b in blocks))
            print(f"  Total chain atoms ~= {total_atoms}")

            graph = block_copolymer_to_mermaid(
                blocks,
                terminus_left=case.get("terminus_left"),
                terminus_right=case.get("terminus_right"),
                name=case["name"],
            )

        print(f"\n```mermaid\n{graph}\n```")

    print("\n" + "=" * 70)
    print("Demo complete. Use polymer_to_mermaid() or block_copolymer_to_mermaid() to integrate into the evaluation pipeline.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Polymer repeat unit -> Mermaid Graph")
    parser.add_argument("--smiles",  type=str, help="Repeat unit SMILES (with two *), e.g. '*CC*'")
    parser.add_argument("--n",       type=int, help="Repeat count")
    parser.add_argument("--label",   type=str, default=None, help="subgraph display label")
    parser.add_argument("--tl",      type=str, default=None, help="Left end-group SMILES")
    parser.add_argument("--tr",      type=str, default=None, help="Right end-group SMILES")
    parser.add_argument("--demo",    action="store_true", help="Run the built-in demo (default behavior)")
    args = parser.parse_args()

    if args.smiles and args.n:
        result = polymer_to_mermaid(
            args.smiles, args.n,
            label=args.label,
            terminus_left=args.tl,
            terminus_right=args.tr,
        )
        print(result)
    else:
        run_demo()
