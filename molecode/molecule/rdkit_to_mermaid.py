"""
Convert RDKit Mol objects to Mermaid Graph format

Key features:
1. Extract atom and bond information from a Mol object
2. Handle stereochemistry (E/Z configuration)
3. Generate Mermaid graph syntax
4. Support custom naming and subgraph organization
"""

from typing import List, Tuple, Dict, Optional
from rdkit import Chem
from rdkit.Chem import Descriptors
from collections import defaultdict


class MolToMermaidConverter:
    """Convert an RDKit Mol object to Mermaid Graph format."""

    # Bond type mapping: RDKit -> Mermaid
    BOND_TYPE_MAP = {
        Chem.BondType.SINGLE: '---',
        Chem.BondType.DOUBLE: '===',
        Chem.BondType.TRIPLE: '-.-',
        Chem.BondType.AROMATIC: '<-->',  # Use Kekule form whenever possible
        Chem.BondType.DATIVE: '-->',
    }

    def __init__(self, subgraph_name: str = "Molecule"):
        """
        Initialize the converter.

        Args:
            subgraph_name: Subgraph name (used in the Mermaid subgraph declaration)
        """
        self.subgraph_name = subgraph_name
        self.atom_id_map: Dict[int, str] = {}  # atom_idx -> atom_id
        self.element_counter: Dict[str, int] = defaultdict(int)  # Element counter

    def _sanitize_identifier(self, name: str) -> str:
        """
        Clean a string so it becomes a valid Mermaid identifier.

        Removes special characters, keeping only letters and digits.
        Example: "(E)-2-Butene" -> "E2Butene"

        Args:
            name: Original name

        Returns:
            A valid identifier (letters and digits only)
        """
        import re
        # Keep only letters and digits; remove everything else
        sanitized = re.sub(r'[^a-zA-Z0-9]', '', name)

        # Ensure the identifier does not start with a digit
        if sanitized and sanitized[0].isdigit():
            sanitized = 'M' + sanitized

        # If the result is empty after cleaning, use a default name
        if not sanitized:
            sanitized = 'Molecule'

        return sanitized

    def convert(self, mol: Chem.Mol, kekulize: bool = True) -> str:
        """
        Convert an RDKit Mol to Mermaid Graph text.

        Args:
            mol: RDKit Mol object
            kekulize: True -> Kekule form (explicit single/double bonds, default);
                      False -> retain aromatic bonds as <--> symbols (suitable for editing tasks)

        Returns:
            Graph text in Mermaid format
        """
        if mol is None:
            return ""

        # Work on a copy to avoid modifying the original molecule
        mol = Chem.Mol(mol)

        # Ensure stereochemistry (including chirality) is correctly assigned
        try:
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        except Exception:
            pass

        # Convert aromatic molecules to Kekule form (explicit single/double bonds)
        # Graph notation can only represent explicit bonds, not aromaticity
        if kekulize:
            try:
                Chem.Kekulize(mol, clearAromaticFlags=True)
            except Exception:
                # If Kekulization fails (unusual structures), continue anyway
                pass

        # Reset counters
        self.atom_id_map = {}
        self.element_counter = defaultdict(int)

        # 1. Extract atom information
        atoms_info = self._extract_atoms(mol)

        # 2. Extract bond information
        bonds_info = self._extract_bonds(mol)

        # 3. Generate Mermaid text
        return self._generate_mermaid(atoms_info, bonds_info)

    def _extract_atoms(self, mol: Chem.Mol) -> List[Tuple[str, str]]:
        """
        Extract atom information.

        Args:
            mol: RDKit Mol object

        Returns:
            List of [(atom_id, label), ...]
        """
        atoms = []

        for atom in mol.GetAtoms():
            idx = atom.GetIdx()

            # Generate atom ID (unique identifier, includes chirality information)
            atom_id = self._generate_atom_id(atom)
            self.atom_id_map[idx] = atom_id

            # Generate display label
            label = self._generate_atom_label(atom)

            atoms.append((atom_id, label))

        return atoms

    def _generate_atom_id(self, atom: Chem.Atom) -> str:
        """
        Generate an atom ID.

        Args:
            atom: RDKit Atom object

        Returns:
            Atom ID in the format: SubgraphName_Symbol_Number
            or SubgraphName_Symbol_Number_Chirality
        """
        symbol = atom.GetSymbol()

        # Count by element type
        self.element_counter[symbol] += 1
        count = self.element_counter[symbol]

        # Clean the subgraph name to produce a valid ID
        clean_name = self._sanitize_identifier(self.subgraph_name)

        # Base ID
        base_id = f"{clean_name}_{symbol}_{count}"

        # Use the absolute CIP configuration computed by RDKit rather than
        # treating CHI_TETRAHEDRAL_CW/CCW directly as R/S. CW/CCW depends on
        # atom ordering; only _CIPCode is a serializable absolute R/S label.
        if atom.HasProp('_CIPCode'):
            cip_code = atom.GetProp('_CIPCode')
            if cip_code in ('R', 'S'):
                return f"{base_id}_{cip_code}"

        # No chirality or unspecified
        return base_id

    def _generate_atom_label(self, atom: Chem.Atom) -> str:
        """
        Generate an atom display label.

        Args:
            atom: RDKit Atom object

        Returns:
            Display label such as 'C', 'OH', 'NH2', 'N(+)', 'O(-)', 'O(2-)'
        """
        symbol = atom.GetSymbol()
        total_h = atom.GetTotalNumHs()
        formal_charge = atom.GetFormalCharge()

        # Build label
        label = symbol

        # Append hydrogens
        if total_h > 0:
            if total_h == 1:
                label += 'H'
            else:
                label += f'H{total_h}'

        # Append charge (new format: wrapped in parentheses, number first, sign last)
        if formal_charge != 0:
            if formal_charge == 1:
                charge_str = '(+)'
            elif formal_charge == -1:
                charge_str = '(-)'
            elif formal_charge > 0:
                charge_str = f'({formal_charge}+)'
            else:
                charge_str = f'({-formal_charge}-)'
            label += charge_str

        return label

    def _extract_bonds(self, mol: Chem.Mol) -> List[Tuple[str, str, str]]:
        """
        Extract bond information.

        Args:
            mol: RDKit Mol object

        Returns:
            List of [(atom1_id, atom2_id, bond_symbol), ...]
        """
        bonds = []

        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()

            # Get atom IDs
            atom1_id = self.atom_id_map[begin_idx]
            atom2_id = self.atom_id_map[end_idx]

            # Determine the bond symbol (including stereochemistry)
            bond_symbol = self._get_bond_symbol(bond)

            bonds.append((atom1_id, atom2_id, bond_symbol))

        return bonds

    def _get_bond_symbol(self, bond: Chem.Bond) -> str:
        """
        Generate the Mermaid bond symbol based on bond type and stereochemistry.

        Args:
            bond: RDKit Bond object

        Returns:
            Mermaid bond symbol, e.g. '---', '===', '===|E|'
        """
        bond_type = bond.GetBondType()
        stereo = bond.GetStereo()

        # Get base bond symbol
        base_symbol = self.BOND_TYPE_MAP.get(bond_type, '---')

        # If the bond is a double bond with stereochemistry information
        if bond_type == Chem.BondType.DOUBLE:
            if stereo == Chem.BondStereo.STEREOE:
                return '===|E|'
            elif stereo == Chem.BondStereo.STEREOZ:
                return '===|Z|'
            elif stereo == Chem.BondStereo.STEREOCIS:
                return '===|cis|'
            elif stereo == Chem.BondStereo.STEREOTRANS:
                return '===|trans|'

        return base_symbol

    def _generate_mermaid(self, atoms: List[Tuple[str, str]],
                         bonds: List[Tuple[str, str, str]]) -> str:
        """
        Generate Mermaid Graph text.

        Args:
            atoms: Atom information list [(atom_id, label), ...]
            bonds: Bond information list [(atom1_id, atom2_id, bond_symbol), ...]

        Returns:
            Text in Mermaid format
        """
        lines = []

        # Graph declaration
        lines.append("graph TB")

        # Add a comment showing the original molecule name
        lines.append(f'    %% Original molecule name: {self.subgraph_name}')

        # Generate a valid subgraph ID
        subgraph_id = self._sanitize_identifier(self.subgraph_name)

        # Subgraph start
        lines.append(f'    subgraph {subgraph_id}["{self.subgraph_name}"]')

        # Add atom definitions
        for atom_id, label in atoms:
            lines.append(f'        {atom_id}[{label}]')

        # Blank line separator
        lines.append('')

        # Add bond connections
        for atom1_id, atom2_id, bond_symbol in bonds:
            lines.append(f'        {atom1_id} {bond_symbol} {atom2_id}')

        # Subgraph end
        lines.append('    end')

        return '\n'.join(lines)


def mol_to_mermaid(mol: Chem.Mol, name: str = "Molecule", kekulize: bool = True) -> str:
    """
    Convert an RDKit Mol to Mermaid Graph (convenience function).

    Args:
        mol: RDKit Mol object
        name: Molecule name (used as the subgraph name)
        kekulize: True -> Kekule form (explicit single/double bonds, default);
                  False -> retain aromatic bonds as <--> symbols (suitable for editing tasks)

    Returns:
        Text in Mermaid format
    """
    converter = MolToMermaidConverter(subgraph_name=name)
    return converter.convert(mol, kekulize=kekulize)


if __name__ == "__main__":
    import os

    print("=" * 60)
    print("RDKit Mol -> Mermaid Graph conversion test")
    print("=" * 60)

    # Create output directory
    output_dir = "mol_graphs"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    # Test cases
    test_cases = [
        ("Aspirin", "Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        # ("Benzene", "Benzene", "c1ccccc1"),
        # ("Acetone", "Acetone", "CC(=O)C"),
        # ("Acetylene", "Acetylene", "C#C"),
        # ("(E)-2-Butene", "E-2-Butene", r"C/C=C/C"),
        # ("(Z)-2-Butene", "Z-2-Butene", r"C/C=C\C"),
        # ("Acetic acid", "AceticAcid", "CC(=O)O"),
        # ("Methylamine", "Methylamine", "CN"),
        # ("Palytoxin", "Palytoxin", "CC1CC2(C(OC(C1)(O2)CCCCCCCC(CC3C(C(C(C(O3)(CC(C(C)C=CC(CCC(C(C4CC(C(C(O4)CC(C(CC5C(C(C(C(O5)CC(C=CC=CCC(C(C(CC=CC(=C)CCC(C(C(C(C)CC6C(C(C(C(O6)C=CC(C(CC7CC8CC(O7)C(O8)CCC9C(CC(O9)CN)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)CC(C)CCCCCC(C(C(C(C(C1C(C(C(C(O1)CC(C(C(=CC(CC(C)C(C(=O)NC=CC(=O)NCCCO)O)O)C)O)O)O)O)O)O)O)O)O)O)C"),
        # ("Mol1", "Example", "O=C1C(O)=C(c2c(C3=COc4c(C3)cccc4)c(O)cc(O)c2)Oc5cc(O)cc(O)c51")
    ]

    # Index file content
    index_content = ["# Molecule Graph Library\n", "## Table of Contents\n\n"]

    for i, (name_en, name_id, smiles) in enumerate(test_cases, 1):
        print(f"\n[{i}/{len(test_cases)}] Processing: {name_en} ({name_id})")

        # Create Mol from SMILES
        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            print(f"  x Cannot parse SMILES: {smiles}")
            continue

        # Convert to Mermaid
        mermaid_text = mol_to_mermaid(mol, name_en)

        # Get molecular formula
        from rdkit.Chem import rdMolDescriptors
        formula = rdMolDescriptors.CalcMolFormula(mol)

        # Generate MD file content
        md_content = f"""# {name_en} ({name_id})

## Basic Information

- **SMILES**: `{smiles}`
- **Formula**: {formula}

## Molecular Structure Diagram

```mermaid
{mermaid_text}
```

## Notes

This diagram was automatically generated by RDKit and shows the molecular topology and bond types.
"""

        # Write file
        filename = f"{name_id}.md"
        filepath = os.path.join(output_dir, filename)

        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(md_content)

        print(f"  + Generated: {filepath}")

        # Add to index
        index_content.append(f"{i}. [{name_en} ({name_id})](./{output_dir}/{filename}) - `{smiles}`\n")

    # Generate index file
    index_content.append("\n---\n\n*Automatically generated by RDKit*\n")
    index_path = os.path.join(output_dir, "INDEX.md")

    with open(index_path, 'w', encoding='utf-8') as f:
        f.writelines(index_content)

    print(f"\n+ Index file generated: {index_path}")

    print("\n" + "=" * 60)
    print(f"All tests complete! Generated {len(test_cases)} file(s).")
    print(f"Output directory: {os.path.abspath(output_dir)}")
    print("=" * 60)
