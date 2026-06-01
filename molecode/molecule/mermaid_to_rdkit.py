"""
Convert Mermaid Graph format to RDKit Mol objects

RDKit molecule representation:
1. Mol object - primary container for molecules
2. Atom object - contains atomic number, element symbol, hydrogen count, formal charge, etc.
3. Bond object - contains begin/end atom indices, bond type (SINGLE, DOUBLE, TRIPLE, etc.)
4. Adjacency list structure - connected via atom indices

Example:
    Ethanol (CH3CH2OH)
    - 3 atoms: C(idx=0), C(idx=1), O(idx=2)
    - 2 bonds: C-C (0-1, SINGLE), C-O (1-2, SINGLE)
    - Implicit H is calculated automatically
"""

import re
from typing import Dict, List, Tuple, Optional
from rdkit import Chem
from rdkit.Chem import AllChem, Draw


class MermaidMolParser:
    """Parse a Mermaid molecule graph and convert it to an RDKit Mol object."""

    # Bond type mapping
    BOND_TYPE_MAP = {
        '---': Chem.BondType.SINGLE,
        '===': Chem.BondType.DOUBLE,
        '-.-': Chem.BondType.TRIPLE,
        '<-->': Chem.BondType.AROMATIC,
        '-->': Chem.BondType.DATIVE,  # Coordinate bond
    }

    def __init__(self):
        self.atoms: Dict[str, int] = {}  # atom_id -> atom_label mapping
        self.bonds: List[Tuple] = []  # (atom1_id, atom2_id, bond_type) or (atom1_id, atom2_id, bond_type, stereo)
        self.chirality: Dict[str, str] = {}  # atom_id -> chirality ('R' or 'S') mapping

    def parse_mermaid_graph(self, mermaid_text: str) -> Optional[Chem.Mol]:
        """
        Parse Mermaid graph text and produce an RDKit Mol object.

        Args:
            mermaid_text: Molecule graph text in Mermaid format

        Returns:
            RDKit Mol object, or None if parsing fails
        """
        self.atoms = {}
        self.bonds = []
        self.chirality = {}

        lines = mermaid_text.strip().split('\n')

        for line in lines:
            line = line.strip()

            # Skip comments, blank lines, graph declarations, subgraph declarations
            if (line.startswith('%%') or
                line.startswith('graph ') or
                line.startswith('subgraph ') or
                line == 'end' or
                not line):
                continue

            # Parse atom definitions and bond connections
            self._parse_line(line)

        # Build the RDKit Mol object
        return self._build_mol()

    def _parse_line(self, line: str):
        """Parse a single line to extract atom and bond information."""

        # First try to match a stereo-labeled double bond: atom1 ===|E| atom2 or atom1 ===|Z| atom2
        # Atom IDs may include a chirality suffix (_R or _S)
        stereo_bond_pattern = r'([\w_]+)\s*===\|([EZez]|cis|trans|CIS|TRANS)\|\s*([\w_]+)'
        stereo_match = re.search(stereo_bond_pattern, line)

        if stereo_match:
            atom1_id = stereo_match.group(1)
            stereo_type = stereo_match.group(2).upper()  # Normalize to uppercase
            atom2_id = stereo_match.group(3)
            self.bonds.append((atom1_id, atom2_id, '===', stereo_type))
            return

        # Try to match a plain bond connection: atom1 bond_type atom2
        # Atom IDs may include a chirality suffix (_R or _S)
        bond_pattern = r'([\w_]+)\s*(<-->|---|\===|-\.-|-->)\s*([\w_]+)'
        bond_match = re.search(bond_pattern, line)

        if bond_match:
            atom1_id = bond_match.group(1)
            bond_type = bond_match.group(2)
            atom2_id = bond_match.group(3)
            self.bonds.append((atom1_id, atom2_id, bond_type))
            return

        # Try to match an atom definition: AtomID[Label] or AtomID_R[Label] or AtomID_S[Label]
        atom_pattern = r'([\w_]+?)(?:_(R|S))?\[([^\]]+)\]'
        atom_match = re.search(atom_pattern, line)

        if atom_match:
            base_id = atom_match.group(1)
            chirality = atom_match.group(2)  # 'R', 'S', or None
            label = atom_match.group(3)

            # Build the full atom ID (including chirality suffix)
            if chirality:
                atom_id = f"{base_id}_{chirality}"
                self.chirality[atom_id] = chirality
            else:
                atom_id = base_id

            if atom_id not in self.atoms:
                self.atoms[atom_id] = label

    def _parse_atom_label(self, label: str) -> Tuple[str, int, int]:
        """
        Parse an atom label to extract element symbol, explicit hydrogen count, and charge.

        Args:
            label: Atom label, e.g. 'C', 'OH', 'NH2', 'N(+)', 'O(-)', 'O(2-)'

        Returns:
            (element symbol, explicit hydrogen count, formal charge)
            Returns ('*', 0, 0) for invalid labels - uses Dummy Atom marker
        """
        label = label.strip()

        # Match element symbol, hydrogen count, and charge (new format: wrapped in parentheses, number first, sign last)
        # Examples: C, OH, NH2, N(+), NH2(+), O(-), O(2-)
        match = re.match(r'^([A-Z][a-z]?)(?:H(\d*))?(?:\((\d*[+-])\))?$', label)

        if match:
            element = match.group(1)
            h_count_str = match.group(2)
            charge_str = match.group(3)

            # Validate element symbol
            try:
                # Attempt to create an atom to verify the element symbol is valid
                test_atom = Chem.Atom(element)
                atomic_num = test_atom.GetAtomicNum()

                # Check for a real element (atomic number > 0)
                if atomic_num == 0 and element != '*':
                    # Not a valid element; return Dummy Atom
                    return '*', 0, 0

            except Exception:
                # Creation failed; return Dummy Atom
                return '*', 0, 0

            # Parse hydrogen count
            if h_count_str is None:
                # No H marker
                h_count = 0
            elif h_count_str == '':
                # H present but no digit means 1 H
                h_count = 1
            else:
                # Explicit hydrogen count
                h_count = int(h_count_str)

            # Parse charge (new format: number first, sign last)
            charge = 0
            if charge_str:
                # charge_str format: "+", "-", "2+", "2-", "3+", etc.
                if charge_str == '+':
                    charge = 1
                elif charge_str == '-':
                    charge = -1
                else:
                    # Extract number and sign
                    sign = charge_str[-1]  # Last character is the sign
                    number = charge_str[:-1]  # Leading characters are the number
                    magnitude = int(number)
                    charge = magnitude if sign == '+' else -magnitude

            return element, h_count, charge

        # If parsing fails, return Dummy Atom marker
        # Dummy Atom: atomic number 0, symbol '*', will not be confused with real elements
        return '*', 0, 0

    def _build_mol(self) -> Optional[Chem.Mol]:
        """Build an RDKit Mol object from the parsed atom and bond information."""

        if not self.atoms:
            return None

        # Create an editable molecule
        mol = Chem.RWMol()

        # Mapping from atom ID to index
        atom_id_to_idx = {}

        # Add atoms
        for atom_id, label in self.atoms.items():
            element, h_count, charge = self._parse_atom_label(label)

            # Create atom
            atom = Chem.Atom(element)

            # Set explicit hydrogen count (if any)
            if h_count > 0:
                atom.SetNumExplicitHs(h_count)

            # Set formal charge (if any)
            if charge != 0:
                atom.SetFormalCharge(charge)

            # Add to molecule
            idx = mol.AddAtom(atom)
            atom_id_to_idx[atom_id] = idx

        # Add bonds. Double bond E/Z and tetrahedral R/S must be set after the full
        # topology is present, so record them here and restore after SanitizeMol.
        stereo_bonds = []
        aromatic_atom_idxs = set()
        for bond_info in self.bonds:
            if len(bond_info) == 3:
                # Plain bond: (atom1_id, atom2_id, bond_type_str)
                atom1_id, atom2_id, bond_type_str = bond_info
                stereo_type = None
            elif len(bond_info) == 4:
                # Stereo bond: (atom1_id, atom2_id, bond_type_str, stereo_type)
                atom1_id, atom2_id, bond_type_str, stereo_type = bond_info
            else:
                continue

            if atom1_id in atom_id_to_idx and atom2_id in atom_id_to_idx:
                idx1 = atom_id_to_idx[atom1_id]
                idx2 = atom_id_to_idx[atom2_id]
                bond_type = self.BOND_TYPE_MAP.get(bond_type_str, Chem.BondType.SINGLE)

                mol.AddBond(idx1, idx2, bond_type)

                if bond_type_str == '<-->':
                    aromatic_atom_idxs.update((idx1, idx2))
                    bond = mol.GetBondBetweenAtoms(idx1, idx2)
                    bond.SetIsAromatic(True)

                if stereo_type:
                    stereo_bonds.append((idx1, idx2, stereo_type))

        for idx in aromatic_atom_idxs:
            mol.GetAtomWithIdx(idx).SetIsAromatic(True)

        # Convert to an immutable Mol object
        mol = mol.GetMol()

        # Sanitize molecule structure (includes aromaticity perception)
        try:
            Chem.SanitizeMol(mol)
            # SanitizeMol performs aromaticity perception automatically;
            # set it explicitly here as well to ensure correct handling
            Chem.SetAromaticity(mol)
        except Exception:
            # If sanitization fails, try without aromaticity
            try:
                Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
            except Exception:
                # Complete failure; return the unsanitized version
                pass

        self._assign_chirality_from_ids(mol, atom_id_to_idx)
        self._assign_double_bond_stereo(mol, stereo_bonds)

        return mol


    def _assign_chirality_from_ids(self, mol: Chem.Mol, atom_id_to_idx: Dict[str, int]):
        """Restore absolute CIP chirality from the _R/_S suffix of atom IDs."""
        if not self.chirality:
            return

        for atom_id, desired_cip in self.chirality.items():
            idx = atom_id_to_idx.get(atom_id)
            if idx is None:
                continue

            atom = mol.GetAtomWithIdx(idx)
            matched = False

            for chiral_tag in (
                Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
            ):
                atom.SetChiralTag(chiral_tag)
                try:
                    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
                except Exception:
                    continue

                if atom.HasProp('_CIPCode') and atom.GetProp('_CIPCode') == desired_cip:
                    matched = True
                    break

            if not matched:
                atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

        try:
            Chem.AssignStereochemistry(mol, cleanIt=False, force=True)
        except Exception:
            pass

    def _assign_double_bond_stereo(self, mol: Chem.Mol, stereo_bonds: List[Tuple[int, int, str]]):
        """Restore ===|E| / ===|Z| double bond geometry."""
        for idx1, idx2, stereo_type in stereo_bonds:
            bond = mol.GetBondBetweenAtoms(idx1, idx2)
            if bond is None:
                continue

            # Get the neighboring atoms on each end of the double bond (used to define stereochemistry)
            atom1 = mol.GetAtomWithIdx(idx1)
            atom2 = mol.GetAtomWithIdx(idx2)

            neighbors1 = [n.GetIdx() for n in atom1.GetNeighbors() if n.GetIdx() != idx2]
            neighbors2 = [n.GetIdx() for n in atom2.GetNeighbors() if n.GetIdx() != idx1]

            if not neighbors1 or not neighbors2:
                continue

            bond.SetStereoAtoms(neighbors1[0], neighbors2[0])

            if stereo_type == 'E':
                bond.SetStereo(Chem.BondStereo.STEREOE)
            elif stereo_type == 'Z':
                bond.SetStereo(Chem.BondStereo.STEREOZ)
            elif stereo_type == 'CIS':
                bond.SetStereo(Chem.BondStereo.STEREOCIS)
            elif stereo_type == 'TRANS':
                bond.SetStereo(Chem.BondStereo.STEREOTRANS)

        try:
            # Do not use cleanIt=True to avoid clearing the E/Z labels
            # explicitly restored from the MoleCode.
            Chem.AssignStereochemistry(mol, cleanIt=False, force=True)
        except Exception:
            pass


def has_invalid_atoms(mol: Chem.Mol) -> bool:
    """
    Check whether a molecule contains invalid atoms (Dummy Atoms).

    Args:
        mol: RDKit Mol object

    Returns:
        True if the molecule contains invalid atoms
    """
    if mol is None:
        return False

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:  # Atomic number 0 indicates a Dummy Atom
            return True
    return False


def get_invalid_atom_labels(mermaid_text: str, mol: Chem.Mol) -> List[str]:
    """
    Get the original labels of all invalid atoms.

    Args:
        mermaid_text: Original Mermaid text
        mol: Parsed Mol object

    Returns:
        List of labels for invalid atoms
    """
    if mol is None:
        return []

    invalid_labels = []
    parser = MermaidMolParser()

    # Re-parse to retrieve original labels
    lines = mermaid_text.strip().split('\n')
    for line in lines:
        line = line.strip()
        atom_pattern = r'([\w_]+?)(?:_(R|S))?\[([^\]]+)\]'
        match = re.search(atom_pattern, line)

        if match:
            label = match.group(3)
            element, _, _ = parser._parse_atom_label(label)
            if element == '*':
                invalid_labels.append(label)

    return invalid_labels


def mermaid_to_mol(mermaid_text: str, strict: bool = True) -> Optional[Chem.Mol]:
    """
    Convert Mermaid graph format to an RDKit Mol object (convenience function).

    Args:
        mermaid_text: Molecule graph text in Mermaid format
        strict: Strict mode. If True, returns None when invalid atoms are found;
                if False, marks invalid atoms with Dummy Atom (*) and continues.

    Returns:
        RDKit Mol object, or None if parsing fails

    Example:
        >>> mermaid_graph = '''
        ... graph TB
        ...     subgraph Ethanol["Ethanol"]
        ...         Ethanol_C_1[C]
        ...         Ethanol_C_2[C]
        ...         Ethanol_O_3[OH]
        ...
        ...         Ethanol_C_1 --- Ethanol_C_2
        ...         Ethanol_C_2 --- Ethanol_O_3
        ...     end
        ... '''
        >>> mol = mermaid_to_mol(mermaid_graph)
        >>> print(Chem.MolToSmiles(mol))
        CCO
    """
    parser = MermaidMolParser()
    mol = parser.parse_mermaid_graph(mermaid_text)

    if mol is None:
        return None

    # Check for invalid atoms
    if strict and has_invalid_atoms(mol):
        invalid_labels = get_invalid_atom_labels(mermaid_text, mol)
        print(f"Warning: invalid atom labels detected: {invalid_labels}")
        print(f"These atoms have been marked as Dummy Atoms (symbol='*', atomic number=0).")
        print(f"Tip: pass strict=False to allow Dummy Atoms.")
        return None

    return mol


def mol_to_smiles(mol: Chem.Mol) -> str:
    """Convert an RDKit Mol to a SMILES string."""
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol)


def mol_to_inchi(mol: Chem.Mol) -> str:
    """Convert an RDKit Mol to an InChI string."""
    if mol is None:
        return ""
    return Chem.MolToInchi(mol)


def visualize_mol(mol: Chem.Mol, title: str = "", save_path: str = None):
    """
    Visualize a single molecule structure.

    Args:
        mol: RDKit Mol object
        title: Image title
        save_path: Save path; if None the image is not saved
    """
    if mol is None:
        print(f"Cannot visualize {title}: Mol object is None")
        return

    img = Draw.MolToImage(mol, size=(300, 300))

    if save_path:
        img.save(save_path)
        print(f"Saved to: {save_path}")

    return img


def visualize_mols_grid(mols: List[Chem.Mol],
                        legends: List[str] = None,
                        mols_per_row: int = 4,
                        sub_img_size: Tuple[int, int] = (300, 300),
                        save_path: str = None):
    """
    Visualize multiple molecules in a grid layout.

    Args:
        mols: List of RDKit Mol objects
        legends: List of labels for each molecule
        mols_per_row: Number of molecules per row
        sub_img_size: Size of each molecule image
        save_path: Save path; if None the image is not saved

    Returns:
        PIL Image object
    """
    if not mols:
        print("No molecules to visualize.")
        return None

    # Filter out None values
    valid_data = [(mol, legends[i] if legends else f"Mol {i+1}")
                  for i, mol in enumerate(mols) if mol is not None]

    if not valid_data:
        print("All molecules are None.")
        return None

    valid_mols, valid_legends = zip(*valid_data)

    img = Draw.MolsToGridImage(
        valid_mols,
        molsPerRow=mols_per_row,
        subImgSize=sub_img_size,
        legends=valid_legends
    )

    if save_path:
        img.save(save_path)
        print(f"Grid image saved to: {save_path}")

    return img


if __name__ == "__main__":
    # Test examples
    print("=" * 60)
    print("Mermaid molecule graph -> RDKit Mol conversion test")
    print("=" * 60)

    # Define test molecules
    test_molecules = [
        ("mol1", """
        graph TB
            %% Original molecule name: Mol1
            subgraph Mol1["Mol1"]
                Mol1_O_1[O]
                Mol1_C_1[C]
                Mol1_C_2[C]
                Mol1_O_2[OH]
                Mol1_C_3[C]
                Mol1_C_4[C]
                Mol1_C_5[C]
                Mol1_C_6[C]
                Mol1_C_7[CH]
                Mol1_O_3[O]
                Mol1_C_8[C]
                Mol1_C_9[C]
                Mol1_C_10[CH2]
                Mol1_C_11[CH]
                Mol1_C_12[CH]
                Mol1_C_13[CH]
                Mol1_C_14[CH]
                Mol1_C_15[C]
                Mol1_O_4[OH]
                Mol1_C_16[CH]
                Mol1_C_17[C]
                Mol1_O_5[OH]
                Mol1_C_18[CH]
                Mol1_O_6[O]
                Mol1_C_19[C]
                Mol1_C_20[CH]
                Mol1_C_21[C]
                Mol1_O_7[OH]
                Mol1_C_22[CH]
                Mol1_C_23[C]
                Mol1_O_8[OH]
                Mol1_C_24[C]

                Mol1_O_1 === Mol1_C_1
                Mol1_C_1 --- Mol1_C_2
                Mol1_C_2 --- Mol1_O_2
                Mol1_C_2 === Mol1_C_3
                Mol1_C_3 --- Mol1_C_4
                Mol1_C_4 --- Mol1_C_5
                Mol1_C_5 --- Mol1_C_6
                Mol1_C_6 === Mol1_C_7
                Mol1_C_7 --- Mol1_O_3
                Mol1_O_3 --- Mol1_C_8
                Mol1_C_8 === Mol1_C_9
                Mol1_C_9 --- Mol1_C_10
                Mol1_C_9 --- Mol1_C_11
                Mol1_C_11 === Mol1_C_12
                Mol1_C_12 --- Mol1_C_13
                Mol1_C_13 === Mol1_C_14
                Mol1_C_5 === Mol1_C_15
                Mol1_C_15 --- Mol1_O_4
                Mol1_C_15 --- Mol1_C_16
                Mol1_C_16 === Mol1_C_17
                Mol1_C_17 --- Mol1_O_5
                Mol1_C_17 --- Mol1_C_18
                Mol1_C_3 --- Mol1_O_6
                Mol1_O_6 --- Mol1_C_19
                Mol1_C_19 === Mol1_C_20
                Mol1_C_20 --- Mol1_C_21
                Mol1_C_21 --- Mol1_O_7
                Mol1_C_21 === Mol1_C_22
                Mol1_C_22 --- Mol1_C_23
                Mol1_C_23 --- Mol1_O_8
                Mol1_C_23 === Mol1_C_24
                Mol1_C_24 --- Mol1_C_1
                Mol1_C_18 === Mol1_C_4
                Mol1_C_10 --- Mol1_C_6
                Mol1_C_14 --- Mol1_C_8
                Mol1_C_24 --- Mol1_C_19
            end
        """),

        ("mol2", """
        graph TB
            %% Original molecule name: Mol1 (furan ring added via ether bond)
            subgraph Mol1["Mol1"]
                Mol1_O_1[O]
                Mol1_C_1[C]
                Mol1_C_2[C]
                Mol1_O_2[O]
                Mol1_C_3[C]
                Mol1_C_4[C]
                Mol1_C_5[C]
                Mol1_C_6[C]
                Mol1_C_7[CH]
                Mol1_O_3[O]
                Mol1_C_8[C]
                Mol1_C_9[C]
                Mol1_C_10[CH2]
                Mol1_C_11[CH]
                Mol1_C_12[CH]
                Mol1_C_13[CH]
                Mol1_C_14[CH]
                Mol1_C_15[C]
                Mol1_O_4[OH]
                Mol1_C_16[CH]
                Mol1_C_17[C]
                Mol1_O_5[OH]
                Mol1_C_18[CH]
                Mol1_O_6[O]
                Mol1_C_19[C]
                Mol1_C_20[CH]
                Mol1_C_21[C]
                Mol1_O_7[OH]
                Mol1_C_22[CH]
                Mol1_C_23[C]
                Mol1_O_8[OH]
                Mol1_C_24[C]

                Mol1_O_1 === Mol1_C_1
                Mol1_C_1 --- Mol1_C_2
                Mol1_C_2 --- Mol1_O_2
                Mol1_C_2 === Mol1_C_3
                Mol1_C_3 --- Mol1_C_4
                Mol1_C_4 --- Mol1_C_5
                Mol1_C_5 --- Mol1_C_6
                Mol1_C_6 === Mol1_C_7
                Mol1_C_7 --- Mol1_O_3
                Mol1_O_3 --- Mol1_C_8
                Mol1_C_8 === Mol1_C_9
                Mol1_C_9 --- Mol1_C_10
                Mol1_C_9 --- Mol1_C_11
                Mol1_C_11 === Mol1_C_12
                Mol1_C_12 --- Mol1_C_13
                Mol1_C_13 === Mol1_C_14
                Mol1_C_5 === Mol1_C_15
                Mol1_C_15 --- Mol1_O_4
                Mol1_C_15 --- Mol1_C_16
                Mol1_C_16 === Mol1_C_17
                Mol1_C_17 --- Mol1_O_5
                Mol1_C_17 --- Mol1_C_18
                Mol1_C_3 --- Mol1_O_6
                Mol1_O_6 --- Mol1_C_19
                Mol1_C_19 === Mol1_C_20
                Mol1_C_20 --- Mol1_C_21
                Mol1_C_21 --- Mol1_O_7
                Mol1_C_21 === Mol1_C_22
                Mol1_C_22 --- Mol1_C_23
                Mol1_C_23 --- Mol1_O_8
                Mol1_C_23 === Mol1_C_24
                Mol1_C_24 --- Mol1_C_1
                Mol1_C_18 === Mol1_C_4
                Mol1_C_10 --- Mol1_C_6
                Mol1_C_14 --- Mol1_C_8
                Mol1_C_24 --- Mol1_C_19
            end

            %% New furan ring substructure
            subgraph furan_ring["furan_ring"]
                furan_ring_C_1[CH]
                furan_ring_C_2[C]
                furan_ring_C_3[CH]
                furan_ring_C_4[CH]
                furan_ring_O_1[O]

                furan_ring_C_1 === furan_ring_C_2
                furan_ring_C_2 --- furan_ring_C_3
                furan_ring_C_3 === furan_ring_C_4
                furan_ring_C_4 --- furan_ring_O_1
                furan_ring_O_1 --- furan_ring_C_1
            end

            %% Connected via ether bond
            Mol1_O_2 --- furan_ring_C_2
                """),

        ("acetone", """
        graph TB
            subgraph Acetone["Acetone"]
                Acetone_C_1[C]
                Acetone_C_2[C]
                Acetone_C_3[C]
                Acetone_O_4[O]

                Acetone_C_1 --- Acetone_C_2
                Acetone_C_2 --- Acetone_C_3
                Acetone_C_2 === Acetone_O_4
            end
        """),

        ("acetylene", """
        graph TB
            subgraph Acetylene["Acetylene"]
                Acetylene_C_1[C]
                Acetylene_C_2[C]

                Acetylene_C_1 -.- Acetylene_C_2
            end
        """),

        ("(E)-2-butene", """
        graph TB
            subgraph EButene["(E)-2-Butene"]
                EB_C1[C]
                EB_C2[C]
                EB_C3[C]
                EB_C4[C]

                EB_C1 --- EB_C2
                EB_C2 ===|E| EB_C3
                EB_C3 --- EB_C4
            end
        """),

        ("(Z)-2-butene", """
        graph TB
            subgraph ZButene["(Z)-2-Butene"]
                ZB_C1[C]
                ZB_C2[C]
                ZB_C3[C]
                ZB_C4[C]

                ZB_C1 --- ZB_C2
                ZB_C2 ===|Z| ZB_C3
                ZB_C3 --- ZB_C4
            end
        """),
    ]

    # Parse all molecules
    mols = []
    legends = []

    for i, (name, graph) in enumerate(test_molecules, 1):
        print(f"\n=== Test {i}: {name} ===")
        mol = mermaid_to_mol(graph)

        if mol:
            smiles = mol_to_smiles(mol)
            formula = Chem.rdMolDescriptors.CalcMolFormula(mol)

            print(f"SMILES:  {smiles}")
            print(f"Formula: {formula}")

            mols.append(mol)
            legends.append(f"{name}\n{smiles}")
        else:
            print(f"Parsing failed!")
            mols.append(None)
            legends.append(f"{name}\n(parsing failed)")

    # Generate grid image
    print("\n" + "=" * 60)
    print("Generating molecule structure visualization...")
    print("=" * 60)

    visualize_mols_grid(
        mols,
        legends=legends,
        mols_per_row=4,
        sub_img_size=(350, 350),
        save_path="test_molecules.png"
    )

    print("\nAll tests complete!")
