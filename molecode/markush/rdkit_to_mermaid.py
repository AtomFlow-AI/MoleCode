"""
将RDKit Mol对象转换为Mermaid Graph格式

主要功能：
1. 从Mol对象提取原子和键信息
2. 处理立体化学（E/Z构型）
3. 生成Mermaid图语法
4. 支持自定义命名和子图组织
"""

from typing import List, Tuple, Dict, Optional
from rdkit import Chem
from rdkit.Chem import Descriptors
from collections import defaultdict


class MolToMermaidConverter:
    """将RDKit Mol对象转换为Mermaid Graph格式"""

    # 键类型映射：RDKit → Mermaid
    BOND_TYPE_MAP = {
        Chem.BondType.SINGLE: '---',
        Chem.BondType.DOUBLE: '===',
        Chem.BondType.TRIPLE: '-.-',
        Chem.BondType.AROMATIC: '<-->',  #非必要不使用，尽量使用凯酷勒式
        Chem.BondType.DATIVE: '-->',
    }

    def __init__(self, subgraph_name: str = "Molecule"):
        """
        初始化转换器

        Args:
            subgraph_name: 子图名称（用于Mermaid subgraph声明）
        """
        self.subgraph_name = subgraph_name
        self.atom_id_map: Dict[int, str] = {}  # atom_idx -> atom_id
        self.element_counter: Dict[str, int] = defaultdict(int)  # 元素计数器

    def _sanitize_identifier(self, name: str) -> str:
        """
        清理字符串使其成为合法的Mermaid标识符

        移除特殊字符，只保留字母、数字、中文
        例如: "(E)-2-丁烯" → "E2丁烯"

        Args:
            name: 原始名称

        Returns:
            合法的标识符（只包含字母、数字、中文）
        """
        import re
        # 只保留字母、数字、中文，移除其他所有字符
        sanitized = re.sub(r'[^a-zA-Z0-9\u4e00-\u9fff]', '', name)

        # 确保不以数字开头
        if sanitized and sanitized[0].isdigit():
            sanitized = 'M' + sanitized

        # 如果清理后为空，使用默认名称
        if not sanitized:
            sanitized = 'Molecule'

        return sanitized

    def convert(self, mol: Chem.Mol) -> str:
        """
        将RDKit Mol转换为Mermaid Graph文本

        Args:
            mol: RDKit Mol对象

        Returns:
            Mermaid格式的图文本
        """
        if mol is None:
            return ""

        # 创建副本避免修改原始分子
        mol = Chem.Mol(mol)

        # 确保立体化学信息被正确分配（包括手性）
        try:
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        except Exception:
            pass

        # 将芳香分子转换为凯库勒式（显式单双键）
        # Graph只能表示显式键，不能表示芳香性
        try:
            Chem.Kekulize(mol, clearAromaticFlags=True)
        except Exception:
            # 如果凯库勒化失败（某些奇怪的结构），继续处理
            pass

        # 重置计数器
        self.atom_id_map = {}
        self.element_counter = defaultdict(int)

        # 1. 提取原子信息
        atoms_info = self._extract_atoms(mol)

        # 2. 提取键信息
        bonds_info = self._extract_bonds(mol)

        # 3. 生成Mermaid文本
        return self._generate_mermaid(atoms_info, bonds_info)

    def _extract_atoms(self, mol: Chem.Mol) -> List[Tuple[str, str, bool]]:
        """
        提取原子信息

        Args:
            mol: RDKit Mol对象

        Returns:
            [(atom_id, label, is_abbreviation), ...] 列表
        """
        atoms = []

        for atom in mol.GetAtoms():
            idx = atom.GetIdx()

            # 生成原子ID（唯一标识符，包含手性信息）
            atom_id = self._generate_atom_id(atom)
            self.atom_id_map[idx] = atom_id

            # Check if this is an abbreviation (Dummy Atom with _abbreviation prop)
            if atom.GetAtomicNum() == 0 and atom.HasProp("_abbreviation"):
                label = atom.GetProp("_abbreviation")
                atoms.append((atom_id, label, True))
            else:
                # 生成显示标签
                label = self._generate_atom_label(atom)
                atoms.append((atom_id, label, False))

        return atoms

    def _generate_atom_id(self, atom: Chem.Atom) -> str:
        """
        生成原子ID

        Args:
            atom: RDKit Atom对象

        Returns:
            原子ID，格式：SubgraphName_Symbol_Number 或 SubgraphName_Symbol_Number_Chirality
        """
        symbol = atom.GetSymbol()

        # Dummy Atom (*) -> use 'X' for valid identifier
        if symbol == '*' or atom.GetAtomicNum() == 0:
            symbol = 'X'

        # 按元素类型计数
        self.element_counter[symbol] += 1
        count = self.element_counter[symbol]

        # 清理subgraph名称以生成合法的ID
        clean_name = self._sanitize_identifier(self.subgraph_name)

        # 基础ID
        base_id = f"{clean_name}_{symbol}_{count}"

        # 使用 RDKit 计算出的绝对 CIP 构型，而不是直接把
        # CHI_TETRAHEDRAL_CW/CCW 当作 R/S。CW/CCW 依赖原子顺序，
        # 只有 _CIPCode 才是可序列化的绝对 R/S 标签。
        if atom.HasProp('_CIPCode'):
            cip_code = atom.GetProp('_CIPCode')
            if cip_code in ('R', 'S'):
                return f"{base_id}_{cip_code}"

        # 无手性或未指定
        return base_id

    def _generate_atom_label(self, atom: Chem.Atom) -> str:
        """
        生成原子显示标签

        Args:
            atom: RDKit Atom对象

        Returns:
            显示标签，如 'C', 'OH', 'NH2', 'N(+)', 'O(-)', 'O(2-)'
        """
        symbol = atom.GetSymbol()
        total_h = atom.GetTotalNumHs()
        formal_charge = atom.GetFormalCharge()

        # 构建标签
        label = symbol

        # 添加氢
        if total_h > 0:
            if total_h == 1:
                label += 'H'
            else:
                label += f'H{total_h}'

        # 添加电荷（新格式：括号包裹，数字在前，符号在后）
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
        提取键信息

        Args:
            mol: RDKit Mol对象

        Returns:
            [(atom1_id, atom2_id, bond_symbol), ...] 列表
        """
        bonds = []

        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()

            # 获取原子ID
            atom1_id = self.atom_id_map[begin_idx]
            atom2_id = self.atom_id_map[end_idx]

            # 确定键符号（包括立体化学）
            bond_symbol = self._get_bond_symbol(bond)

            bonds.append((atom1_id, atom2_id, bond_symbol))

        return bonds

    def _get_bond_symbol(self, bond: Chem.Bond) -> str:
        """
        根据键类型和立体化学生成Mermaid键符号

        Args:
            bond: RDKit Bond对象

        Returns:
            Mermaid键符号，如 '---', '===', '===|E|'
        """
        bond_type = bond.GetBondType()
        stereo = bond.GetStereo()

        # 获取基础键符号
        base_symbol = self.BOND_TYPE_MAP.get(bond_type, '---')

        # 如果是双键且有立体化学信息
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

    def _generate_mermaid(self, atoms: List[Tuple[str, str, bool]],
                         bonds: List[Tuple[str, str, str]]) -> str:
        """
        生成Mermaid Graph文本

        Args:
            atoms: 原子信息列表 [(atom_id, label, is_abbreviation), ...]
            bonds: 键信息列表 [(atom1_id, atom2_id, bond_symbol), ...]

        Returns:
            Mermaid格式的文本
        """
        lines = []

        # 图声明
        lines.append("graph TB")

        # 添加注释显示原始分子名称
        lines.append(f'    %% 原始分子名称: {self.subgraph_name}')

        # 生成合法的subgraph ID
        subgraph_id = self._sanitize_identifier(self.subgraph_name)

        # 子图开始
        lines.append(f'    subgraph {subgraph_id}["{self.subgraph_name}"]')

        # 添加原子定义
        for atom_id, label, is_abbrev in atoms:
            if is_abbrev:
                # 缩写基团用花括号
                lines.append(f'        {atom_id}{{{label}}}')
            else:
                # 普通原子用方括号
                lines.append(f'        {atom_id}[{label}]')

        # 空行分隔
        lines.append('')

        # 添加键连接
        for atom1_id, atom2_id, bond_symbol in bonds:
            lines.append(f'        {atom1_id} {bond_symbol} {atom2_id}')

        # 子图结束
        lines.append('    end')

        return '\n'.join(lines)


def mol_to_mermaid(mol: Chem.Mol, name: str = "Molecule") -> str:
    """
    将RDKit Mol转换为Mermaid Graph（便捷函数）

    Args:
        mol: RDKit Mol对象
        name: 分子名称（用作subgraph名称）

    Returns:
        Mermaid格式的文本


    """
    converter = MolToMermaidConverter(subgraph_name=name)
    return converter.convert(mol)


if __name__ == "__main__":
    import os

    print("=" * 60)
    print("RDKit Mol -> Mermaid Graph 转换测试")
    print("=" * 60)

    # 创建输出目录
    output_dir = "mol_graphs"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"创建输出目录: {output_dir}")

    # 测试用例
    test_cases = [
        ("Aspirin", "Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        # ("苯", "Benzene", "c1ccccc1"),
        # ("丙酮", "Acetone", "CC(=O)C"),
        # ("乙炔", "Acetylene", "C#C"),
        # ("(E)-2-丁烯", "E-2-Butene", r"C/C=C/C"),
        # ("(Z)-2-丁烯", "Z-2-Butene", r"C/C=C\C"),
        # ("乙酸", "AceticAcid", "CC(=O)O"),
        # ("甲胺", "Methylamine", "CN"),
        # ("毒素", "Palytoxin", "CC1CC2(C(OC(C1)(O2)CCCCCCCC(CC3C(C(C(C(O3)(CC(C(C)C=CC(CCC(C(C4CC(C(C(O4)CC(C(CC5C(C(C(C(O5)CC(C=CC=CCC(C(C(CC=CC(=C)CCC(C(C(C(C)CC6C(C(C(C(O6)C=CC(C(CC7CC8CC(O7)C(O8)CCC9C(CC(O9)CN)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)CC(C)CCCCCC(C(C(C(C(C1C(C(C(C(O1)CC(C(C(=CC(CC(C)C(C(=O)NC=CC(=O)NCCCO)O)O)C)O)O)O)O)O)O)O)O)O)O)C"),
        # ("Mol1","Example","O=C1C(O)=C(c2c(C3=COc4c(C3)cccc4)c(O)cc(O)c2)Oc5cc(O)cc(O)c51")
    ]

    # 索引文件内容
    index_content = ["# 分子图库\n", "## 目录\n\n"]

    for i, (name_cn, name_en, smiles) in enumerate(test_cases, 1):
        print(f"\n[{i}/{len(test_cases)}] 处理: {name_cn} ({name_en})")

        # 从SMILES创建Mol
        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            print(f"  ✗ 无法解析SMILES: {smiles}")
            continue

        # 转换为Mermaid
        mermaid_text = mol_to_mermaid(mol, name_cn)

        # 获取分子式
        from rdkit.Chem import rdMolDescriptors
        formula = rdMolDescriptors.CalcMolFormula(mol)

        # 生成MD文件内容
        md_content = f"""# {name_cn} ({name_en})

## 基本信息

- **SMILES**: `{smiles}`
- **分子式**: {formula}

## 分子结构图

```mermaid
{mermaid_text}
```

## 说明

此图由RDKit自动生成，展示了分子的拓扑结构和键类型。
"""

        # 写入文件
        filename = f"{name_en}.md"
        filepath = os.path.join(output_dir, filename)

        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(md_content)

        print(f"  ✓ 已生成: {filepath}")

        # 添加到索引
        index_content.append(f"{i}. [{name_cn} ({name_en})](./{output_dir}/{filename}) - `{smiles}`\n")

    # 生成索引文件
    index_content.append("\n---\n\n*由 RDKit 自动生成*\n")
    index_path = os.path.join(output_dir, "INDEX.md")

    with open(index_path, 'w', encoding='utf-8') as f:
        f.writelines(index_content)

    print(f"\n✓ 索引文件已生成: {index_path}")

    print("\n" + "=" * 60)
    print(f"所有测试完成！共生成 {len(test_cases)} 个文件")
    print(f"输出目录: {os.path.abspath(output_dir)}")
    print("=" * 60)
