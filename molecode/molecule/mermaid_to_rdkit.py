"""
将Mermaid Graph格式转换为RDKit Mol对象

RDKit分子表示方式：
1. Mol对象 - 分子的主要容器
2. Atom对象 - 包含原子序号、元素符号、氢数、形式电荷等
3. Bond对象 - 包含起始/终止原子索引、键类型（SINGLE, DOUBLE, TRIPLE等）
4. 邻接表结构 - 通过原子索引连接

示例：
    乙醇 (CH3CH2OH)
    - 3个原子: C(idx=0), C(idx=1), O(idx=2)
    - 2个键: C-C (0-1, SINGLE), C-O (1-2, SINGLE)
    - 隐式H会自动计算
"""

import re
from typing import Dict, List, Tuple, Optional
from rdkit import Chem
from rdkit.Chem import AllChem, Draw


class MermaidMolParser:
    """解析Mermaid分子图并转换为RDKit Mol对象"""

    # 键类型映射
    BOND_TYPE_MAP = {
        '---': Chem.BondType.SINGLE,
        '===': Chem.BondType.DOUBLE,
        '-.-': Chem.BondType.TRIPLE,
        '-->': Chem.BondType.DATIVE,  # 配位键
    }

    def __init__(self):
        self.atoms: Dict[str, int] = {}  # atom_id -> atom_label映射
        self.bonds: List[Tuple] = []  # (atom1_id, atom2_id, bond_type) 或 (atom1_id, atom2_id, bond_type, stereo)
        self.chirality: Dict[str, str] = {}  # atom_id -> chirality ('R' or 'S') 映射

    def parse_mermaid_graph(self, mermaid_text: str) -> Optional[Chem.Mol]:
        """
        解析Mermaid图文本并生成RDKit Mol对象

        Args:
            mermaid_text: Mermaid格式的分子图文本

        Returns:
            RDKit Mol对象，解析失败返回None
        """
        self.atoms = {}
        self.bonds = []
        self.chirality = {}

        lines = mermaid_text.strip().split('\n')

        for line in lines:
            line = line.strip()

            # 跳过注释、空行、graph声明、subgraph声明
            if (line.startswith('%%') or
                line.startswith('graph ') or
                line.startswith('subgraph ') or
                line == 'end' or
                not line):
                continue

            # 解析原子定义和键连接
            self._parse_line(line)

        # 构建RDKit Mol对象
        return self._build_mol()

    def _parse_line(self, line: str):
        """解析单行，提取原子和键信息"""

        # 先尝试匹配带立体化学标签的双键: atom1 ===|E| atom2 或 atom1 ===|Z| atom2
        # 原子ID可能包含手性后缀 (_R 或 _S)
        stereo_bond_pattern = r'([\w_]+)\s*===\|([EZez]|cis|trans|CIS|TRANS)\|\s*([\w_]+)'
        stereo_match = re.search(stereo_bond_pattern, line)

        if stereo_match:
            atom1_id = stereo_match.group(1)
            stereo_type = stereo_match.group(2).upper()  # 转大写统一处理
            atom2_id = stereo_match.group(3)
            self.bonds.append((atom1_id, atom2_id, '===', stereo_type))
            return

        # 尝试匹配普通键连接: atom1 bond_type atom2
        # 原子ID可能包含手性后缀 (_R 或 _S)
        bond_pattern = r'([\w_]+)\s*(---|\===|-\.-|-->)\s*([\w_]+)'
        bond_match = re.search(bond_pattern, line)

        if bond_match:
            atom1_id = bond_match.group(1)
            bond_type = bond_match.group(2)
            atom2_id = bond_match.group(3)
            self.bonds.append((atom1_id, atom2_id, bond_type))
            return

        # 尝试匹配原子定义: AtomID[Label] 或 AtomID_R[Label] 或 AtomID_S[Label]
        atom_pattern = r'([\w_]+?)(?:_(R|S))?\[([^\]]+)\]'
        atom_match = re.search(atom_pattern, line)

        if atom_match:
            base_id = atom_match.group(1)
            chirality = atom_match.group(2)  # 'R', 'S', 或 None
            label = atom_match.group(3)

            # 构建完整的原子ID（包含手性后缀）
            if chirality:
                atom_id = f"{base_id}_{chirality}"
                self.chirality[atom_id] = chirality
            else:
                atom_id = base_id

            if atom_id not in self.atoms:
                self.atoms[atom_id] = label

    def _parse_atom_label(self, label: str) -> Tuple[str, int, int]:
        """
        解析原子标签，提取元素符号、显式氢数和电荷

        Args:
            label: 原子标签，如 'C', 'OH', 'NH2', 'N(+)', 'O(-)', 'O(2-)'

        Returns:
            (元素符号, 显式氢数, 形式电荷)
            无效标签返回 ('*', 0, 0) - 使用Dummy Atom标记
        """
        label = label.strip()

        # 匹配元素符号、氢数和电荷（新格式：括号包裹，数字在前，符号在后）
        # 例如: C, OH, NH2, N(+), NH2(+), O(-), O(2-)
        match = re.match(r'^([A-Z][a-z]?)(?:H(\d*))?(?:\((\d*[+-])\))?$', label)

        if match:
            element = match.group(1)
            h_count_str = match.group(2)
            charge_str = match.group(3)

            # 验证是否为合法元素符号
            try:
                # 尝试创建原子以验证元素符号合法性
                test_atom = Chem.Atom(element)
                atomic_num = test_atom.GetAtomicNum()

                # 检查是否为真实元素（原子序数>0）
                if atomic_num == 0 and element != '*':
                    # 不是合法元素，返回Dummy Atom
                    return '*', 0, 0

            except Exception:
                # 创建失败，返回Dummy Atom
                return '*', 0, 0

            # 解析氢数
            if h_count_str is None:
                # 没有H标记
                h_count = 0
            elif h_count_str == '':
                # 有H但没有数字，表示1个H
                h_count = 1
            else:
                # 明确指定H的数量
                h_count = int(h_count_str)

            # 解析电荷（新格式：数字在前，符号在后）
            charge = 0
            if charge_str:
                # charge_str 格式: "+", "-", "2+", "2-", "3+" 等
                if charge_str == '+':
                    charge = 1
                elif charge_str == '-':
                    charge = -1
                else:
                    # 提取数字和符号
                    sign = charge_str[-1]  # 最后一个字符是符号
                    number = charge_str[:-1]  # 前面的是数字
                    magnitude = int(number)
                    charge = magnitude if sign == '+' else -magnitude

            return element, h_count, charge

        # 如果无法解析，返回Dummy Atom标记
        # Dummy Atom: 原子序数为0，符号为'*'，不会与真实元素混淆
        return '*', 0, 0

    def _build_mol(self) -> Optional[Chem.Mol]:
        """根据解析的原子和键信息构建RDKit Mol对象"""

        if not self.atoms:
            return None

        # 创建可编辑的分子
        mol = Chem.RWMol()

        # 原子ID到索引的映射
        atom_id_to_idx = {}

        # 添加原子
        for atom_id, label in self.atoms.items():
            element, h_count, charge = self._parse_atom_label(label)

            # 创建原子
            atom = Chem.Atom(element)

            # 设置显式氢数（如果有）
            if h_count > 0:
                atom.SetNumExplicitHs(h_count)

            # 设置形式电荷（如果有）
            if charge != 0:
                atom.SetFormalCharge(charge)

            # 添加到分子中
            idx = mol.AddAtom(atom)
            atom_id_to_idx[atom_id] = idx

            # 设置手性（如果有）
            if atom_id in self.chirality:
                chirality_type = self.chirality[atom_id]
                atom_obj = mol.GetAtomWithIdx(idx)

                if chirality_type == 'R':
                    atom_obj.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                elif chirality_type == 'S':
                    atom_obj.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)

        # 添加键
        for bond_info in self.bonds:
            if len(bond_info) == 3:
                # 普通键: (atom1_id, atom2_id, bond_type_str)
                atom1_id, atom2_id, bond_type_str = bond_info
                stereo_type = None
            elif len(bond_info) == 4:
                # 带立体化学的键: (atom1_id, atom2_id, bond_type_str, stereo_type)
                atom1_id, atom2_id, bond_type_str, stereo_type = bond_info
            else:
                continue

            if atom1_id in atom_id_to_idx and atom2_id in atom_id_to_idx:
                idx1 = atom_id_to_idx[atom1_id]
                idx2 = atom_id_to_idx[atom2_id]
                bond_type = self.BOND_TYPE_MAP.get(bond_type_str, Chem.BondType.SINGLE)

                mol.AddBond(idx1, idx2, bond_type)

                # 设置立体化学（稍后统一处理，需要先添加所有键）
                if stereo_type:
                    # 记录需要设置立体化学的键
                    if not hasattr(mol, '_stereo_bonds'):
                        mol._stereo_bonds = []
                    mol._stereo_bonds.append((idx1, idx2, stereo_type))

        # 在转换为不可编辑的Mol之前，设置立体化学
        if hasattr(mol, '_stereo_bonds'):
            for idx1, idx2, stereo_type in mol._stereo_bonds:
                bond = mol.GetBondBetweenAtoms(idx1, idx2)

                # 获取双键两端原子的邻接原子（用于定义立体化学）
                atom1 = mol.GetAtomWithIdx(idx1)
                atom2 = mol.GetAtomWithIdx(idx2)

                # 找到idx1的邻居（除了idx2）
                neighbors1 = [n.GetIdx() for n in atom1.GetNeighbors() if n.GetIdx() != idx2]
                # 找到idx2的邻居（除了idx1）
                neighbors2 = [n.GetIdx() for n in atom2.GetNeighbors() if n.GetIdx() != idx1]

                # 如果两端都有邻居，设置立体化学
                if neighbors1 and neighbors2:
                    # 使用第一个邻居作为参考原子
                    bond.SetStereoAtoms(neighbors1[0], neighbors2[0])

                    if stereo_type == 'E':
                        bond.SetStereo(Chem.BondStereo.STEREOE)
                    elif stereo_type == 'Z':
                        bond.SetStereo(Chem.BondStereo.STEREOZ)
                    elif stereo_type == 'CIS':
                        bond.SetStereo(Chem.BondStereo.STEREOCIS)
                    elif stereo_type == 'TRANS':
                        bond.SetStereo(Chem.BondStereo.STEREOTRANS)

        # 转换为不可编辑的Mol对象
        mol = mol.GetMol()

        # 清理分子结构（包括芳香性感知）
        try:
            Chem.SanitizeMol(mol)
            # SanitizeMol会自动进行芳香性感知
            # 但我们也可以显式设置以确保正确处理
            Chem.SetAromaticity(mol)
        except Exception:
            # 如果清理失败，尝试不带芳香性的清理
            try:
                Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
            except Exception:
                # 完全失败，返回未清理的版本
                pass

        return mol


def has_invalid_atoms(mol: Chem.Mol) -> bool:
    """
    检查分子是否包含无效原子（Dummy Atom）

    Args:
        mol: RDKit Mol对象

    Returns:
        True表示包含无效原子
    """
    if mol is None:
        return False

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:  # 原子序数0表示Dummy Atom
            return True
    return False


def get_invalid_atom_labels(mermaid_text: str, mol: Chem.Mol) -> List[str]:
    """
    获取所有无效原子的原始标签

    Args:
        mermaid_text: 原始Mermaid文本
        mol: 解析后的Mol对象

    Returns:
        无效原子的标签列表
    """
    if mol is None:
        return []

    invalid_labels = []
    parser = MermaidMolParser()

    # 重新解析以获取原始标签
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
    将Mermaid图格式转换为RDKit Mol对象（便捷函数）

    Args:
        mermaid_text: Mermaid格式的分子图文本
        strict: 严格模式。如果为True，遇到无效原子时返回None；
                如果为False，使用Dummy Atom (*)标记无效原子并继续

    Returns:
        RDKit Mol对象，解析失败返回None

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

    # 检查是否包含无效原子
    if strict and has_invalid_atoms(mol):
        invalid_labels = get_invalid_atom_labels(mermaid_text, mol)
        print(f"警告: 检测到无效原子标签: {invalid_labels}")
        print(f"这些原子已被标记为Dummy Atom (符号='*', 原子序数=0)")
        print(f"提示: 使用 strict=False 参数可以允许Dummy Atom")
        return None

    return mol


def mol_to_smiles(mol: Chem.Mol) -> str:
    """将RDKit Mol转换为SMILES字符串"""
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol)


def mol_to_inchi(mol: Chem.Mol) -> str:
    """将RDKit Mol转换为InChI字符串"""
    if mol is None:
        return ""
    return Chem.MolToInchi(mol)


def visualize_mol(mol: Chem.Mol, title: str = "", save_path: str = None):
    """
    可视化单个分子结构

    Args:
        mol: RDKit Mol对象
        title: 图像标题
        save_path: 保存路径，如果为None则不保存
    """
    if mol is None:
        print(f"无法可视化 {title}: Mol对象为None")
        return

    img = Draw.MolToImage(mol, size=(300, 300))

    if save_path:
        img.save(save_path)
        print(f"已保存到: {save_path}")

    return img


def visualize_mols_grid(mols: List[Chem.Mol],
                        legends: List[str] = None,
                        mols_per_row: int = 4,
                        sub_img_size: Tuple[int, int] = (300, 300),
                        save_path: str = None):
    """
    以网格形式可视化多个分子

    Args:
        mols: RDKit Mol对象列表
        legends: 每个分子的标签列表
        mols_per_row: 每行显示的分子数
        sub_img_size: 每个分子图像的大小
        save_path: 保存路径，如果为None则不保存

    Returns:
        PIL Image对象
    """
    if not mols:
        print("没有分子可以可视化")
        return None

    # 过滤掉None值
    valid_data = [(mol, legends[i] if legends else f"Mol {i+1}")
                  for i, mol in enumerate(mols) if mol is not None]

    if not valid_data:
        print("所有分子都是None")
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
        print(f"网格图已保存到: {save_path}")

    return img


if __name__ == "__main__":
    # 测试示例
    print("=" * 60)
    print("Mermaid分子图 -> RDKit Mol 转换测试")
    print("=" * 60)

    # 定义测试分子
    test_molecules = [
        ("mol1", """
        graph TB
            %% 原始分子名称: Mol1
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
            %% 原始分子名称: Mol1 (通过醚键添加呋喃环)
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
        
            %% 新增呋喃环子结构
            subgraph 呋喃环["呋喃环"]
                呋喃环_C_1[CH]
                呋喃环_C_2[C]
                呋喃环_C_3[CH]
                呋喃环_C_4[CH]
                呋喃环_O_1[O]
        
                呋喃环_C_1 === 呋喃环_C_2
                呋喃环_C_2 --- 呋喃环_C_3
                呋喃环_C_3 === 呋喃环_C_4
                呋喃环_C_4 --- 呋喃环_O_1
                呋喃环_O_1 --- 呋喃环_C_1
            end
        
            %% 通过醚键连接
            Mol1_O_2 --- 呋喃环_C_2
                """),

        ("丙酮", """
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

        ("乙炔", """
        graph TB
            subgraph Acetylene["Acetylene"]
                Acetylene_C_1[C]
                Acetylene_C_2[C]

                Acetylene_C_1 -.- Acetylene_C_2
            end
        """),

        ("(E)-2-丁烯", """
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

        ("(Z)-2-丁烯", """
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

    # 解析所有分子
    mols = []
    legends = []

    for i, (name, graph) in enumerate(test_molecules, 1):
        print(f"\n=== 测试{i}: {name} ===")
        mol = mermaid_to_mol(graph)

        if mol:
            smiles = mol_to_smiles(mol)
            formula = Chem.rdMolDescriptors.CalcMolFormula(mol)

            print(f"SMILES:  {smiles}")
            print(f"分子式:  {formula}")

            mols.append(mol)
            legends.append(f"{name}\n{smiles}")
        else:
            print(f"解析失败！")
            mols.append(None)
            legends.append(f"{name}\n(解析失败)")

    # 生成网格图
    print("\n" + "=" * 60)
    print("生成分子结构可视化...")
    print("=" * 60)

    visualize_mols_grid(
        mols,
        legends=legends,
        mols_per_row=4,
        sub_img_size=(350, 350),
        save_path="test_molecules.png"
    )

    print("\n所有测试完成！")
