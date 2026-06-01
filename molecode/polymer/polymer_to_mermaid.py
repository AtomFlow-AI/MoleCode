"""
polymer_to_mermaid.py

将高分子（重复单元）转换为 Mermaid Graph 格式，使用 subgraph 表示 Repeat Unit，
支持均聚物和嵌段共聚物。

与 rdkit_to_mermaid.py 风格保持一致，可直接集成到 ChemIQ 评测体系。

用法示例（从 mgl_code/ 目录运行）：
    python3 polymer_to_mermaid.py                    # 运行内置演示
    python3 polymer_to_mermaid.py --smiles "*CC*" --n 100   # 自定义输入

重复单元 SMILES 约定：
    - 使用 * 标记两个连接点，例如 "*CC*" (聚乙烯 repeat unit)
    - 第一个 * 对应左侧连接，第二个 * 对应右侧连接
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


# ── 键类型映射（与 rdkit_to_mermaid.py 保持一致）──────────────────────────────

BOND_TYPE_MAP = {
    Chem.BondType.SINGLE:   "---",
    Chem.BondType.DOUBLE:   "===",
    Chem.BondType.TRIPLE:   "-.-",
    Chem.BondType.AROMATIC: "<-->",
}


# ── 数据结构 ───────────────────────────────────────────────────────────────────

@dataclass
class BlockSpec:
    """描述共聚物中一个嵌段。"""
    smiles: str               # 含两个 * 连接点的 SMILES，如 "*CC*"
    n: int                    # 重复次数
    label: Optional[str] = None  # 显示标签（None 则自动生成）


@dataclass
class PolymerSpec:
    """完整高分子描述。"""
    blocks: List[BlockSpec]
    terminus_left:  Optional[str] = None  # 左端基 SMILES（如 "C" 表示甲基）
    terminus_right: Optional[str] = None  # 右端基 SMILES
    name: str = "Polymer"


# ── 原子/键标签工具函数（与 rdkit_to_mermaid._generate_atom_label 对应）────────

def _sanitize_id(s: str) -> str:
    """生成合法的 Mermaid 节点 ID（只保留字母数字）。"""
    cleaned = re.sub(r"[^a-zA-Z0-9]", "", s)
    return ("M" + cleaned) if cleaned and cleaned[0].isdigit() else (cleaned or "M")


def _atom_label(atom: Chem.Atom) -> str:
    """生成原子显示标签，与 rdkit_to_mermaid._generate_atom_label 逻辑一致。"""
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
    """获取键的 Mermaid 符号（含立体化学）。"""
    bt = bond.GetBondType()
    stereo = bond.GetStereo()
    base = BOND_TYPE_MAP.get(bt, "---")
    if bt == Chem.BondType.DOUBLE:
        if stereo == Chem.BondStereo.STEREOE:
            return "===|E|"
        elif stereo == Chem.BondStereo.STEREOZ:
            return "===|Z|"
    return base


# ── 核心转换类 ─────────────────────────────────────────────────────────────────

class RepeatUnitConverter:
    """
    将单个重复单元 SMILES 转换为 Mermaid subgraph 片段。

    返回：
        - subgraph 文本行列表
        - entry_node_id：连接左侧的节点 ID
        - exit_node_id：连接右侧的节点 ID
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
        # 使用 RDKit 计算的绝对 CIP 构型（_CIPCode），而不是把
        # CHI_TETRAHEDRAL_CW/CCW 直接当 R/S —— CW/CCW 依赖原子顺序，
        # 只有 _CIPCode 才是可序列化的绝对 R/S 标签。
        # AssignStereochemistry 已在 convert() 中调用，确保 _CIPCode 存在。
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
            raise ValueError(f"无法解析 SMILES: {self.block.smiles!r}")

        # Kekulize（显式单/双键）— optional; skip for aromatic <--> notation
        if self.kekulize:
            try:
                Chem.Kekulize(mol, clearAromaticFlags=True)
            except Exception:
                pass

        try:
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        except Exception:
            pass

        # 找到连接点 (* = 原子序数 0)
        ap_idxs = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        if len(ap_idxs) != 2:
            raise ValueError(
                f"重复单元 SMILES 必须含恰好 2 个 * 连接点，"
                f"当前含 {len(ap_idxs)} 个：{self.block.smiles!r}"
            )

        ap_set = set(ap_idxs)

        # 找到与每个 * 相连的真实原子（entry / exit）
        def _real_neighbor(ap_idx: int) -> int:
            for nbr in mol.GetAtomWithIdx(ap_idx).GetNeighbors():
                if nbr.GetAtomicNum() != 0:
                    return nbr.GetIdx()
            raise ValueError(f"连接点 * (idx={ap_idx}) 没有非 * 邻居原子")

        entry_real = _real_neighbor(ap_idxs[0])
        exit_real  = _real_neighbor(ap_idxs[1])

        # 为真实原子分配节点 ID
        for atom in mol.GetAtoms():
            if atom.GetIdx() in ap_set:
                continue
            self._node_map[atom.GetIdx()] = self._new_node_id(atom)

        # 生成 subgraph 文本
        sg_label = self.block.label or f"×{self.block.n}"
        lines: List[str] = [f'    subgraph {self.prefix}["{sg_label}"]']

        # 原子节点定义
        for atom in mol.GetAtoms():
            if atom.GetIdx() in ap_set:
                continue
            nid = self._node_map[atom.GetIdx()]
            lbl = _atom_label(atom)
            lines.append(f"        {nid}[{lbl}]")

        lines.append("")

        # 键连接（跳过包含 * 的键）
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
    """将 PolymerSpec 转换为完整 Mermaid Graph。"""

    def __init__(self, kekulize: bool = True):
        self.kekulize = kekulize

    def convert(self, spec: PolymerSpec) -> str:
        lines: List[str] = ["graph LR"]
        lines.append(f"    %% 高分子: {spec.name}")
        lines.append("")

        # 1. 转换各嵌段，收集 subgraph 文本和 entry/exit 节点 ID
        block_data: List[Tuple[List[str], str, str]] = []
        for i, block in enumerate(spec.blocks):
            prefix = f"B{i}"
            conv = RepeatUnitConverter(block, prefix=prefix, kekulize=self.kekulize)
            sg_lines, entry_nid, exit_nid = conv.convert()
            block_data.append((sg_lines, entry_nid, exit_nid))

        # 2. 左端基节点
        tl_nid = "TL"
        tl_label = self._terminus_label(spec.terminus_left, default="H")
        lines.append(f"    {tl_nid}[{tl_label}]")

        # 3. 插入所有 subgraph
        for sg_lines, _, _ in block_data:
            lines.extend(sg_lines)
            lines.append("")

        # 4. 右端基节点
        tr_nid = "TR"
        tr_label = self._terminus_label(spec.terminus_right, default="H")
        lines.append(f"    {tr_nid}[{tr_label}]")
        lines.append("")

        # 5. 连接：TL --- entry0, exit0 --- entry1, ..., exitN --- TR
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
        """从端基 SMILES 生成显示标签。"""
        if smiles is None:
            return default
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() > 0:
                return _atom_label(atom)
        return default


# ── 便捷函数 ──────────────────────────────────────────────────────────────────

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
    均聚物转换便捷接口。

    Args:
        repeat_smiles:  含两个 * 的重复单元 SMILES，如 "*CC*"
        n:              重复次数
        label:          subgraph 显示标签（默认 "×n"）
        terminus_left:  左端基 SMILES（如 "C" 甲基，None 则显示 H）
        terminus_right: 右端基 SMILES
        name:           图注释名称
        kekulize:       True  → Kekulize芳香键为交替单/双键 (default, for counting tasks)
                        False → 保留芳香键为 <--> (for editing tasks)

    Returns:
        Mermaid 格式字符串
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
    嵌段共聚物转换便捷接口。

    Args:
        blocks:         有序 BlockSpec 列表，按链序排列
        terminus_left:  左端基 SMILES
        terminus_right: 右端基 SMILES
        name:           图注释名称

    Returns:
        Mermaid 格式字符串
    """
    spec = PolymerSpec(
        blocks=blocks,
        terminus_left=terminus_left,
        terminus_right=terminus_right,
        name=name,
    )
    return PolymerToMermaidConverter().convert(spec)


# ── 演示 / CLI ─────────────────────────────────────────────────────────────────

DEMO_CASES = [
    {
        "name": "聚乙烯 (PE)",
        "type": "homopolymer",
        "smiles": "*CC*",
        "n": 1000,
        "terminus_left": "C",
        "terminus_right": "C",
        "description": "最简单线性均聚物；C=100 以内的小分子等价于直链烷烃",
    },
    {
        "name": "聚丙烯 (PP)",
        "type": "homopolymer",
        "smiles": "*CC(C)*",
        "n": 500,
        "terminus_left": "C",
        "terminus_right": "C",
        "description": "含甲基侧链的均聚物",
    },
    {
        "name": "聚乙二醇 (PEG)",
        "type": "homopolymer",
        "smiles": "*CCO*",
        "n": 200,
        "terminus_left": "CO",
        "terminus_right": "O",
        "description": "生物医用高分子，重复单元含氧",
    },
    {
        "name": "聚苯乙烯 (PS)",
        "type": "homopolymer",
        "smiles": "*CC(c1ccccc1)*",
        "n": 300,
        "terminus_left": "C",
        "terminus_right": "C",
        "description": "苯环侧链均聚物（芳香体系）",
    },
    {
        "name": "尼龙-6,6 (Nylon-6,6)",
        "type": "homopolymer",
        "smiles": "*C(=O)CCCCCC(=O)NCCCCCCN*",
        "n": 100,
        "terminus_left": None,
        "terminus_right": None,
        "description": "酰胺键高分子，重复单元较复杂",
    },
    {
        "name": "PEG-b-PPO 嵌段共聚物",
        "type": "block_copolymer",
        "blocks": [
            BlockSpec("*CCO*",     n=50,  label="PEG ×50"),
            BlockSpec("*CC(C)O*",  n=30,  label="PPO ×30"),
        ],
        "terminus_left": "CO",
        "terminus_right": "O",
        "description": "两嵌段共聚物，亲水 PEG + 疏水 PPO",
    },
]


def run_demo():
    print("=" * 70)
    print("polymer_to_mermaid  —  高分子 → Mermaid Graph 演示")
    print("=" * 70)

    for case in DEMO_CASES:
        print(f"\n{'─'*60}")
        print(f"【{case['name']}】")
        print(f"  描述: {case['description']}")

        if case["type"] == "homopolymer":
            mol = Chem.MolFromSmiles(case["smiles"].replace("*", "[*]"))
            formula_note = ""
            if mol:
                # 对于均聚物，统计单元内的原子数
                ap_cnt = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 0)
                heavy_cnt = mol.GetNumAtoms() - ap_cnt
                formula_note = f"  单元重原子数: {heavy_cnt}  |  链总原子 ≈ {heavy_cnt * case['n']}"

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
            print(f"  嵌段: " + " | ".join(f"{b.label or b.smiles} n={b.n}" for b in blocks))
            print(f"  链总原子 ≈ {total_atoms}")

            graph = block_copolymer_to_mermaid(
                blocks,
                terminus_left=case.get("terminus_left"),
                terminus_right=case.get("terminus_right"),
                name=case["name"],
            )

        print(f"\n```mermaid\n{graph}\n```")

    print("\n" + "=" * 70)
    print("演示完成。使用 polymer_to_mermaid() 或 block_copolymer_to_mermaid() 集成到评测流程。")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="高分子 repeat unit → Mermaid Graph")
    parser.add_argument("--smiles",  type=str, help="重复单元 SMILES（含两个 *），如 '*CC*'")
    parser.add_argument("--n",       type=int, help="重复次数")
    parser.add_argument("--label",   type=str, default=None, help="subgraph 显示标签")
    parser.add_argument("--tl",      type=str, default=None, help="左端基 SMILES")
    parser.add_argument("--tr",      type=str, default=None, help="右端基 SMILES")
    parser.add_argument("--demo",    action="store_true", help="运行内置演示（默认行为）")
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
