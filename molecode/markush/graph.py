"""
MoleCode Graph -- Pure text-level MoleCode graph parsing and graph-isomorphism comparison (no RDKit dependency)

Core features:
1. Parse node + edge graph structure from MoleCode mermaid text
2. Graph-isomorphism comparison based on networkx (supports node/edge attribute matching)
3. Normalized comparison after abbreviation expansion
"""

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import networkx as nx


@dataclass
class NodeInfo:
    """Node information in an MoleCode graph."""
    node_id: str          # Raw ID (e.g. Mol_C_1)
    label: str            # Display label (e.g. CH3, Boc, R1)
    label_type: str       # "atom" (square brackets []) or "abbrev" (curly braces {})
    chirality: str = ""   # R, S, or empty


@dataclass
class EdgeInfo:
    """Edge information in an MoleCode graph."""
    src: str              # Source node ID
    dst: str              # Destination node ID
    bond_type: str        # ---, ===, -.- , -->
    stereo: str = ""      # E, Z, CIS, TRANS, or empty


class MoleCodeGraph:
    """Graph structure parsed from MoleCode mermaid text."""

    def __init__(self):
        self.nodes: Dict[str, NodeInfo] = {}
        self.edges: List[EdgeInfo] = []

    @classmethod
    def from_text(cls, text: str) -> 'MoleCodeGraph':
        """Parse graph structure from MoleCode mermaid text."""
        graph = cls()
        if not text:
            return graph

        for line in text.strip().split('\n'):
            line = line.strip()

            # Skip comments, blank lines, graph/subgraph/end
            if (not line or
                line.startswith('%%') or
                line.startswith('graph ') or
                line.startswith('subgraph ') or
                line == 'end'):
                continue

            graph._parse_line(line)

        return graph

    def _parse_line(self, line: str):
        """Parse a single MoleCode line."""

        # Double bond with stereochemistry: atom1 ===|E| atom2
        m = re.search(r'([\w_]+)\s*===\|([EZez]|cis|trans|CIS|TRANS)\|\s*([\w_]+)', line)
        if m:
            self.edges.append(EdgeInfo(
                src=m.group(1), dst=m.group(3),
                bond_type='===', stereo=m.group(2).upper(),
            ))
            return

        # Regular bond: atom1 bond_type atom2
        m = re.search(r'([\w_]+)\s*(---|\===|-\.-|-->)\s*([\w_]+)', line)
        if m:
            self.edges.append(EdgeInfo(
                src=m.group(1), dst=m.group(3),
                bond_type=m.group(2),
            ))
            return

        # Abbreviation definition (curly braces) -- must be matched before square brackets because {R[5]} contains []
        m = re.search(r'([\w_]+?)(?:_(R|S))?\{([^}]+)\}', line)
        if m:
            base_id = m.group(1)
            chirality = m.group(2) or ""
            label = m.group(3)
            node_id = f"{base_id}_{chirality}" if chirality else base_id
            if node_id not in self.nodes:
                self.nodes[node_id] = NodeInfo(
                    node_id=node_id, label=label,
                    label_type="abbrev", chirality=chirality,
                )
            return

        # Atom definition (square brackets): AtomID[Label] or AtomID_R[Label]
        m = re.search(r'([\w_]+?)(?:_(R|S))?\[([^\]]+)\]', line)
        if m:
            base_id = m.group(1)
            chirality = m.group(2) or ""
            label = m.group(3)
            node_id = f"{base_id}_{chirality}" if chirality else base_id
            if node_id not in self.nodes:
                self.nodes[node_id] = NodeInfo(
                    node_id=node_id, label=label,
                    label_type="atom", chirality=chirality,
                )
            return
            return

    def to_networkx(self) -> nx.Graph:
        """Convert to a networkx graph."""
        G = nx.Graph()
        for nid, info in self.nodes.items():
            G.add_node(nid, label=info.label, label_type=info.label_type,
                       chirality=info.chirality)
        for edge in self.edges:
            if edge.src in self.nodes and edge.dst in self.nodes:
                G.add_edge(edge.src, edge.dst,
                           bond_type=edge.bond_type, stereo=edge.stereo)
        return G

    @property
    def num_nodes(self):
        return len(self.nodes)

    @property
    def num_edges(self):
        return len(self.edges)

    def get_abbrev_labels(self) -> List[str]:
        """Return a list of labels for all abbreviation nodes."""
        return [n.label for n in self.nodes.values() if n.label_type == "abbrev"]


def normalize_abbrev_name(name: str) -> str:
    """Normalize abbreviation names: R1 = R[1] = R1 (superscript), CH[2]?n = (CH2)n, etc."""
    name = name.strip()
    # R^1 -> R1, R^a -> Ra, R^f -> Rf (caret removal)
    name = name.replace('^', '')
    # R[1] -> R1, B[5] -> B5, etc.
    name = re.sub(r'\[(\d+)\]', r'\1', name)
    # R¹ R² etc -> R1 R2 (unicode superscript)
    superscripts = str.maketrans('⁰¹²³⁴⁵⁶⁷⁸⁹', '0123456789')
    name = name.translate(superscripts)
    # R₁ R₂ etc -> R1 R2 (unicode subscript)
    subscripts = str.maketrans('₀₁₂₃₄₅₆₇₈₉', '0123456789')
    name = name.translate(subscripts)
    # CH[2]?n -> (CH2)n, CH[2]n -> (CH2)n, CH2?n -> (CH2)n, CH2n -> (CH2)n
    name = re.sub(r'CH\[2\]\??n', '(CH2)n', name)
    name = re.sub(r'^CH2\??n$', '(CH2)n', name)
    # CH2?x -> (CH2)x
    name = re.sub(r'^CH2\??x$', '(CH2)x', name)
    # R1?n -> (R1)n, R[5]?t -> (R5)t
    # (these are already handled by bracket removal above)
    # COOCH[3] -> COOCH3, etc.
    name = re.sub(r'\[(\d+)\]', r'\1', name)
    # R[2]' -> R2', X[1Z] -> X1Z
    name = re.sub(r'\[([^\]]+)\]', r'\1', name)
    # Normalize case for common abbreviations
    lower = name.lower()
    # boc -> Boc, etc.
    case_map = {
        'boc': 'Boc', 'cbz': 'Cbz', 'fmoc': 'Fmoc', 'tbdms': 'TBDMS',
        'tms': 'TMS', 'ac': 'Ac', 'ts': 'Ts', 'tf': 'Tf',
        'me': 'Me', 'et': 'Et', 'ph': 'Ph', 'bn': 'Bn',
        'ome': 'OMe', 'oet': 'OEt', 'oac': 'OAc', 'obn': 'OBn',
    }
    if lower in case_map:
        name = case_map[lower]

    # Synonym normalization: chemically equivalent abbreviation names
    synonym_map = {
        'CO2R': 'COOR', 'COOMe': 'CO2Me', 'COOEt': 'CO2Et',
        'COOCH3': 'CO2Me', 'CO2H': 'COOH',
        'MeO': 'OMe', 'OCH3': 'OMe', 'EtO': 'OEt',
        'AcHN': 'NHAc', 'MeO2C': 'CO2Me',
        'O2N': 'NO2', 'Tos': 'Ts',
        'BOC': 'Boc', 'boc': 'Boc',
    }
    if name in synonym_map:
        name = synonym_map[name]

    # Trailing comma/punctuation cleanup: B5, -> B5
    name = name.rstrip(',. ')

    return name


def _detect_aromatic_ring_edges(G: nx.Graph) -> set:
    """Detect edges on aromatic rings.

    Criterion: all edges around the ring strictly alternate === and ---,
    and the ring size is 5, 6, or 7.
    Returns the set of frozenset(src, dst) pairs marked as aromatic ring edges.
    """
    aromatic_edges = set()

    try:
        cycles = nx.cycle_basis(G)
    except Exception:
        return aromatic_edges

    for cycle in cycles:
        if len(cycle) not in (5, 6, 7):
            continue

        # Get bond types around the cycle
        edge_types = []
        for i in range(len(cycle)):
            u = cycle[i]
            v = cycle[(i + 1) % len(cycle)]
            edata = G.edges.get((u, v))
            if edata is None:
                edata = G.edges.get((v, u))
            if edata is None:
                break
            edge_types.append(edata.get('bond_type', ''))
        else:
            # Check alternating single/double pattern
            if all(t in ('===', '---') for t in edge_types):
                is_alternating = all(
                    edge_types[i] != edge_types[(i + 1) % len(edge_types)]
                    for i in range(len(edge_types))
                )
                if is_alternating:
                    for i in range(len(cycle)):
                        u = cycle[i]
                        v = cycle[(i + 1) % len(cycle)]
                        aromatic_edges.add(frozenset((u, v)))

    return aromatic_edges


def _mark_aromatic_edges(G: nx.Graph) -> nx.Graph:
    """Mark aromatic ring edges on the graph by setting the is_aromatic=True attribute."""
    G = G.copy()
    aromatic = _detect_aromatic_ring_edges(G)
    for u, v, data in G.edges(data=True):
        if frozenset((u, v)) in aromatic:
            data['is_aromatic'] = True
        else:
            data['is_aromatic'] = False
    return G


def molecode_isomorphic(g1: MoleCodeGraph, g2: MoleCodeGraph,
                   ignore_stereo: bool = True,
                   normalize_abbrevs: bool = True,
                   abbrev_expand_map: Optional[dict] = None) -> Tuple[bool, dict]:
    """MoleCode-level graph isomorphism comparison.

    Args:
        g1, g2: Two MoleCode graphs
        ignore_stereo: Whether to ignore stereochemistry (E/Z, R/S)
        normalize_abbrevs: Whether to normalize abbreviation names
        abbrev_expand_map: Abbreviation expansion map (for matching expanded vs. unexpanded cases)

    Returns:
        (is_isomorphic, details_dict)
        details_dict contains:
            - "match": bool
            - "reason": str (reason for match/mismatch)
            - "unmatched_abbrevs": list (abbreviations that could not be matched)
    """
    details = {"match": False, "reason": "", "unmatched_abbrevs": []}

    G1 = _mark_aromatic_edges(g1.to_networkx())
    G2 = _mark_aromatic_edges(g2.to_networkx())

    # Quick check
    if G1.number_of_nodes() == 0 or G2.number_of_nodes() == 0:
        details["reason"] = "empty graph"
        return False, details

    def _normalize_label(label, label_type):
        if normalize_abbrevs and label_type == "abbrev":
            return normalize_abbrev_name(label)
        return label

    def node_match(n1_attrs, n2_attrs):
        l1 = _normalize_label(n1_attrs['label'], n1_attrs['label_type'])
        l2 = _normalize_label(n2_attrs['label'], n2_attrs['label_type'])
        t1 = n1_attrs['label_type']
        t2 = n2_attrs['label_type']

        # Same type, same label
        if t1 == t2 and l1 == l2:
            if not ignore_stereo:
                return n1_attrs.get('chirality', '') == n2_attrs.get('chirality', '')
            return True

        # Same type, case-insensitive match for abbreviations
        if t1 == t2 == "abbrev" and l1.lower() == l2.lower():
            return True

        # atom [OH] vs abbrev {OH} — same label, different type notation
        if l1 == l2 and {t1, t2} == {"atom", "abbrev"}:
            return True

        # One is abbrev, other is atom — check if abbrev can expand to match
        if abbrev_expand_map and t1 != t2:
            abbrev_label = l1 if t1 == "abbrev" else l2
            atom_label = l1 if t1 == "atom" else l2
            expanded = abbrev_expand_map.get(normalize_abbrev_name(abbrev_label))
            if expanded and expanded.get("single_atom_label") == atom_label:
                return True

        return False

    def edge_match(e1_attrs, e2_attrs):
        bt1 = e1_attrs['bond_type']
        bt2 = e2_attrs['bond_type']
        ar1 = e1_attrs.get('is_aromatic', False)
        ar2 = e2_attrs.get('is_aromatic', False)

        # Both aromatic ring edges: === and --- are equivalent (Kekulé ambiguity)
        if ar1 and ar2:
            if bt1 in ('===', '---') and bt2 in ('===', '---'):
                return True

        if bt1 != bt2:
            return False
        if not ignore_stereo:
            return e1_attrs.get('stereo', '') == e2_attrs.get('stereo', '')
        return True

    # Direct isomorphism check
    is_iso = nx.is_isomorphic(G1, G2, node_match=node_match, edge_match=edge_match)

    if is_iso:
        details["match"] = True
        details["reason"] = "isomorphic"
        return True, details

    # If not matching, try with abbreviation expansion on BOTH sides
    # This handles: model kept {Boc} but GT expanded it, or vice versa
    if abbrev_expand_map:
        G1_expanded = _expand_graph(G1, abbrev_expand_map)
        G2_expanded = _expand_graph(G2, abbrev_expand_map)

        # Mark aromatic edges on expanded graphs
        G1_expanded = _mark_aromatic_edges(G1_expanded)
        G2_expanded = _mark_aromatic_edges(G2_expanded)

        # Use the full node_match (with abbrev handling) on expanded graphs
        is_iso_expanded = nx.is_isomorphic(
            G1_expanded, G2_expanded,
            node_match=node_match,
            edge_match=edge_match,
        )
        if is_iso_expanded:
            details["match"] = True
            details["reason"] = "isomorphic after expansion"
            return True, details

        # Also try: expand only one side (asymmetric case)
        for Ga, Gb, label in [(G1_expanded, G2, "expand pred only"),
                               (G1, G2_expanded, "expand gt only")]:
            Gb_marked = _mark_aromatic_edges(Gb) if Gb is not G2 else _mark_aromatic_edges(Gb)
            Ga_marked = _mark_aromatic_edges(Ga) if Ga is not G1 else Ga
            if Ga_marked.number_of_nodes() == Gb_marked.number_of_nodes():
                is_iso_asym = nx.is_isomorphic(
                    Ga_marked, Gb_marked,
                    node_match=node_match,
                    edge_match=edge_match,
                )
                if is_iso_asym:
                    details["match"] = True
                    details["reason"] = f"isomorphic after {label}"
                    return True, details

    # Collect unmatched abbreviations
    abbrevs1 = {_normalize_label(n.label, n.label_type)
                for n in g1.nodes.values() if n.label_type == "abbrev"}
    abbrevs2 = {_normalize_label(n.label, n.label_type)
                for n in g2.nodes.values() if n.label_type == "abbrev"}
    all_abbrevs = abbrevs1 | abbrevs2
    if abbrev_expand_map:
        unmatched = [a for a in all_abbrevs if normalize_abbrev_name(a) not in abbrev_expand_map]
    else:
        unmatched = list(all_abbrevs)
    details["unmatched_abbrevs"] = unmatched

    # Diagnosis
    if G1.number_of_nodes() != G2.number_of_nodes():
        details["reason"] = f"node count mismatch ({G1.number_of_nodes()} vs {G2.number_of_nodes()})"
    elif G1.number_of_edges() != G2.number_of_edges():
        details["reason"] = f"edge count mismatch ({G1.number_of_edges()} vs {G2.number_of_edges()})"
    else:
        details["reason"] = "topology mismatch"

    return False, details


def _expand_graph(G: nx.Graph, abbrev_map: dict) -> nx.Graph:
    """Expand expandable abbreviation nodes in the graph."""
    G_new = G.copy()
    nodes_to_remove = []
    node_counter = [0]

    for node, attrs in list(G_new.nodes(data=True)):
        if attrs.get('label_type') != 'abbrev':
            continue
        label = normalize_abbrev_name(attrs['label'])
        expansion = abbrev_map.get(label)
        if not expansion:
            continue

        neighbors = list(G_new.neighbors(node))

        if expansion.get("single_atom_label"):
            # Simple 1:1 replacement (e.g., Me -> CH3)
            G_new.nodes[node]['label'] = expansion["single_atom_label"]
            G_new.nodes[node]['label_type'] = "atom"
        elif expansion.get("subgraph"):
            # Multi-atom expansion
            sub = expansion["subgraph"]
            prefix = f"_exp{node_counter[0]}_"
            node_counter[0] += 1

            # Add subgraph nodes
            attach_node = None
            for sub_id, sub_label in sub["atoms"]:
                new_id = f"{prefix}{sub_id}"
                G_new.add_node(new_id, label=sub_label, label_type="atom", chirality="")
                if sub_id == sub.get("attach"):
                    attach_node = new_id

            # Add subgraph edges
            for s, d, bt in sub["bonds"]:
                G_new.add_edge(f"{prefix}{s}", f"{prefix}{d}", bond_type=bt, stereo="")

            # Connect attachment point to original neighbors
            if attach_node and neighbors:
                for nb in neighbors:
                    edge_data = G_new.edges[node, nb]
                    G_new.add_edge(attach_node, nb, **edge_data)

            nodes_to_remove.append(node)

    for node in nodes_to_remove:
        G_new.remove_node(node)

    return G_new


if __name__ == "__main__":
    # Quick test
    mc1 = """
    graph TB
        subgraph Mol["test"]
            Mol_C_1[C]
            Mol_C_2[CH]
            Mol_C_3[CH]
            Mol_C_1 === Mol_C_2
            Mol_C_2 --- Mol_C_3
            Mol_C_3 === Mol_C_1
        end
    """
    mc2 = """
    graph TB
        subgraph M["test"]
            M_C_1[CH]
            M_C_2[C]
            M_C_3[CH]
            M_C_1 --- M_C_2
            M_C_2 === M_C_3
            M_C_3 === M_C_1
        end
    """
    g1 = MoleCodeGraph.from_text(mc1)
    g2 = MoleCodeGraph.from_text(mc2)
    print(f"g1: {g1.num_nodes} nodes, {g1.num_edges} edges")
    print(f"g2: {g2.num_nodes} nodes, {g2.num_edges} edges")
    is_iso, details = molecode_isomorphic(g1, g2)
    print(f"Isomorphic: {is_iso}, reason: {details['reason']}")

    # Test with abbreviation
    mc3 = """
    graph TB
        subgraph Mol["test"]
            Mol_C_1[C]
            Mol_X_1{Me}
            Mol_C_1 --- Mol_X_1
        end
    """
    mc4 = """
    graph TB
        subgraph Mol["test"]
            Mol_C_1[C]
            Mol_X_1{CH3}
            Mol_C_1 --- Mol_X_1
        end
    """
    g3 = MoleCodeGraph.from_text(mc3)
    g4 = MoleCodeGraph.from_text(mc4)
    abbrev_map = {"Me": {"single_atom_label": "CH3"}, "CH3": {"single_atom_label": "CH3"}}
    is_iso, details = molecode_isomorphic(g3, g4, abbrev_expand_map=abbrev_map)
    print(f"Me vs CH3: {is_iso}, reason: {details['reason']}")
