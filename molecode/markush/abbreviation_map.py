"""
Abbreviation group -> MoleCode subgraph expansion map.

Two expansion types:
1. single_atom_label: the abbreviation is equivalent to a single atom label (e.g. Me -> [CH3])
2. subgraph: the abbreviation expands to a multi-atom + bond subgraph (e.g. Boc -> OC(=O)C(C)(C)C)

R groups (R, R1, R2...) and generic placeholders (X, Y, Z, Ar) are not expanded;
they are kept as abbreviation labels because they represent variable groups with no unique expansion.
"""

# ============================================================
# Single-atom expansion (abbreviation <-> single atom label)
# ============================================================

SINGLE_ATOM_MAP = {
    # Alkyl
    "Me": "CH3",
    "CH3": "CH3",
    "Et": "CH2",       # Et as a substituent: the connecting end is CH2
    "CH2CH3": "CH2",

    # Halogens / small groups (already atom-like)
    "F": "F",
    "Cl": "Cl",
    "Br": "Br",
    "I": "I",

    # Functional groups that map to single atom with properties
    "CN": "C",          # -C=N: the connecting carbon (subgraph expansion is more accurate, but single-atom is a useful approximation)
    "NC": "N",
    "N3": "N",
    "NO": "N",
    "CHO": "CH",        # Aldehyde carbon

    # Charged atoms
    "NH3+": "NH3(+)",
    "NH2+": "NH2(+)",
    "O-": "O(-)",
    "N+": "N(+)",
    "SO3-": "S",

    # Additional from Markush data
    "H": "H",
    "N": "N",
    "C": "C",
    "O": "O",
    "B": "B",
    "OH": "OH",
}

# ============================================================
# Subgraph expansion (abbreviation -> multi-atom subgraph)
# ============================================================
# Format: {"atoms": [(id, label), ...], "bonds": [(id1, id2, type), ...], "attach": id}
# attach is the attachment point (the atom that connects to the parent graph)

SUBGRAPH_MAP = {
    "Boc": {
        "atoms": [("C1", "C"), ("O1", "O"), ("O2", "O"), ("C2", "C"), ("C3", "CH3"), ("C4", "CH3"), ("C5", "CH3")],
        "bonds": [("C1", "O1", "==="), ("C1", "O2", "---"), ("O2", "C2", "---"),
                  ("C2", "C3", "---"), ("C2", "C4", "---"), ("C2", "C5", "---")],
        "attach": "C1",
    },
    "BOC": {  # alias
        "atoms": [("C1", "C"), ("O1", "O"), ("O2", "O"), ("C2", "C"), ("C3", "CH3"), ("C4", "CH3"), ("C5", "CH3")],
        "bonds": [("C1", "O1", "==="), ("C1", "O2", "---"), ("O2", "C2", "---"),
                  ("C2", "C3", "---"), ("C2", "C4", "---"), ("C2", "C5", "---")],
        "attach": "C1",
    },
    "boc": {  # alias
        "atoms": [("C1", "C"), ("O1", "O"), ("O2", "O"), ("C2", "C"), ("C3", "CH3"), ("C4", "CH3"), ("C5", "CH3")],
        "bonds": [("C1", "O1", "==="), ("C1", "O2", "---"), ("O2", "C2", "---"),
                  ("C2", "C3", "---"), ("C2", "C4", "---"), ("C2", "C5", "---")],
        "attach": "C1",
    },
    "Cbz": {
        "atoms": [("C1", "C"), ("O1", "O"), ("O2", "O"), ("C2", "CH2"),
                  ("C3", "C"), ("C4", "CH"), ("C5", "CH"), ("C6", "CH"), ("C7", "CH"), ("C8", "CH")],
        "bonds": [("C1", "O1", "==="), ("C1", "O2", "---"), ("O2", "C2", "---"), ("C2", "C3", "---"),
                  ("C3", "C4", "==="), ("C4", "C5", "---"), ("C5", "C6", "==="),
                  ("C6", "C7", "---"), ("C7", "C8", "==="), ("C8", "C3", "---")],
        "attach": "C1",
    },
    "Ac": {
        "atoms": [("C1", "C"), ("O1", "O"), ("C2", "CH3")],
        "bonds": [("C1", "O1", "==="), ("C1", "C2", "---")],
        "attach": "C1",
    },
    "Ts": {
        "atoms": [("S1", "S"), ("O1", "O"), ("O2", "O"), ("C1", "C"), ("C2", "CH"), ("C3", "CH"),
                  ("C4", "C"), ("C5", "CH"), ("C6", "CH"), ("C7", "CH3")],
        "bonds": [("S1", "O1", "==="), ("S1", "O2", "==="), ("S1", "C1", "---"),
                  ("C1", "C2", "==="), ("C2", "C3", "---"), ("C3", "C4", "==="),
                  ("C4", "C5", "---"), ("C5", "C6", "==="), ("C6", "C1", "---"), ("C4", "C7", "---")],
        "attach": "S1",
    },
    "Tos": {  # alias for Ts
        "atoms": [("S1", "S"), ("O1", "O"), ("O2", "O"), ("C1", "C"), ("C2", "CH"), ("C3", "CH"),
                  ("C4", "C"), ("C5", "CH"), ("C6", "CH"), ("C7", "CH3")],
        "bonds": [("S1", "O1", "==="), ("S1", "O2", "==="), ("S1", "C1", "---"),
                  ("C1", "C2", "==="), ("C2", "C3", "---"), ("C3", "C4", "==="),
                  ("C4", "C5", "---"), ("C5", "C6", "==="), ("C6", "C1", "---"), ("C4", "C7", "---")],
        "attach": "S1",
    },
    "Ph": {
        "atoms": [("C1", "C"), ("C2", "CH"), ("C3", "CH"), ("C4", "CH"), ("C5", "CH"), ("C6", "CH")],
        "bonds": [("C1", "C2", "==="), ("C2", "C3", "---"), ("C3", "C4", "==="),
                  ("C4", "C5", "---"), ("C5", "C6", "==="), ("C6", "C1", "---")],
        "attach": "C1",
    },
    "Bn": {
        "atoms": [("C0", "CH2"), ("C1", "C"), ("C2", "CH"), ("C3", "CH"),
                  ("C4", "CH"), ("C5", "CH"), ("C6", "CH")],
        "bonds": [("C0", "C1", "---"), ("C1", "C2", "==="), ("C2", "C3", "---"),
                  ("C3", "C4", "==="), ("C4", "C5", "---"), ("C5", "C6", "==="), ("C6", "C1", "---")],
        "attach": "C0",
    },
    "CF3": {
        "atoms": [("C1", "C"), ("F1", "F"), ("F2", "F"), ("F3", "F")],
        "bonds": [("C1", "F1", "---"), ("C1", "F2", "---"), ("C1", "F3", "---")],
        "attach": "C1",
    },
    "CHF2": {
        "atoms": [("C1", "CH"), ("F1", "F"), ("F2", "F")],
        "bonds": [("C1", "F1", "---"), ("C1", "F2", "---")],
        "attach": "C1",
    },
    "OCF3": {
        "atoms": [("O1", "O"), ("C1", "C"), ("F1", "F"), ("F2", "F"), ("F3", "F")],
        "bonds": [("O1", "C1", "---"), ("C1", "F1", "---"), ("C1", "F2", "---"), ("C1", "F3", "---")],
        "attach": "O1",
    },
    "OMe": {
        "atoms": [("O1", "O"), ("C1", "CH3")],
        "bonds": [("O1", "C1", "---")],
        "attach": "O1",
    },
    "OCH3": {  # alias for OMe
        "atoms": [("O1", "O"), ("C1", "CH3")],
        "bonds": [("O1", "C1", "---")],
        "attach": "O1",
    },
    "OEt": {
        "atoms": [("O1", "O"), ("C1", "CH2"), ("C2", "CH3")],
        "bonds": [("O1", "C1", "---"), ("C1", "C2", "---")],
        "attach": "O1",
    },
    "OAc": {
        "atoms": [("O1", "O"), ("C1", "C"), ("O2", "O"), ("C2", "CH3")],
        "bonds": [("O1", "C1", "---"), ("C1", "O2", "==="), ("C1", "C2", "---")],
        "attach": "O1",
    },
    "OBn": {
        "atoms": [("O1", "O"), ("C0", "CH2"), ("C1", "C"), ("C2", "CH"), ("C3", "CH"),
                  ("C4", "CH"), ("C5", "CH"), ("C6", "CH")],
        "bonds": [("O1", "C0", "---"), ("C0", "C1", "---"),
                  ("C1", "C2", "==="), ("C2", "C3", "---"), ("C3", "C4", "==="),
                  ("C4", "C5", "---"), ("C5", "C6", "==="), ("C6", "C1", "---")],
        "attach": "O1",
    },
    "OTs": {
        "atoms": [("O1", "O"), ("S1", "S"), ("O2", "O"), ("O3", "O"),
                  ("C1", "C"), ("C2", "CH"), ("C3", "CH"), ("C4", "C"), ("C5", "CH"), ("C6", "CH"), ("C7", "CH3")],
        "bonds": [("O1", "S1", "---"), ("S1", "O2", "==="), ("S1", "O3", "==="), ("S1", "C1", "---"),
                  ("C1", "C2", "==="), ("C2", "C3", "---"), ("C3", "C4", "==="),
                  ("C4", "C5", "---"), ("C5", "C6", "==="), ("C6", "C1", "---"), ("C4", "C7", "---")],
        "attach": "O1",
    },
    "OTf": {
        "atoms": [("O1", "O"), ("S1", "S"), ("O2", "O"), ("O3", "O"),
                  ("C1", "C"), ("F1", "F"), ("F2", "F"), ("F3", "F")],
        "bonds": [("O1", "S1", "---"), ("S1", "O2", "==="), ("S1", "O3", "==="), ("S1", "C1", "---"),
                  ("C1", "F1", "---"), ("C1", "F2", "---"), ("C1", "F3", "---")],
        "attach": "O1",
    },
    "Tf": {
        "atoms": [("S1", "S"), ("O1", "O"), ("O2", "O"), ("C1", "C"), ("F1", "F"), ("F2", "F"), ("F3", "F")],
        "bonds": [("S1", "O1", "==="), ("S1", "O2", "==="), ("S1", "C1", "---"),
                  ("C1", "F1", "---"), ("C1", "F2", "---"), ("C1", "F3", "---")],
        "attach": "S1",
    },
    "NO2": {
        "atoms": [("N1", "N(+)"), ("O1", "O"), ("O2", "O(-)")],
        "bonds": [("N1", "O1", "==="), ("N1", "O2", "---")],
        "attach": "N1",
    },
    "COOH": {
        "atoms": [("C1", "C"), ("O1", "O"), ("O2", "OH")],
        "bonds": [("C1", "O1", "==="), ("C1", "O2", "---")],
        "attach": "C1",
    },
    "CO2H": {  # alias
        "atoms": [("C1", "C"), ("O1", "O"), ("O2", "OH")],
        "bonds": [("C1", "O1", "==="), ("C1", "O2", "---")],
        "attach": "C1",
    },
    "COOMe": {
        "atoms": [("C1", "C"), ("O1", "O"), ("O2", "O"), ("C2", "CH3")],
        "bonds": [("C1", "O1", "==="), ("C1", "O2", "---"), ("O2", "C2", "---")],
        "attach": "C1",
    },
    "CO2Me": {  # alias
        "atoms": [("C1", "C"), ("O1", "O"), ("O2", "O"), ("C2", "CH3")],
        "bonds": [("C1", "O1", "==="), ("C1", "O2", "---"), ("O2", "C2", "---")],
        "attach": "C1",
    },
    "COOCH3": {  # alias
        "atoms": [("C1", "C"), ("O1", "O"), ("O2", "O"), ("C2", "CH3")],
        "bonds": [("C1", "O1", "==="), ("C1", "O2", "---"), ("O2", "C2", "---")],
        "attach": "C1",
    },
    "COOEt": {
        "atoms": [("C1", "C"), ("O1", "O"), ("O2", "O"), ("C2", "CH2"), ("C3", "CH3")],
        "bonds": [("C1", "O1", "==="), ("C1", "O2", "---"), ("O2", "C2", "---"), ("C2", "C3", "---")],
        "attach": "C1",
    },
    "CO2Et": {  # alias
        "atoms": [("C1", "C"), ("O1", "O"), ("O2", "O"), ("C2", "CH2"), ("C3", "CH3")],
        "bonds": [("C1", "O1", "==="), ("C1", "O2", "---"), ("O2", "C2", "---"), ("C2", "C3", "---")],
        "attach": "C1",
    },
    "NHBoc": {
        "atoms": [("N1", "NH"), ("C1", "C"), ("O1", "O"), ("O2", "O"),
                  ("C2", "C"), ("C3", "CH3"), ("C4", "CH3"), ("C5", "CH3")],
        "bonds": [("N1", "C1", "---"), ("C1", "O1", "==="), ("C1", "O2", "---"), ("O2", "C2", "---"),
                  ("C2", "C3", "---"), ("C2", "C4", "---"), ("C2", "C5", "---")],
        "attach": "N1",
    },
    "NHCbz": {
        "atoms": [("N1", "NH"), ("C1", "C"), ("O1", "O"), ("O2", "O"), ("C2", "CH2"),
                  ("C3", "C"), ("C4", "CH"), ("C5", "CH"), ("C6", "CH"), ("C7", "CH"), ("C8", "CH")],
        "bonds": [("N1", "C1", "---"), ("C1", "O1", "==="), ("C1", "O2", "---"), ("O2", "C2", "---"),
                  ("C2", "C3", "---"), ("C3", "C4", "==="), ("C4", "C5", "---"), ("C5", "C6", "==="),
                  ("C6", "C7", "---"), ("C7", "C8", "==="), ("C8", "C3", "---")],
        "attach": "N1",
    },
    "NHAc": {
        "atoms": [("N1", "NH"), ("C1", "C"), ("O1", "O"), ("C2", "CH3")],
        "bonds": [("N1", "C1", "---"), ("C1", "O1", "==="), ("C1", "C2", "---")],
        "attach": "N1",
    },
    "NHMe": {
        "atoms": [("N1", "NH"), ("C1", "CH3")],
        "bonds": [("N1", "C1", "---")],
        "attach": "N1",
    },
    "NMe2": {
        "atoms": [("N1", "N"), ("C1", "CH3"), ("C2", "CH3")],
        "bonds": [("N1", "C1", "---"), ("N1", "C2", "---")],
        "attach": "N1",
    },
    "NHOH": {
        "atoms": [("N1", "NH"), ("O1", "OH")],
        "bonds": [("N1", "O1", "---")],
        "attach": "N1",
    },
    "SMe": {
        "atoms": [("S1", "S"), ("C1", "CH3")],
        "bonds": [("S1", "C1", "---")],
        "attach": "S1",
    },
    "SO2Me": {
        "atoms": [("S1", "S"), ("O1", "O"), ("O2", "O"), ("C1", "CH3")],
        "bonds": [("S1", "O1", "==="), ("S1", "O2", "==="), ("S1", "C1", "---")],
        "attach": "S1",
    },
    "SO3H": {
        "atoms": [("S1", "S"), ("O1", "O"), ("O2", "O"), ("O3", "OH")],
        "bonds": [("S1", "O1", "==="), ("S1", "O2", "==="), ("S1", "O3", "---")],
        "attach": "S1",
    },
    "ONO2": {
        "atoms": [("O1", "O"), ("N1", "N(+)"), ("O2", "O"), ("O3", "O(-)")],
        "bonds": [("O1", "N1", "---"), ("N1", "O2", "==="), ("N1", "O3", "---")],
        "attach": "O1",
    },
    "TMS": {
        "atoms": [("Si1", "Si"), ("C1", "CH3"), ("C2", "CH3"), ("C3", "CH3")],
        "bonds": [("Si1", "C1", "---"), ("Si1", "C2", "---"), ("Si1", "C3", "---")],
        "attach": "Si1",
    },
    "tBu": {
        "atoms": [("C1", "C"), ("C2", "CH3"), ("C3", "CH3"), ("C4", "CH3")],
        "bonds": [("C1", "C2", "---"), ("C1", "C3", "---"), ("C1", "C4", "---")],
        "attach": "C1",
    },
    "iPr": {
        "atoms": [("C1", "CH"), ("C2", "CH3"), ("C3", "CH3")],
        "bonds": [("C1", "C2", "---"), ("C1", "C3", "---")],
        "attach": "C1",
    },
    "Bz": {
        "atoms": [("C1", "C"), ("O1", "O"),
                  ("C2", "C"), ("C3", "CH"), ("C4", "CH"), ("C5", "CH"), ("C6", "CH"), ("C7", "CH")],
        "bonds": [("C1", "O1", "==="), ("C1", "C2", "---"),
                  ("C2", "C3", "==="), ("C3", "C4", "---"), ("C4", "C5", "==="),
                  ("C5", "C6", "---"), ("C6", "C7", "==="), ("C7", "C2", "---")],
        "attach": "C1",
    },
    "NCO": {
        "atoms": [("N1", "N"), ("C1", "C"), ("O1", "O")],
        "bonds": [("N1", "C1", "==="), ("C1", "O1", "===")],
        "attach": "N1",
    },
    "Fmoc": {
        # Simplified: treat as single abbreviation node (too complex to expand)
        "atoms": [("X1", "Fmoc")],
        "bonds": [],
        "attach": "X1",
    },
    "SEM": {
        # 2-(Trimethylsilyl)ethoxymethyl - complex, keep as abbreviation
        "atoms": [("X1", "SEM")],
        "bonds": [],
        "attach": "X1",
    },
    "Bpin": {
        # Pinacol boronate - complex, keep as abbreviation
        "atoms": [("X1", "Bpin")],
        "bonds": [],
        "attach": "X1",
    },
    # Additional from Markush data
    "MeO": {  # alias for OMe (different order)
        "atoms": [("O1", "O"), ("C1", "CH3")],
        "bonds": [("O1", "C1", "---")],
        "attach": "O1",
    },
    "CONH2": {
        "atoms": [("C1", "C"), ("O1", "O"), ("N1", "NH2")],
        "bonds": [("C1", "O1", "==="), ("C1", "N1", "---")],
        "attach": "C1",
    },
    "CONHCH3": {
        "atoms": [("C1", "C"), ("O1", "O"), ("N1", "NH"), ("C2", "CH3")],
        "bonds": [("C1", "O1", "==="), ("C1", "N1", "---"), ("N1", "C2", "---")],
        "attach": "C1",
    },
    "CONR1R2": {
        "atoms": [("C1", "C"), ("O1", "O"), ("N1", "N"), ("X1", "R1"), ("X2", "R2")],
        "bonds": [("C1", "O1", "==="), ("C1", "N1", "---"), ("N1", "X1", "---"), ("N1", "X2", "---")],
        "attach": "C1",
    },
    "AcHN": {  # = NHAc reversed
        "atoms": [("N1", "NH"), ("C1", "C"), ("O1", "O"), ("C2", "CH3")],
        "bonds": [("N1", "C1", "---"), ("C1", "O1", "==="), ("C1", "C2", "---")],
        "attach": "N1",
    },
    "MeO2C": {  # = CO2Me reversed attachment
        "atoms": [("C1", "C"), ("O1", "O"), ("O2", "O"), ("C2", "CH3")],
        "bonds": [("C1", "O1", "==="), ("C1", "O2", "---"), ("O2", "C2", "---")],
        "attach": "C1",
    },
    "O2N": {  # = NO2 reversed attachment
        "atoms": [("N1", "N(+)"), ("O1", "O"), ("O2", "O(-)")],
        "bonds": [("N1", "O1", "==="), ("N1", "O2", "---")],
        "attach": "N1",
    },
    "CF2": {
        "atoms": [("C1", "C"), ("F1", "F"), ("F2", "F")],
        "bonds": [("C1", "F1", "---"), ("C1", "F2", "---")],
        "attach": "C1",
    },
    "CCCH3": {  # propyl-like
        "atoms": [("C1", "CH2"), ("C2", "CH2"), ("C3", "CH3")],
        "bonds": [("C1", "C2", "---"), ("C2", "C3", "---")],
        "attach": "C1",
    },
    "CH(C6H5)2": {
        # Diphenylmethyl - complex, keep as abbreviation
        "atoms": [("X1", "CH(C6H5)2")],
        "bonds": [],
        "attach": "X1",
    },
    "DMT": {
        # Dimethoxytrityl - complex, keep as abbreviation
        "atoms": [("X1", "DMT")],
        "bonds": [],
        "attach": "X1",
    },
    "DMTO": {
        "atoms": [("X1", "DMT"), ("O1", "O")],
        "bonds": [("X1", "O1", "---")],
        "attach": "O1",
    },
    "TBS": {  # = TBDMS alias
        "atoms": [("Si1", "Si"), ("C1", "C"), ("C2", "CH3"), ("C3", "CH3"), ("C4", "CH3"), ("C5", "CH3")],
        "bonds": [("Si1", "C1", "---"), ("Si1", "C2", "---"), ("Si1", "C3", "---"),
                  ("C1", "C4", "---"), ("C1", "C5", "---")],
        "attach": "Si1",
    },
    "TBSO": {
        "atoms": [("O1", "O"), ("Si1", "Si"), ("C1", "C"), ("C2", "CH3"), ("C3", "CH3"), ("C4", "CH3"), ("C5", "CH3")],
        "bonds": [("O1", "Si1", "---"), ("Si1", "C1", "---"), ("Si1", "C2", "---"), ("Si1", "C3", "---"),
                  ("C1", "C4", "---"), ("C1", "C5", "---")],
        "attach": "O1",
    },
    "TIPS": {
        "atoms": [("Si1", "Si"), ("C1", "CH"), ("C2", "CH3"), ("C3", "CH3"),
                  ("C4", "CH"), ("C5", "CH3"), ("C6", "CH3"),
                  ("C7", "CH"), ("C8", "CH3"), ("C9", "CH3")],
        "bonds": [("Si1", "C1", "---"), ("C1", "C2", "---"), ("C1", "C3", "---"),
                  ("Si1", "C4", "---"), ("C4", "C5", "---"), ("C4", "C6", "---"),
                  ("Si1", "C7", "---"), ("C7", "C8", "---"), ("C7", "C9", "---")],
        "attach": "Si1",
    },
    "PMB": {
        # p-Methoxybenzyl
        "atoms": [("C0", "CH2"), ("C1", "C"), ("C2", "CH"), ("C3", "CH"),
                  ("C4", "C"), ("C5", "CH"), ("C6", "CH"), ("O1", "O"), ("C7", "CH3")],
        "bonds": [("C0", "C1", "---"), ("C1", "C2", "==="), ("C2", "C3", "---"),
                  ("C3", "C4", "==="), ("C4", "C5", "---"), ("C5", "C6", "==="),
                  ("C6", "C1", "---"), ("C4", "O1", "---"), ("O1", "C7", "---")],
        "attach": "C0",
    },
    "Tr": {
        # Trityl - complex, keep as abbreviation
        "atoms": [("X1", "Tr")],
        "bonds": [],
        "attach": "X1",
    },
    "EOM": {
        # Ethoxymethyl
        "atoms": [("C1", "CH2"), ("O1", "O"), ("C2", "CH2"), ("C3", "CH3")],
        "bonds": [("C1", "O1", "---"), ("O1", "C2", "---"), ("C2", "C3", "---")],
        "attach": "C1",
    },
    "OEOM": {
        "atoms": [("O1", "O"), ("C1", "CH2"), ("O2", "O"), ("C2", "CH2"), ("C3", "CH3")],
        "bonds": [("O1", "C1", "---"), ("C1", "O2", "---"), ("O2", "C2", "---"), ("C2", "C3", "---")],
        "attach": "O1",
    },
    "OBoc": {
        "atoms": [("O1", "O"), ("C1", "C"), ("O2", "O"), ("O3", "O"),
                  ("C2", "C"), ("C3", "CH3"), ("C4", "CH3"), ("C5", "CH3")],
        "bonds": [("O1", "C1", "---"), ("C1", "O2", "==="), ("C1", "O3", "---"),
                  ("O3", "C2", "---"), ("C2", "C3", "---"), ("C2", "C4", "---"), ("C2", "C5", "---")],
        "attach": "O1",
    },
    "CO2R": {
        # Ester with generic R - keep as abbreviation
        "atoms": [("X1", "CO2R")],
        "bonds": [],
        "attach": "X1",
    },
    "COOR": {
        "atoms": [("X1", "COOR")],
        "bonds": [],
        "attach": "X1",
    },
    "SO2": {
        "atoms": [("S1", "S"), ("O1", "O"), ("O2", "O")],
        "bonds": [("S1", "O1", "==="), ("S1", "O2", "===")],
        "attach": "S1",
    },
    # Fold patterns: model may merge atom+abbrev into one abbreviation
    "OR": {
        "atoms": [("O1", "O"), ("X1", "R")],
        "bonds": [("O1", "X1", "---")],
        "attach": "O1",
    },
    "RO": {
        "atoms": [("X1", "R"), ("O1", "O")],
        "bonds": [("X1", "O1", "---")],
        "attach": "O1",
    },
    "CONR1R2": {
        "atoms": [("C1", "C"), ("O1", "O"), ("N1", "N"), ("X1", "R1"), ("X2", "R2")],
        "bonds": [("C1", "O1", "==="), ("C1", "N1", "---"), ("N1", "X1", "---"), ("N1", "X2", "---")],
        "attach": "C1",
    },
    "DMTO": {
        "atoms": [("X1", "DMT"), ("O1", "O")],
        "bonds": [("X1", "O1", "---")],
        "attach": "O1",
    },
    "ZBCHN": {
        "atoms": [("X1", "ZBC"), ("N1", "NH")],
        "bonds": [("X1", "N1", "---")],
        "attach": "N1",
    },
    "NPG2": {
        "atoms": [("N1", "N"), ("X1", "PG2")],
        "bonds": [("N1", "X1", "---")],
        "attach": "N1",
    },
    "PG1N": {
        "atoms": [("X1", "PG1"), ("N1", "N")],
        "bonds": [("X1", "N1", "---")],
        "attach": "N1",
    },
    "R27NH": {
        "atoms": [("X1", "R27"), ("N1", "NH")],
        "bonds": [("X1", "N1", "---")],
        "attach": "N1",
    },
    "OCH3": {
        "atoms": [("O1", "O"), ("C1", "CH3")],
        "bonds": [("O1", "C1", "---")],
        "attach": "O1",
    },
    "OR^d": {
        "atoms": [("O1", "O"), ("X1", "Rd")],
        "bonds": [("O1", "X1", "---")],
        "attach": "O1",
    },
    "OPG1": {
        "atoms": [("O1", "O"), ("X1", "PG1")],
        "bonds": [("O1", "X1", "---")],
        "attach": "O1",
    },
    "R^a": {
        "atoms": [("X1", "Ra")],
        "bonds": [],
        "attach": "X1",
    },
}

# ============================================================
# Non-expandable labels (R groups, generic placeholders)
# These are matched directly by name during graph isomorphism comparison
# ============================================================

NON_EXPANDABLE = {
    "R", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9",
    "R10", "R11", "R12", "R13", "R14", "R15", "R16", "R17",
    "R1a", "R1b", "R2a", "R2b", "R3a", "R4a", "R5a",
    "R'", "Ra", "Rb", "Rf", "Rg",
    "R51", "R0",
    "R1-1",
    "R15'", "R16'", "R16\"", "R17'",
    "R^1", "R^2", "R^f", "R^g", "R^N1",
    "R27NH", "Y3R27",
    "X", "X1", "X1z", "Y", "Z", "Z1", "Z2", "W",
    "A", "A1", "A2", "B", "B1", "D", "E", "Q", "Q1", "Q2",
    "Ar", "Ar1", "Ar2",
    "PG", "FG",
    "Het 1", "Het 2",
    "Alkyl", "Cyc", "Chx",
    "dum",
    # Variable chains
    "(CH2)n", "(CH2)x", "(CH2)f11", "(CH2)m-1",
    "CH[2]n", "CH[2]?n", "CH2n", "CH2?n", "CH2?x", "CH2",
    "(R1)n", "(R5)t", "()n", "R1?n",
    "X(CH2)n", "[OCH2CH2]2", "OCH2CH22",
    # Generic / model-invented
    "HX", "NR", "ORd", "Rd",
    "G2", "G3", "G4",
    "x",
    "O(N)", "X/B(OH)2",
    "PG1", "PG1N", "PG2", "NPG2", "RN1", "R27", "Y8", "X1Z",
    # Round 2 additions
    "R^a", "R^0", "R2'",
    "A4", "B3", "B4", "B5", "B5,", "D1",
    "ZBC",
    "PD3",
    "CONHCH",
}


def build_expand_map():
    """Build the complete expansion map for use by graph.molecode_isomorphic."""
    expand_map = {}

    # Single atom expansions
    for name, label in SINGLE_ATOM_MAP.items():
        expand_map[name] = {"single_atom_label": label}

    # Subgraph expansions
    for name, sub in SUBGRAPH_MAP.items():
        expand_map[name] = {"subgraph": sub}

    # Non-expandable: map to themselves (no expansion needed, just name matching)
    for name in NON_EXPANDABLE:
        expand_map[name] = {"keep": True}

    return expand_map


# Singleton
EXPAND_MAP = build_expand_map()


if __name__ == "__main__":
    print(f"Single atom mappings: {len(SINGLE_ATOM_MAP)}")
    print(f"Subgraph mappings: {len(SUBGRAPH_MAP)}")
    print(f"Non-expandable: {len(NON_EXPANDABLE)}")
    print(f"Total expand map: {len(EXPAND_MAP)}")
