BASE_INSTRUCTION = '''
# Mermaid Molecular Graph Syntax Specification

> Complete syntax definition and parser requirements for Mermaid molecular graphs based on the MolToGraph project

## Syntax Overview

**Design Principles**:
- Syntax must be explicit and unambiguous for regex parsing
- Atom identifiers must be unique and traceable to chemical elements
- Bond representation must map to RDKit bond types
- Support encoding of stereochemical information

---

## Standard Example: para-Chlorophenol

**Molecular Formula**: C₆H₅ClO
**Structure**: Benzene ring + OH at C1 + Cl at C4

```mermaid
graph TB
    %% para-chlorophenol molecule
    subgraph chlorophenol["para-chlorophenol"]
        %% C1: connected to -OH, must use [C]
        chlorophenol_C_1[C]
        chlorophenol_O_1[OH]

        %% C2, C3: only ring connections, must use [CH]
        chlorophenol_C_2[CH]
        chlorophenol_C_3[CH]

        %% C4: connected to -Cl, must use [C]
        chlorophenol_C_4[C]
        chlorophenol_Cl_1[Cl]

        %% C5, C6: only ring connections, must use [CH]
        chlorophenol_C_5[CH]
        chlorophenol_C_6[CH]

        %% Kekulé form benzene ring (alternating single/double bonds)
        chlorophenol_C_1 === chlorophenol_C_2
        chlorophenol_C_2 --- chlorophenol_C_3
        chlorophenol_C_3 === chlorophenol_C_4
        chlorophenol_C_4 --- chlorophenol_C_5
        chlorophenol_C_5 === chlorophenol_C_6
        chlorophenol_C_6 --- chlorophenol_C_1

        %% Substituent single bonds
        chlorophenol_C_1 --- chlorophenol_O_1
        chlorophenol_C_4 --- chlorophenol_Cl_1
    end
```

---

## Complete Syntax Definition

### 1. Document Structure (Required)

**BNF Representation**:
```bnf
<document>       ::= <graph-header> <content>
<graph-header>   ::= "graph" <orientation>
<orientation>    ::= "TB" | "LR"
<content>        ::= (<comment> | <subgraph> | <connection>)*
```

**Constraints**:
1. First non-comment line must be `graph TB` or `graph LR`
2. Each graph must contain at least one subgraph
3. Inter-subgraph connections must be defined after all subgraph definitions

**Parser Requirements**:
- Parser skips lines starting with `graph`
- Parser skips lines starting with `%%` (comments)

---

### 2. Atom Identifier Specification (Mandatory)

**Format Definition**:
```bnf
<atom-id> ::= <subgraph-prefix> "_" <element> "_" <number> [<chirality>]
```

**Components**:

| Part | Regex | Required | Constraints |
|------|-------|----------|-------------|
| `<subgraph-prefix>` | `[A-Za-z][A-Za-z0-9]*` | Yes | Start with ASCII letter, alphanumeric only |
| `<element>` | `[A-Z][a-z]?` | Yes | Valid element symbol (C, O, N, Cl, Br...) |
| `<number>` | `[1-9][0-9]*` | Yes | Positive integer, numbered independently per element type |
| `<chirality>` | `_(R\|S)` | Optional | Chirality identifier |

**Valid Examples**:
```
chlorophenol_C_1        # Valid: standard format
ethanol_O_2             # Valid: 2nd oxygen atom
alanine_C_2_R           # Valid: R-configuration chiral carbon
E2butene_C_3            # Valid: digit can be in prefix
```

**Invalid Examples**:
```
123mol_C_1              # Invalid: prefix cannot start with digit
mol-name_C_1            # Invalid: prefix contains hyphen
mol_C_0                 # Invalid: numbering cannot be 0
mol_Xy_1                # Invalid: Xy is not a valid element symbol
mol_C_1_X               # Invalid: X is not a valid chirality identifier
```

**Parser Behavior**:
```python
# Regex for parsing atom definitions
atom_pattern = r'([\w_]+?)(?:_(R|S))?\[([^\]]+)\]'
# Capture groups:
#   group(1): base_id (prefix_element_number)
#   group(2): chirality (R or S, optional)
#   group(3): label (display label)
```

---

### 3. Atom Display Label Specification (Required)

**Format Definition**:
```bnf
<atom-node>      ::= <atom-id> "[" <display-label> "]"
<display-label>  ::= <element> [<h-count>] [<charge>]
<h-count>        ::= "H" [<digit>]
<charge>         ::= "(" [<digit>] ("+" | "-") ")"
```

**Parsing Rules** :
```python
# Regex: ^([A-Z][a-z]?)(?:H(\d*))?(?:\((\d*[+-])\))?$
# Capture groups:
#   group(1): element (element symbol)
#   group(2): h_count (hydrogen count, optional, empty string means 1)
#   group(3): charge (charge, optional, in parentheses with number+sign, e.g., "+", "2-")
# Examples:
#   "C"       → (element='C', H=0, charge=0)
#   "CH3"     → (element='C', H=3, charge=0)
#   "OH"      → (element='O', H=1, charge=0)
#   "NH2"     → (element='N', H=2, charge=0)
#   "NH2(+)"  → (element='N', H=2, charge=+1)
#   "O(-)"    → (element='O', H=0, charge=-1)
#   "N(+)"    → (element='N', H=0, charge=+1)
#   "NH3(+)"  → (element='N', H=3, charge=+1)
#   "O(2-)"   → (element='O', H=0, charge=-2)
```

**Bonding Capacity Constraints (Key Constraints)**:

Note: Apparent bonding number = Σ(bond orders) - formal charge

| Element | Typical Bonds | Label Examples | Chemical Meaning |
|---------|---------------|----------------|------------------|
| C | 4 bonds | `[CH3]`, `[CH2]`, `[CH]`, `[C]` | Tetravalent carbon |
| N | 3 bonds | `[NH2]`, `[NH]`, `[N]` | Trivalent neutral nitrogen |
| N | 4 bonds | `[NH3(+)]`, `[NH2(+)]`, `[NH(+)]`, `[N(+)]` | Tetravalent cationic nitrogen |
| O | 2 bonds | `[OH]`, `[O]` | Divalent oxygen |
| O | 1 bond | `[O(-)]` | Monovalent anionic oxygen |
| S | 2/4/6 bonds | `[SH]`, `[S]` | Multivalent sulfur |
| Cl, Br, I, F | 1 bond | `[Cl]`, `[Br]`, `[I]`, `[F]` | Monovalent halogens |

**Aromatic Carbon Label Rules**:

**Definition**: Aromatic carbon refers to carbon atoms in aromatic ring systems (e.g., benzene, naphthalene, pyridine rings)

**Rules**:
1. **Aromatic carbon with external substituent**: Must use `[C]` (zero H)
   - External substituent: Any non-ring atom or group (–OH, –Cl, –NH2, –CH3, –NO2, –COOH, etc.)
   - Chemical principle: Substituent occupies one bond position, remaining H count is zero

2. **Aromatic carbon with only ring connections**: Must use `[CH]` (one H)
   - Ring connections: Only bonded to other ring atoms (2 σ bonds)
   - Chemical principle: sp² hybridization, one H perpendicular to ring plane

3. **Bridgehead carbon or fused ring shared carbon**: Use `[C]` (zero H)
   - Chemical principle: Connected to 3 or more ring atoms, no H

**Error Example and Correction**:
```
# Error: Phenol C1 connected to OH but labeled [CH]
phenol_C_1[CH]          # ❌ Should be [C]

# Correct:
phenol_C_1[C]           # ✓ Has OH substituent, use [C]
```

---

### 4. Bond Type Mapping Table (Parser Mapping)

**Syntax Definition**:
```bnf
<connection>  ::= <atom-id> <bond-symbol> <atom-id>
<bond-symbol> ::= "---" | "===" | "-.-" | "-->" | <stereo-bond>
<stereo-bond> ::= "===|" <stereo-label> "|"
<stereo-label> ::= "E" | "Z" | "cis" | "trans"
```

**RDKit Bond Type Mapping**:

| Mermaid Symbol | RDKit.BondType | Chemical Meaning | Notes |
|----------------|----------------|------------------|-------|
| `---` | `Chem.BondType.SINGLE` | σ bond (single) | Most common |
| `===` | `Chem.BondType.DOUBLE` | σ+π bond (double) | Alkenes, carbonyls |
| `-.-` | `Chem.BondType.TRIPLE` | σ+2π bonds (triple) | Alkynes, nitriles |
| `-->` | `Chem.BondType.DATIVE` | Coordinate bond | **Only for complexes** |
| `===\|E\|` | `DOUBLE` + `BondStereo.STEREOE` | E-configuration (trans) | Stereochemistry |
| `===\|Z\|` | `DOUBLE` + `BondStereo.STEREOZ` | Z-configuration (cis) | Stereochemistry |

**Constraints**:
1. Atom IDs in bond definitions must be declared in atom definitions
2. Duplicate bond definitions for the same atom pair are not allowed
3. Coordinate bond `-->` only for metal complexes, forbidden for organic functional groups

**Aromaticity Representation Rules**:

**Principle**: All aromatic rings must use **Kekulé form** (explicit alternating single/double bonds)

**Valid Example** (benzene ring):
```mermaid
benzene_C_1 === benzene_C_2
benzene_C_2 --- benzene_C_3
benzene_C_3 === benzene_C_4
benzene_C_4 --- benzene_C_5
benzene_C_5 === benzene_C_6
benzene_C_6 --- benzene_C_1
```

**Invalid Examples**:
```mermaid
# ❌ Aromatic bond symbols not supported
benzene_C_1 <--> benzene_C_2    # Parser cannot recognize
benzene_C_1:benzene_C_2         # Not valid syntax
```

**Parser Behavior**:
- RDKit automatically recognizes aromaticity during `SanitizeMol()` stage
- Kekulé form is only for graph representation, does not affect final aromaticity recognition

---

### 5. Stereochemistry Encoding Specification

**Double Bond Stereochemistry (E/Z)**:

**Syntax**:
```bnf
<stereo-edge> ::= <atom-id> "===|" ("E"|"Z"|"cis"|"trans") "|" <atom-id>
```

**Parsing Rules** :
```python
stereo_bond_pattern = r'([\w_]+)\s*===\|([EZez]|cis|trans|CIS|TRANS)\|\s*([\w_]+)'
# Capture groups:
#   group(1): atom1_id
#   group(2): stereo_type (case-insensitive)
#   group(3): atom2_id
```

**Examples**:
```mermaid
# (E)-2-butene
E2butene_C_2 ===|E| E2butene_C_3

# (Z)-2-butene
Z2butene_C_2 ===|Z| Z2butene_C_3
```

**Chiral Centers (R/S)**:

**Syntax**: Encoded via atom ID suffix
```bnf
<chiral-atom-id> ::= <subgraph-prefix> "_" <element> "_" <number> "_" ("R"|"S")
```

**Parsing Rules** :
```python
atom_pattern = r'([\w_]+?)(?:_(R|S))?\[([^\]]+)\]'
# Chirality identifier extracted during atom definition
# R → Chem.ChiralType.CHI_TETRAHEDRAL_CW
# S → Chem.ChiralType.CHI_TETRAHEDRAL_CCW
```

**Examples**:
```mermaid
# (R)-2-butanol
R2butanol_C_2_R[CH]     # R-configuration chiral carbon

# (S)-2-butanol
S2butanol_C_2_S[CH]     # S-configuration chiral carbon
```

---

### 6. Subgraph Definition Specification

**Syntax**:
```bnf
<subgraph>     ::= "subgraph" <subgraph-id> "[" <display-name> "]" <subgraph-body> "end"
<subgraph-id>  ::= [A-Za-z][A-Za-z0-9_]*
<display-name> ::= <any Unicode characters>
<subgraph-body> ::= (<comment> | <atom-node> | <connection>)*
```

**Constraints**:
1. `<subgraph-id>` must follow identifier naming rules (ASCII letter start)
2. `<subgraph-id>` must be unique within the document
3. All atom IDs within a subgraph should have prefix matching `<subgraph-id>` (recommended but not enforced)
4. Each subgraph must have corresponding `end` keyword

**Example**:
```mermaid
subgraph chlorophenol["para-chlorophenol"]
    chlorophenol_C_1[C]
    chlorophenol_O_1[OH]
    chlorophenol_C_1 --- chlorophenol_O_1
end
```

---

### 7. Comment Specification (Strict)

**Syntax**:
```bnf
<comment> ::= "%%" <any text> <newline>
```


**Valid Examples**:
```mermaid
%% This is a hydroxyl oxygen
chlorophenol_O_1[OH]

%% Define C-O bond below
chlorophenol_C_1 --- chlorophenol_O_1
```


**Parser Behavior** :
```python
# Parser skips the following lines:
if (line.startswith('%%') or
    line.startswith('graph ') or
    line.startswith('subgraph ') or
    line == 'end' or
    not line):
    continue
```

---

## Special Functional Group Encoding Specification

This section defines standard encoding rules for common organic functional groups to ensure the parser can correctly reconstruct molecular structures.

### 1. Nitro Group (–NO₂) - Mandatory Specification

**Chemical Structure**: R–N⁺(=O)–O⁻

**Encoding Rules**:
1. Nitrogen atom must be labeled `[N(+)]` (positive charge)
2. One oxygen must be labeled `[O]` (neutral, double bond)
3. Another oxygen must be labeled `[O(-)]` (negative charge, single bond)
4. C–N bond must use single bond `---`
5. N=O bond must use double bond `===`
6. N–O⁻ bond must use single bond `---`
7. **Coordinate bonds forbidden** `-->`

**Standard Template**:
```mermaid
graph TB
    subgraph nitromethane["nitromethane"]
        %% Methyl group
        nitromethane_C_1[CH3]

        %% Nitro group three atoms
        nitromethane_N_1[N(+)]
        nitromethane_O_1[O]
        nitromethane_O_2[O(-)]

        %% Bond definitions
        nitromethane_C_1 --- nitromethane_N_1
        nitromethane_N_1 === nitromethane_O_1
        nitromethane_N_1 --- nitromethane_O_2
    end
```

**Bonding Capacity Verification**:
- Nitrogen: 3 bonds (1 to C, 1 double to O, 1 single to O) + positive charge = 3 ✓
- Double-bonded oxygen: 1 double bond = 2 bonds ✓
- Single-bonded oxygen: 1 single bond + negative charge = 2 bonds ✓

**Common Errors**:
```mermaid
# ❌ Error 1: Using coordinate bond
nitromethane_N_1 --> nitromethane_O_1  # Forbidden

# ❌ Error 2: Incorrect charge notation
nitromethane_N_1[N]      # Should be [N(+)]
nitromethane_O_2[O]      # Should be [O(-)]

# ❌ Error 3: Wrong bond type
nitromethane_C_1 === nitromethane_N_1  # C-N should be single bond ---
```

---

### 2. Carboxyl Group (–COOH)

**Chemical Structure**: R–C(=O)–OH

**Encoding Rules**:
1. Carbonyl carbon connected to two oxygen atoms
2. Carbonyl oxygen labeled `[O]` (double bond)
3. Hydroxyl oxygen labeled `[OH]` (single bond)
4. C=O uses double bond `===`
5. C–OH uses single bond `---`

**Standard Template**:
```mermaid
graph TB
    subgraph aceticacid["acetic acid"]
        %% Methyl group
        aceticacid_C_1[CH3]

        %% Carboxyl group
        aceticacid_C_2[C]
        aceticacid_O_1[O]
        aceticacid_O_2[OH]

        %% Bond definitions
        aceticacid_C_1 --- aceticacid_C_2
        aceticacid_C_2 === aceticacid_O_1
        aceticacid_C_2 --- aceticacid_O_2
    end
```

---

### 3. Amines (–NHₓ)

**Encoding Rules**: Nitrogen hydrogen count depends on substitution degree

| Type | Label | Connections | Example |
|------|-------|-------------|---------|
| Primary amine | `[NH2]` | 1 C bond | CH₃NH₂ |
| Secondary amine | `[NH]` | 2 C bonds | (CH₃)₂NH |
| Tertiary amine | `[N]` | 3 C bonds | (CH₃)₃N |

**Example (Methylamine)**:
```mermaid
graph TB
    subgraph methylamine["methylamine"]
        methylamine_C_1[CH3]
        methylamine_N_1[NH2]

        methylamine_C_1 --- methylamine_N_1
    end
```

---

### Functional Group Encoding Reference Table

| Functional Group | Structure | Atom Labels | Bond Types | Notes |
|------------------|-----------|-------------|------------|-------|
| Nitro | –NO₂ | `[N(+)]`, `[O]`, `[O(-)]` | `---`, `===`, `---` | No coordinate bonds, charges in parentheses |
| Carboxyl | –COOH | `[C]`, `[O]`, `[OH]` | `===`, `---` | Carbonyl + hydroxyl |
| Carbonyl | C=O | `[C]`, `[O]` | `===` | Aldehyde/ketone |
| Hydroxyl | –OH | `[OH]` | `---` | Alcohol/phenol |
| Ether | –O– | `[O]` | `---`, `---` | Single bonds on both sides |
| Primary amine | –NH₂ | `[NH2]` | `---` | Monosubstituted |
| Secondary amine | >NH | `[NH]` | `---`, `---` | Disubstituted |
| Tertiary amine | >N– | `[N]` | `---` ×3 | Trisubstituted |
| Quaternary ammonium | >N⁺– | `[N(+)]` | `---` ×4 | Tetrasubstituted, charged |
| Nitrile | –C≡N | `[C]`, `[N]` | `-.-` | Triple bond |
| Halogen | –X | `[F/Cl/Br/I]` | `---` | Single bond |

---

## Quick Reference for Display Labels

| Atom Type | Label | Description |
|-----------|-------|-------------|
| **Carbon** | `[CH3]` | Methyl (1 C bond) |
| | `[CH2]` | Methylene (2 C bonds) |
| | `[CH]` | Methine (3 C bonds) or unsubstituted aromatic C |
| | `[C]` | Quaternary C (4 C bonds) or substituted aromatic C |
| **Nitrogen** | `[NH2]` | Primary amine |
| | `[NH]` | Secondary amine / unsubstituted saturated ring N |
| | `[N]` | Tertiary amine / aromatic N / substituted saturated ring N |
| | `[N(+)]` | Ammonium / nitro nitrogen (positively charged) |
| | `[NH3(+)]` | Protonated primary amine |
| **Oxygen** | `[OH]` | Hydroxyl (alcohol/phenol) |
| | `[O]` | Carbonyl / ether / nitro double-bonded O |
| | `[O(-)]` | Nitro single-bonded O (negatively charged) |
| **Halogens** | `[Cl]` `[Br]` `[I]` `[F]` | Halogen substituents |

---

## Parser Requirements Summary

### Must Satisfy Constraints

#### 1. Structural Constraints
- [x] First line must be `graph TB` or `graph LR`
- [x] Each subgraph must have `subgraph ID["name"]` start and `end` finish
- [x] Atom ID format: `prefix_Element_Number[_Chirality]`
- [x] Atom ID prefix must start with letter, alphanumeric only
- [x] Element numbering starts from 1, counted independently per element type

#### 2. Label Constraints
- [x] Display label format: `[Element][H-digits][(charge)]`
- [x] Aromatic C with substituent must use `[C]`, without must use `[CH]`
- [x] Nitro group must use `[N(+)]`, `[O]`, `[O(-)]`
- [x] Charges in parentheses: `(+)`, `(-)`, `(+digit)`, `(-digit)`

#### 3. Bond Constraints
- [x] Single bond `---`, double bond `===`, triple bond `-.-`
- [x] Coordinate bond `-->` only for complexes, forbidden for organic functional groups
- [x] Aromatic rings must use Kekulé form (alternating single/double bonds)
- [x] Stereo double bonds use `===|E|` or `===|Z|`

#### 4. Syntax Constraints
- [x] Atom definitions before bond definitions
- [x] Bonds must reference defined atom IDs

#### 5. Chemical Constraints
- [x] All atoms follow bonding capacity rules (C=4, N=3/4, O=2, Cl=1)
- [x] Hydrogen count matches connection count
- [x] Stereochemistry correctly encoded (E/Z via edge labels, R/S via ID suffix)

---

## Common Errors Reference Table

| Error | ❌ Wrong | ✅ Correct |
|-------|---------|-----------|
| Aromatic C H-count | Phenol C1 has OH but writes `[CH]` | `chlorophenol_C_1[C]` |
| Aromatic bond symbol | `C_1 <--> C_2` | `C_1 === C_2` |
| Nitro coordinate bond | `N_1 --> O_1` | `N_1 === O_1` |
| Nitro charge missing | `N_1[N]`, `O_2[O]` | `N_1[N(+)]`, `O_2[O(-)]` |
| Old charge format | `N_1[N+]`, `O_2[O-]` | `N_1[N(+)]`, `O_2[O(-)]` |

---

## Verification Checklist

Before generating graphs, verify:

**Structural Integrity**:
- [ ] First line `graph TB/LR`
- [ ] All subgraphs have `subgraph` and `end`
- [ ] All atom IDs follow naming convention
- [ ] All bonds reference defined atoms

**Label Correctness**:
- [ ] Aromatic C H-count correct (substituted `[C]`, unsubstituted `[CH]`)
- [ ] Nitro charges correct (`[N(+)]`, `[O(-)]`, in parentheses)
- [ ] Other charges explicitly in parentheses: `(+)`, `(-)`, `(+digit)`, `(-digit)`
- [ ] H-count follows bonding capacity

**Syntax Compliance**:
- [ ] Bond types correctly mapped
- [ ] Kekulé form aromatic rings
- [ ] Stereochemistry correctly encoded

**Chemical Reasonability**:
- [ ] Bonding capacity satisfied (C=4, N=3/4, O=2)
- [ ] No isolated atoms (except monoatomic molecules)
- [ ] Ring structures closed
- [ ] Functional groups encoded properly
'''