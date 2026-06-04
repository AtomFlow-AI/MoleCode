# MoleCode — API reference

> The public API of the installable `molecode` package (`pip install molecode`).

`pip install molecode` gives you the **library** only. The runnable
[`examples/`](../examples) and the [Agent Skill](../.claude/skills/molecode/) live
in the repository, not in the PyPI package. Dependencies: `rdkit` (all
conversion) and `networkx` (Markush isomorphism), pulled in automatically.

```python
import molecode
molecode.__version__          # '0.1.0'
```

The package has five submodules:

| Module | Purpose |
| --- | --- |
| [`molecode.molecule`](#molecodemolecule) | small molecules ⇄ MoleCode graph |
| [`molecode.polymer`](#molecodepolymer) | polymers (repeat unit + `×n`) ⇄ MoleCode |
| [`molecode.markush`](#molecodemarkush) | Markush `{}` abbreviation nodes + graph isomorphism |
| [`molecode.prompts`](#molecodeprompts) | LLM system prompts |
| [`molecode.llm`](#molecodellm) | optional OpenAI-compatible client |

The top-level `molecode` namespace re-exports the most common small-molecule
helpers and the client:

```python
from molecode import mol_to_mermaid, mermaid_to_mol, mol_to_smiles, mol_to_inchi, LLMClient
```

---

## `molecode.molecule`

Deterministic, lossless conversion between RDKit molecules / SMILES and the
MoleCode graph for small molecules.

**Functions**

```python
mol_to_mermaid(mol: Chem.Mol, name: str = "Molecule", kekulize: bool = True) -> str
mermaid_to_mol(mermaid_text: str, strict: bool = True) -> Optional[Chem.Mol]
mol_to_smiles(mol: Chem.Mol) -> str
mol_to_inchi(mol: Chem.Mol) -> str
has_invalid_atoms(mol: Chem.Mol) -> bool
get_invalid_atom_labels(mermaid_text: str, mol: Chem.Mol) -> List[str]
visualize_mol(mol, title="", save_path=None)
visualize_mols_grid(mols, legends=None, mols_per_row=4, sub_img_size=(300, 300), save_path=None)
```

- `strict=True` makes `mermaid_to_mol` return `None` if the graph contains
  unrecognized/dummy atoms; use `strict=False` to keep them as `*` placeholders.
- `kekulize=False` emits aromatic bonds as `<-->` instead of Kekulé single/double.

**Classes** — `MolToMermaidConverter(subgraph_name="Molecule").convert(mol, kekulize=True)`,
`MermaidMolParser().parse_mermaid_graph(text)` (the functions above wrap these).

```python
from rdkit import Chem
from molecode.molecule import mol_to_mermaid, mermaid_to_mol, mol_to_smiles

graph = mol_to_mermaid(Chem.MolFromSmiles("CCO"), name="Ethanol")
mol_to_smiles(mermaid_to_mol(graph))        # -> 'CCO'
```

```
graph TB
    subgraph Ethanol["Ethanol"]
        Ethanol_C_1[CH3]
        Ethanol_C_2[CH2]
        Ethanol_O_1[OH]
        Ethanol_C_1 --- Ethanol_C_2
        Ethanol_C_2 --- Ethanol_O_1
    end
```

---

## `molecode.polymer`

Polymers as an explicit repeat unit carrying a symbolic `×n` count. R/S and E/Z
stereochemistry round-trip.

**Functions**

```python
polymer_to_mermaid(repeat_smiles: str, n: int, label=None,
                   terminus_left=None, terminus_right=None,
                   name="Polymer", kekulize=True) -> str
block_copolymer_to_mermaid(blocks: List[BlockSpec],
                           terminus_left=None, terminus_right=None,
                           name="BlockCopolymer") -> str
mermaid_to_psmiles(mermaid: str) -> Optional[str]
parse_element(label: str) -> tuple[str, int]
```

The repeat-unit SMILES must contain exactly two `*` attachment points (first `*`
= left/entry, second `*` = right/exit).

**Classes / dataclasses** — `BlockSpec(smiles, n, label=None)`,
`PolymerSpec(blocks, terminus_left=None, terminus_right=None, name="Polymer")`,
`RepeatUnitConverter`, `PolymerToMermaidConverter`.

```python
from molecode.polymer import polymer_to_mermaid, mermaid_to_psmiles, block_copolymer_to_mermaid, BlockSpec

graph = polymer_to_mermaid("*CC(C)*", n=100, name="PP")
mermaid_to_psmiles(graph)                   # -> '*CC(*)C'

block_copolymer_to_mermaid(
    [BlockSpec("*CCO*", n=20, label="PEG"),
     BlockSpec("*CC(C)O*", n=10, label="PPO")],
    name="PEG-b-PPO",
)
```

---

## `molecode.markush`

Markush / generic structures with `{}` abbreviation nodes, plus an RDKit-free
graph-isomorphism comparator.

**Functions**

```python
mermaid_to_mol(mermaid_text, strict=True) -> Optional[Chem.Mol]   # use strict=False to keep {abbrev} as '*'
mol_to_mermaid(mol, name="Molecule") -> str
mol_to_smiles(mol) -> str ; mol_to_inchi(mol) -> str
has_invalid_atoms(mol) -> bool ; get_invalid_atom_labels(text, mol) -> List[str]
normalize_abbrev_name(name: str) -> str
molecode_isomorphic(g1: MoleCodeGraph, g2: MoleCodeGraph,
                    ignore_stereo=True, normalize_abbrevs=True,
                    abbrev_expand_map=None) -> Tuple[bool, dict]
build_expand_map() -> dict
```

**Classes** — `MoleCodeGraph` (`.from_text(text)` classmethod, `.to_networkx()`,
`.get_abbrev_labels()`, `.num_nodes`, `.num_edges`), `NodeInfo`, `EdgeInfo`,
`MolToMermaidConverter`, `MermaidMolParser`.

**Data tables** — `SINGLE_ATOM_MAP` (24), `SUBGRAPH_MAP` (78),
`NON_EXPANDABLE` (120), `EXPAND_MAP` (217; the combined expansion lookup).

```python
from molecode.markush import MoleCodeGraph, molecode_isomorphic, EXPAND_MAP

a = MoleCodeGraph.from_text('graph TB\n subgraph M["m"]\n M_N_1[NH]\n M_X_1{Boc}\n M_N_1 --- M_X_1\n end')
b = MoleCodeGraph.from_text('graph TB\n subgraph M["m"]\n M_X_1{NHBoc}\n end')
molecode_isomorphic(a, b, abbrev_expand_map=EXPAND_MAP)
# -> (True, {'reason': 'isomorphic after expansion', ...})    # [NH]-{Boc} == {NHBoc}
```

---

## `molecode.prompts`

The authoritative grammar specifications to feed an LLM as the system prompt.

```python
from molecode.prompts import MOLECULE_SYSTEM_PROMPT   # small-molecule grammar (~18k chars)
from molecode.prompts import MARKUSH_SYSTEM_PROMPT     # Markush + image→MoleCode / OCSR (~6k chars)
```

---

## `molecode.llm`

An optional, dependency-free, OpenAI-compatible chat client. You supply the API
key and base URL — nothing is hard-coded.

```python
class LLMClient(api_key=None, base_url=None, model=None, *, timeout=120.0, default_temperature=0.0)
    .chat(user, system=None, *, images=None, model=None, temperature=None, **extra) -> str
    .complete(messages, *, model=None, temperature=None, **extra) -> str

call_llm(system, user, *, temperature=0.0, **client_kwargs) -> Optional[str]
image_to_data_uri(path_or_url: str) -> str
DEFAULT_MODEL      # 'gemini-3.1-pro-preview'
DEFAULT_BASE_URL   # 'https://api.openai.com/v1'
```

Credentials fall back to env vars `MOLECODE_API_KEY` / `OPENAI_API_KEY`,
`MOLECODE_BASE_URL`, `MOLECODE_MODEL`. Pass `images=[...]` (file paths or URLs) to
a vision model for OCSR (molecule image → MoleCode).

```python
from molecode import LLMClient
from molecode.prompts import MOLECULE_SYSTEM_PROMPT, MARKUSH_SYSTEM_PROMPT

client = LLMClient(api_key="sk-...", base_url="https://api.openai.com/v1", model="gpt-4o")

# text task
client.chat("How many carbons are in this molecule? ...", system=MOLECULE_SYSTEM_PROMPT)

# OCSR: molecule image -> MoleCode (needs a vision-capable model)
client.chat("Output this molecule as a MoleCode graph.",
            system=MARKUSH_SYSTEM_PROMPT, images=["molecule.png"])
```

---

## Task cheat-sheet

| Goal | How |
| --- | --- |
| Let an LLM read / edit a molecule | `mol_to_mermaid` → give `MOLECULE_SYSTEM_PROMPT` → `mermaid_to_mol` to validate the reply |
| Polymers | `polymer_to_mermaid` ⇄ `mermaid_to_psmiles` |
| Markush authoring / scoring | hand-write a `{R}` graph → `MoleCodeGraph` + `molecode_isomorphic` |
| Image → structure (OCSR) | `LLMClient.chat(..., system=MARKUSH_SYSTEM_PROMPT, images=[...])` |

Runnable end-to-end scripts for all of these are in [`../examples/`](../examples)
(`01`–`08`), available when you clone the repository.
