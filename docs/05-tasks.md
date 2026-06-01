# MoleCode — LLM task families

> How to drive understanding, generation, editing, and reasoning with any LLM.

MoleCode is a drop-in representation for an LLM. The pattern is always the same:

1. Give the model the grammar as a **system prompt** —
   [`molecode.prompts.MOLECULE_SYSTEM_PROMPT`](../molecode/prompts) (or
   `MARKUSH_SYSTEM_PROMPT`).
2. Convert your molecule to a graph with `mol_to_mermaid(...)` and put it in the
   **user message**.
3. Parse and **validate** the model's reply with `mermaid_to_mol(...)` — the
   output *is* a structure, so validation is a deterministic RDKit round-trip,
   not fragile string parsing.

Runnable demos for every family live in [`../examples/`](../examples) (offline by
default; set `MOLECODE_API_KEY` to call a model).

> **Calling a model.** Use any LLM SDK — the prompts are plain strings — or the
> bundled dependency-free, OpenAI-compatible client. You supply the key and URL:
>
> ```python
> from molecode import LLMClient
> client = LLMClient(api_key="sk-...", base_url="https://api.openai.com/v1", model="gpt-4o-mini")
> reply = client.chat(user_prompt, system=MOLECULE_SYSTEM_PROMPT)
> ```
>
> Credentials can instead come from `MOLECODE_API_KEY` / `MOLECODE_BASE_URL` /
> `MOLECODE_MODEL`, so you never commit a key.

## 1. Understanding — [`examples/04_understanding.py`](../examples/04_understanding.py)

Read a structure, answer a structural question. Because every atom is an explicit
node, the model reads structure directly instead of reconstructing it.

Covered tasks: carbon/atom counting, molecular-formula prediction, ring counting,
functional-group identification, ring-system / scaffold identification, structural
equivalence, IUPAC ↔ structure, reaction-product prediction, NMR elucidation.

```text
Task: Count the number of carbon atoms in the molecule below.
Molecule (Mermaid): ```mermaid … ```
Answer:
```

## 2. Generation — [`examples/05_generation.py`](../examples/05_generation.py)

Design a molecule under constraints; the model emits a MoleCode graph that you
parse and validate.

Covered tasks: atom-number-constrained generation, functional-group-guided
generation, and text-based *de novo* design. For polymers, *de novo* repeat-unit
generation (constraint satisfaction and fragment assembly).

```python
from molecode.molecule import mermaid_to_mol, has_invalid_atoms, mol_to_smiles
mol = mermaid_to_mol(model_output, strict=False)
ok  = mol is not None and not has_invalid_atoms(mol)   # deterministic validation
```

## 3. Editing — [`examples/06_editing.py`](../examples/06_editing.py)

Modify a structure: add, delete, or substitute atoms/groups. In MoleCode an edit
is a **local graph operation** — adding a methyl is literally one new node plus
one new edge — so edits are localized, auditable, and diff-able rather than a
whole-string rewrite. (Remember the explicit-H bookkeeping: bonding onto a
`[CH3]` makes it `[CH2]`.)

Covered tasks: atom/group addition, deletion, substitution (small molecules and
polymer repeat units).

## 4. Reasoning — [`examples/07_reasoning.py`](../examples/07_reasoning.py)

Multi-step chemical reasoning, e.g. reaction-product or reaction-mechanism
prediction. Persistent atom IDs let the model **track atom identity across a
transformation**, and `%%` comments let it write its reasoning *inside* the graph.
MoleCode also represents mechanisms as ordered graph states with explicit
electron-transfer edits.

## Task coverage by domain

| Domain | Understanding | Generation | Editing | Reasoning |
| --- | :---: | :---: | :---: | :---: |
| Molecules | ✅ | ✅ | ✅ | ✅ |
| Polymers | ✅ | ✅ | ✅ | — |
| Markush | ✅ | — | — | — |

Next: [06-why-it-works.md](06-why-it-works.md)
