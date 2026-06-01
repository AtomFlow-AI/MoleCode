# MoleCode examples

**English** | [中文](README.zh-CN.md)

Runnable demonstrations of MoleCode. Run any of them from the repository root:

```bash
pip install -e .
python examples/01_molecule_roundtrip.py
```

| Script | What it shows | Needs an LLM? |
| --- | --- | --- |
| `01_molecule_roundtrip.py` | SMILES → MoleCode graph → SMILES (lossless) | no |
| `02_polymer_roundtrip.py` | PSMILES → MoleCode graph (`×n`) → PSMILES; block copolymers | no |
| `03_markush_roundtrip.py` | Markush `{}` abbreviation nodes; abbreviation-aware isomorphism | no |
| `04_understanding.py` | **Understanding** — count atoms / formula / rings | optional |
| `05_generation.py` | **Generation** — de novo design under constraints | optional |
| `06_editing.py` | **Editing** — local graph edits (add/delete/substitute) | optional |
| `07_reasoning.py` | **Reasoning** — reaction-product prediction | optional |
| `08_image_to_molecode.py` | **OCSR** — molecule image → MoleCode graph | vision model |

`08_image_to_molecode.py` performs OCSR (optical chemical structure recognition):
it renders a demo molecule to a PNG (or takes an image path you pass), sends the
image plus the Markush MoleCode system prompt to a **vision-capable** model, and
parses the returned MoleCode graph back to a structure. Set `MOLECODE_MODEL` to a
vision model (e.g. `gpt-4o-mini`, `gpt-4o`); with no API key it dry-runs.

## The four LLM task families

Examples 04–07 each build the exact prompt you would send to any LLM and pass it
through `examples/_llm.py` (a thin wrapper over the library's
`molecode.llm.LLMClient`). They run **offline by default** — with no API key set
they print the assembled system + user prompt (a "dry run"), so you can see
precisely what MoleCode sends. To actually call a model, set environment
variables (any OpenAI-compatible endpoint works):

```bash
export MOLECODE_API_KEY="sk-..."
export MOLECODE_BASE_URL="https://api.openai.com/v1"   # or your provider
export MOLECODE_MODEL="gpt-4o-mini"                    # or any chat model
python examples/04_understanding.py
```

In your own code, use the client directly (api key + url are yours to supply):

```python
from molecode import LLMClient
client = LLMClient(api_key="sk-...", base_url="https://api.openai.com/v1", model="gpt-4o-mini")
reply = client.chat(user_prompt, system=system_prompt)
```

MoleCode itself never calls an LLM — the client is an optional, dependency-free
convenience. The reusable ingredients are:

* `molecode.prompts.MOLECULE_SYSTEM_PROMPT` / `MARKUSH_SYSTEM_PROMPT` — the
  grammar specification you give the model as a system prompt;
* `molecode.molecule.mol_to_mermaid(...)` — turn your molecule into the graph
  the model reads;
* `molecode.molecule.mermaid_to_mol(...)` — parse and **validate** the model's
  output deterministically with RDKit.
