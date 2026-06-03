# 🧬 MoleCode

### An LLM-native, graph-explicit molecular language

*MoleCode presents molecules as code so LLMs can read, write, edit, and reason on chemistry directly — instead of reconstructing structure from cryptic strings.*

[Paper (arXiv:2605.16480)](https://arxiv.org/abs/2605.16480) · [GitHub](https://github.com/AtomFlow-AI/MoleCode) · [AtomFlow](https://atomflow-ai.com/)

<img src="https://raw.githubusercontent.com/AtomFlow-AI/MoleCode/main/docs/assets/overview.png" alt="MoleCode overview" width="100%">

---

## What is MoleCode?

A molecule **is** a graph: atoms are nodes, bonds are edges. Yet LLMs are usually
fed molecules as *linear strings* (SMILES) where the graph is **implicit** — so
the model must rebuild the structure from syntax before it can reason.

**MoleCode makes the structure the language.** Every atom and bond is a typed
declaration with a persistent identifier, serialized as a
[Mermaid](https://mermaid.js.org/) graph, and **deterministically, losslessly
inter-convertible with SMILES / MOL via RDKit** (no learned model). One grammar
spans **small molecules, polymers, and Markush structures**.

## Install

```bash
pip install molecode
```

## Quick start

```python
from rdkit import Chem
from molecode import mol_to_mermaid, mermaid_to_mol, mol_to_smiles

# SMILES -> MoleCode graph
graph = mol_to_mermaid(Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O"), name="Aspirin")
print(graph)

# MoleCode graph -> SMILES (lossless round-trip)
assert mol_to_smiles(mermaid_to_mol(graph)) == Chem.CanonSmiles("CC(=O)Oc1ccccc1C(=O)O")
```

Three domains, one grammar:

```python
from molecode.polymer import polymer_to_mermaid, mermaid_to_psmiles   # polymers (×n repeat unit)
from molecode.markush import MoleCodeGraph, molecode_isomorphic        # Markush ({R1}/{Boc} nodes)
from molecode.prompts import MOLECULE_SYSTEM_PROMPT, MARKUSH_SYSTEM_PROMPT  # LLM system prompts
from molecode import LLMClient                                         # optional OpenAI-compatible client
```

## Why it matters

- **Generalization, not memorization.** On novel molecules SMILES accuracy falls
  to ~20% while MoleCode holds ~76–80% (same model, same knowledge).
- **Cheaper reasoning.** Sub-linear chain-of-thought growth — ~5× lower total
  token cost per query than SMILES.
- **Scales to big, repetitive objects.** Full-chain SMILES collapses toward 0% as
  polymer chains grow; MoleCode stays flat.

## For coding agents

MoleCode also ships as an **Agent Skill** (Claude Code auto-discovers it; Codex
reads `AGENTS.md`) with a CLI for SMILES/PSMILES/Markush ↔ MoleCode conversion,
validation, isomorphism comparison, and OCSR (molecule image → MoleCode). See the
[repository](https://github.com/AtomFlow-AI/MoleCode).

## Citation

```bibtex
@article{yan2026molecode,
  title={MoleCode unlocks structural intelligence in large language models},
  author={Yan, Zhiyuan and Liu, Chen and Zhao, Boxuan and Lin, Kaiqing and Zhao, Jixiang and Wang, Yimi and Lv, Liuzhenghao and Li, Hao and Zhang, Shanzhuo and Yuan, Li and others},
  journal={arXiv preprint arXiv:2605.16480},
  year={2026}
}
```

## License

[MIT](https://github.com/AtomFlow-AI/MoleCode/blob/main/LICENSE) © 2026 AtomFlow-AI
