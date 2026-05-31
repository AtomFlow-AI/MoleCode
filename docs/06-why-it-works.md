# MoleCode — Why it works

> The empirical case for explicit structure, summarized from the MoleCode paper.

MoleCode is more verbose than SMILES on *input* — so why does it win? Because it
**reallocates inference from structural reconstruction to productive reasoning**.
The findings below are reported in the MoleCode technical report; the figures live
in [`assets/`](assets).

## Generalization, not memorization

When molecules are split by familiarity (popular / classical / novel PubChem
compounds), SMILES accuracy falls from **~42% on familiar molecules to ~20% on
novel ones** — a memorization signature. MoleCode holds **~76–80% across all
tiers**. On molecules that are both complex *and* novel, SMILES collapses to ~16%
while MoleCode stays ~72%. A canonical-vs-permuted-SMILES probe shows a +17% gap
for SMILES (it has memorized specific strings) and **no gap** for MoleCode.

![Generalization and reasoning](assets/results_1_main.png)

## Robustness to molecular size

As molecular weight rises from 50–150 to 350–500 Da, SMILES accuracy drops from
~65% to ~32%; MoleCode stays ~78% → ~76%.

## Token economics: long input, short reasoning

MoleCode's *input* scales as ~`38.4·C + 144` tokens versus SMILES' ~`2.0·C + 86`
(C = carbon count) — longer. **But** the chain-of-thought scales **sub-linearly**
for MoleCode (~C^0.52) and **super-linearly** for SMILES (~C^1.65) and SELFIES
(~C^1.36). Net **total** cost per query:

| Representation | Total tokens / query |
| --- | --- |
| **MoleCode** | **~2,005** |
| SELFIES | ~7,178 |
| SMILES | ~10,370 |

About a **5× reduction** — MoleCode operates in a *long-input / short-CoT /
high-accuracy* regime. And under MoleCode, more thinking correlates with better
results (ρ=0.34, p=0.03); under SMILES it does not (ρ=0.03, p=0.87) — extra SMILES
tokens go into reconstructing connectivity in the wrong direction.

![Inference scaling](assets/results_3_scaling.png)

## Goal-directed design induces interpretable edits

On LogP and aqueous-solubility optimization, MoleCode makes **localized,
property-aligned, chemically interpretable** edits while preserving Tanimoto
similarity (e.g. fewer H-bond donors/acceptors and more halogenation to raise
LogP) — exactly the moves a medicinal chemist would make.

![Goal-directed design](assets/results_2_chemistry.png)

## The benefit grows with size and repetition

For polymers, full-chain SMILES accuracy falls toward **0%** as the chain
lengthens, and PSMILES degrades as the repeat unit gets more complex. MoleCode —
explicit repeat unit + symbolic `×n` — stays accurate on **both** axes, and also
wins on polymer editing and *de novo* generation.

![Long molecules](assets/results_4_long_molecules.png)

## One grammar for scientific structure

The same Subgraph–Node–Edge grammar extends beyond atoms to **Markush families**,
**reaction mechanisms**, and **multimodal document parsing** — and beats
domain-specialized formats. Markush understanding rises from **38.1% → 84.0%**;
mechanism prediction improves on accuracy (+16.1%), atom conservation (+31.1%),
electron-transfer correctness (+9.6%), and continuous-path validity (+15.1%).

![General language](assets/results_5_extension.png)

## Limitations

MoleCode does not add knowledge a model lacks — a weak model can still emit an
invalid structure. It is more verbose on input than SMILES. Highly specialized
domains (organometallics, inorganic solids, biological macromolecules) may need
extra primitives for coordination geometry or periodicity.

Back to [01-overview.md](01-overview.md).
