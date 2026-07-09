---
type: concept
title: "Codon optimization (OptimizeSequence)"
tags: [annotation, algorithm]
sources:
  - docs/Evidence/CODON-OPT-001-Evidence.md
source_commit: 8d1b85e321fa52d6dea20205e2d0a4d2f28d1dbc
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: codon-opt-001-evidence
      evidence: "ID: CODON-OPT-001 ... Algorithm: Sequence Optimization (OptimizeSequence)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:codon-adaptation-index
      source: codon-opt-001-evidence
      evidence: "MaximizeCAI: Use most frequent codons for each amino acid — Source: Sharp & Li (1987), CAI definition; optimization drives the CAI score toward its maximum"
      confidence: high
      status: current
---

# Codon optimization (OptimizeSequence)

The one **rewriting** operation in the codon-usage family: given a protein-coding sequence and a
target-organism codon-usage table, **replace codons with synonymous alternatives** to improve
heterologous (cross-species) expression, while leaving the encoded protein unchanged. Where
[[relative-synonymous-codon-usage|RSCU]], [[codon-adaptation-index|CAI]] and
[[effective-number-of-codons|ENC/Nc]] *measure* a gene's existing codon bias, codon optimization
*produces a new gene* that better matches the host's tRNA abundances. Validated as
[[codon-opt-001-evidence|CODON-OPT-001]]; see [[test-unit-registry]] for how the unit is tracked.

## Why it exists

Rare codons in a host slow translation, can deplete charged tRNAs and stall ribosomes, and — because
translation rate influences **cotranslational folding** — can also perturb the folded product. Codon
optimization tunes codon choice to the host to raise yield. The literature (Plotkin & Kudla 2011)
catalogues several strategies beyond "use the most frequent codon": local mRNA folding, codon-pair
bias, the 5' **codon ramp**, and **codon harmonization**.

## Optimization strategies

Seqeron's `OptimizeSequence` exposes five strategies, each tracing to a source point in the artifact:

- **MaximizeCAI** — use the single most-frequent synonymous codon for each amino acid, driving the
  [[codon-adaptation-index|CAI]] toward 1 (Sharp & Li 1987).
- **BalancedOptimization** — balance CAI against a **GC-content constraint** (40–60% target), because
  mRNA secondary structure and GC affect translation. After GC balancing it **rebuilds the `Changes`
  list** so every modification is reflected.
- **HarmonizeExpression** — match the *distribution* of host codon usage rather than always picking the
  optimum, i.e. codon harmonization (Mignon et al. 2018).
- **AvoidRareCodons** — replace only codons whose frequency falls below a threshold; leaves the rest
  untouched (consecutive rare codons inhibit translation). This is the actuator for exactly what
  [[rare-codon-analysis]] detects.
- **MinimizeSecondary** — avoid mRNA secondary structures (5' structure especially inhibits
  initiation). In codon *selection* it delegates to BalancedOptimization; a dedicated
  `ReduceSecondaryStructure` method handles structure reduction separately.

## Invariants

- **Protein preservation.** Optimization must not change the encoded amino-acid sequence — every
  substitution is synonymous. This is the defining correctness invariant, verified across *all*
  strategies.
- **Met and Trp are fixed points.** AUG (Met) and UGG (Trp) are the only codons for their amino acids,
  so they are never rewritten — the same single-codon-family fact that gives them `w ≡ 1` in
  [[codon-adaptation-index|CAI]].
- **Stop codons preserved**, not optimized.
- **CAI range** of any result stays in **(0, 1]** (geometric mean of `w ∈ (0, 1]`).

## Implementation and edge behaviour

- **RNA notation.** Works in RNA (`U`, not `T`); DNA input is auto-converted `T → U`; processing is
  case-insensitive.
- **Codon framing.** Trims to complete codons (`length % 3 == 0`); empty input → empty result.
- **Zero-frequency clamp.** The internal CAI reporting clamps zero-frequency codons to `1e-6` per the
  Sharp & Li prescription — the same guard documented on [[codon-adaptation-index]].
- **Failure modes.** An unknown codon translates to `X` (or errors); a stop codon mid-sequence can
  truncate the protein prematurely; non-RNA characters are handled gracefully.

Because different organisms have different preferred codons, the *same* input optimizes differently
per host — the artifact carries E. coli K12 (Kazusa 316407), S. cerevisiae (4932) and H. sapiens
(9606) preferred-codon tables (e.g. Leu → CUG in E. coli/human but UUA/UUG in yeast; Arg → CGC/CGU in
E. coli but AGA in yeast) as the organism-specific test fixtures.

## Place in the codon-usage family

Codon optimization is the **actuator** of the family that RSCU/CAI/ENC only measure: it consumes a
codon-usage table (the object [[relative-synonymous-codon-usage|RSCU]] normalizes) and optimizes toward
a high [[codon-adaptation-index|CAI]]. Its **AvoidRareCodons** strategy removes exactly the codons
[[rare-codon-analysis]] flags. Other siblings still in `docs/Evidence/` (raw frequency/usage tables)
share the same synonymous-family machinery. See
[[algorithm-validation-evidence]] for the shared evidence-artifact pattern.
</content>
</invoke>
