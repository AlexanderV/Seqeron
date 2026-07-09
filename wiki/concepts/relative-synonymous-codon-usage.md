---
type: concept
title: "Relative Synonymous Codon Usage (RSCU)"
tags: [annotation, algorithm]
sources:
  - docs/Evidence/ANNOT-CODONUSAGE-001-Evidence.md
source_commit: 1a7037743423e6db365bbc5460f2d2d04f9384a5
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: annot-codonusage-001-evidence
      evidence: "Test Unit ID: ANNOT-CODONUSAGE-001 ... Algorithm: Relative Synonymous Codon Usage (RSCU)"
      confidence: high
      status: current
---

# Relative Synonymous Codon Usage (RSCU)

A per-codon **codon-usage-bias** measure: within each amino acid's set of synonymous codons,
how over- or under-represented one codon is relative to uniform (unbiased) usage. Introduced by
**Sharp & Li (1986)**, it is the canonical normalization for comparing codon preference across
genes, genomes, and organisms independent of amino-acid composition. Validated as
[[annot-codonusage-001-evidence|ANNOT-CODONUSAGE-001]] against the LIRMM methods page,
PMC2528880, and the CodonU reference implementation. See [[test-unit-registry]] for how the
unit is tracked.

## The measure

For codon `j` of amino acid `i`, with `n_i` codons in that amino acid's synonymous family and
`x_{i,j}` occurrences of codon `j`:

    RSCU_{i,j} = n_i · x_{i,j} / Σ_j x_{i,j}

The denominator sums occurrences over **all** synonymous codons of amino acid `i`.
Equivalently, RSCU is the observed codon frequency divided by the frequency expected if all
synonymous codons were used equally.

**Reading the value:** RSCU = 1.0 ⇒ no bias (used exactly as often as expected); > 1.0 ⇒
preferred / over-represented; < 1.0 ⇒ under-represented. Values are bounded in **[0, n_i]**.

**Invariant:** the RSCU values over an observed synonymous family sum to `n_i` (an algebraic
identity of the formula, and a good differential-test oracle).

## Edge behaviour

- **Aggregation.** Codon counts are **pooled across all input reference sequences** before the
  ratio is computed — RSCU is a property of the aggregate codon pool, not an average of
  per-sequence values.
- **Sense codons only.** Stop codons (TAA, TAG, TGA in the standard code) are excluded from the
  families (Biopython `forward_table` convention).
- **Single-codon amino acids.** Met (ATG) and Trp (TGG) have `n_i = 1`, so their RSCU is always
  exactly **1.0** regardless of count.
- **Unobserved family.** If a whole synonymous family has zero occurrences the denominator is
  0 and the base definition is undefined; Seqeron reports **0.0** for each member.
- **Genetic code.** Defaults to NCBI Standard table 1; an overload accepts a non-standard
  `GeneticCode` for alternate tables (the six-codon families are Leu, Arg, Ser).

## Relation to other codon-usage statistics

RSCU is the shared normalization at the base of the codon-usage-bias family. The **[[codon-adaptation-index|CAI]]**
(Codon Adaptation Index) reuses RSCU-style relative adaptiveness — but referenced to the family
**maximum** codon and reduced to a single geometric-mean gene score — and guards `log(0)` for
unobserved codons with a pseudocount-style adjustment that plain RSCU does **not** apply (Sharp &
Li 1987 use a **0.5 pseudocount** at reference-table build time; Seqeron's CAI implementation uses
a **`1e-6` clamp** at score time — see [[codon-adaptation-index]]). The
**[[effective-number-of-codons|ENC / Nc]]** (Wright 1990) is the reference-free sibling: it reduces a
gene's codon bias to a single number in [20, 61] using codon homozygosity `F̂` built from the same
synonymous-codon frequencies `p_i`. **[[codon-optimization]]** is the family's *rewriting* operation:
it consumes a codon-usage table and substitutes synonymous codons to improve host expression. Other
siblings in `docs/Evidence/` include rare-codon analysis and raw codon-frequency/usage tables; each
builds on the same synonymous-family counting that RSCU formalizes.
