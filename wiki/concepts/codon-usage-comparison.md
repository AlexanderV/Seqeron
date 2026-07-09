---
type: concept
title: "Codon usage table and TVD comparison"
tags: [annotation, algorithm]
sources:
  - docs/Evidence/CODON-USAGE-001-Evidence.md
source_commit: dda05cc3a1a2b473ced5663ee7beeeb062b11907
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: codon-usage-001-evidence
      evidence: "Test Unit ID: CODON-USAGE-001 ... Methods Under Test: CalculateCodonUsage, CompareCodonUsage"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:relative-synonymous-codon-usage
      source: codon-usage-001-evidence
      evidence: "CalculateCodonUsage splits the sequence into codons and counts occurrences (raw counts) — the same non-overlapping-triplet tally RSCU normalizes"
      confidence: high
      status: current
---

# Codon usage table and TVD comparison

The **raw** end of the codon-usage family (`CodonOptimizer.CalculateCodonUsage` +
`CompareCodonUsage`, test unit [[codon-usage-001-evidence|CODON-USAGE-001]]): produce an
unnormalized codon-count table, and compare two such tables with a distribution-distance
similarity. Distinct from the codon-usage-*bias* measures — it does no per-amino-acid
normalization and computes no single per-gene score. See [[test-unit-registry]] for tracking.

## The raw usage table

`CalculateCodonUsage(sequence)` splits a coding sequence into non-overlapping triplets from
offset 0, counts each codon, and returns `Dictionary<string, int>` of **raw counts**:

    Count(c) = |{ i : sequence[i : i+3] = c }|

It converts `T → U` internally (RNA representation), uppercases, and drops an incomplete trailing
codon. This is the same counting primitive behind [[relative-synonymous-codon-usage|RSCU]]'s
`CountCodons`, but **left unnormalized** — hence the raw base of the family rather than a bias
measure. Invariant: Σ counts = total codons.

## TVD-based comparison

`CompareCodonUsage(seq1, seq2)` normalizes both count tables to frequency distributions
(`f(c) = count(c) / total`) and returns a **Total Variation Distance (TVD) similarity** in
`[0, 1]`:

    Similarity = 1 − ( Σ_c |f₁(c) − f₂(c)| ) / 2

The `Σ|f₁−f₂|` is the L¹ distance between the two frequency vectors; halving it gives the TVD of
the two probability distributions, and `1 − TVD` turns distance into similarity.

**Proven properties** (from TVD theory, and used as differential-test oracles):

- **Identity** — `sim(s, s) = 1.0` (zero distance for identical distributions).
- **Symmetry** — `sim(a, b) = sim(b, a)` (since `|x−y| = |y−x|`).
- **Range** — `sim ∈ [0, 1]` (TVD of probability distributions ∈ [0, 1]).
- **Disjoint → 0** — non-overlapping codon sets give `Σ|f₁−f₂| = 1 + 1 = 2`, so `sim = 0`.
- **Partial overlap** — exact value derivable for any input (e.g. 2/3 shared codons → `sim = 2/3`).

**Edge behaviour:** empty input → similarity `0` (convention: no data → 0), not NaN or an
exception; `T`/`U` treated as equivalent; incomplete trailing codons ignored.

## Why TVD (a deliberate metric choice)

Wikipedia's codon-usage-bias article lists cosine similarity and correlation coefficients among
codon-table comparison metrics; this unit uses TVD-based similarity instead. Every expected test
value is analytically derivable from the TVD formula, and the four proven properties above follow
directly from TVD theory — so the choice is defensible and fully oracle-backed rather than a
departure from a fixed reference. The predefined organism tables (E. coli K12, S. cerevisiae,
H. sapiens) were verified against the **Kazusa Codon Usage Database** (all 64 relative fractions
match to 2 dp).

## Relation to other codon-usage statistics

This is the unnormalized base of the family. **[[relative-synonymous-codon-usage|RSCU]]** normalizes
the same codon counts per synonymous family; **[[codon-adaptation-index|CAI]]** reduces them to a
single geometric-mean gene score; **[[effective-number-of-codons|ENC/Nc]]** measures reference-free
bias via codon homozygosity; **[[rare-codon-analysis]]** thresholds per-family frequencies; and
**[[codon-optimization]]** consumes a usage table to rewrite a CDS. The aggregation/reporting view
is [[codon-stats-001-evidence|CODON-STATS-001]] (`GetStatistics`). What is unique here is the
**distribution-comparison** operation (`CompareCodonUsage`): a TVD similarity between two whole
codon-frequency distributions, which none of the bias measures provide.
