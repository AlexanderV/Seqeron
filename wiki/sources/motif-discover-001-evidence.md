---
type: source
title: "Evidence: MOTIF-DISCOVER-001 (De novo motif discovery via over-represented k-mers, O/E enrichment)"
tags: [validation, motif]
doc_path: docs/Evidence/MOTIF-DISCOVER-001-Evidence.md
sources:
  - docs/Evidence/MOTIF-DISCOVER-001-Evidence.md
source_commit: ca2ab6940b31b0df620d2b47f5d84f4cc49652d5
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: MOTIF-DISCOVER-001

The validation-evidence artifact for test unit **MOTIF-DISCOVER-001** — **de novo motif
discovery via over-represented k-mers** (`MotifFinder.DiscoverMotifs`): enumerate every
length-`k` k-mer of one DNA sequence, count occurrences, and rank by the observed/expected
(O/E) enrichment ratio. It is one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the model, contract,
invariants, worked oracles, and corner cases are synthesized in
[[overrepresented-kmer-discovery]]. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **Compeau & Pevzner, *Bioinformatics Algorithms*** (rank 1, textbook content reproduced
    verbatim on the wikiselev GitHub wiki) — the **zero-order i.i.d. uniform** background
    ("each nucleotide with probability 0.25"), the probability formula
    `Pr(N,A,Pattern,t) ≈ ( N−t·(k−1) | t )/A^(t·k)`, and the **expected count** at `t=1`:
    `(N − k + 1)/4^k` (window count `N−k+1` over `4^k` distinct k-mers). Worked example
    `Pr(1000,4,9,1)·500 ≈ 1.9` (a random 9-mer expected ≈2 times across 500 length-1000
    sequences).
  - **O/E-ratio corroboration** (rank 3, WebSearch → monaLisa `getKmerFreq` + PeerJ
    supplemental) — the standard overrepresentation statistic is the **ratio of observed to
    expected** k-mer frequency (enrichment > 1 ⇒ over-represented), with
    **expected = N/4^k** under the zero-order uniform model (consistent with `(N−k+1)/4^k`
    once the exact window count is used).
- **Datasets (documented oracles):**
  - *Tandem repeat* `ATGCATGCATGC` (N=12, k=4): windows 9, `E = 9/256 = 0.03515625`,
    `ATGC` at positions **0,4,8** (Count 3) → enrichment **768/9 ≈ 85.333**.
  - *Homopolymer* `AAAAAAAAAA` (N=10, k=3): windows 8, `E = 8/64 = 0.125`, `AAA` at
    positions **0..7** (Count 8) → enrichment **64.0**.
- **Corner cases / failure modes:**
  - *Self-overlap approximation* — the closed-form *probability* statistic warns it ignores
    self-overlap; this affects only the probability, **not** the deterministic Count or the
    O/E expected-count denominator used here.
  - *`k > N`* — zero length-`k` windows (`N − k + 1 ≤ 0`); discovery returns nothing.

## Deviations and assumptions

**Deviations: none** — the exact expected-count `(N − k + 1)/4^k`, the O/E enrichment ratio,
and overlapping window enumeration are implemented verbatim from the cited theory. Documented
**intentional simplifications** (not deviations from the validated statistic): zero-order
uniform background only (no higher-order Markov model, so O/E can over/under-state enrichment
on compositionally biased sequences), and no closed-form p-value / E-value (rank by Count /
Enrichment or an external significance tool).

One **assumption**: the `minCount` parameter (default 2) is a **presentation threshold**, not
part of the published statistic — the formal O/E value is defined for every k-mer regardless of
the cutoff. Changing `minCount` only includes/excludes rows; it never alters a returned record's
Count, Positions, or Enrichment, so it is not correctness-affecting.

Recommended coverage (MUST): exact O/E for `ATGC` in `ATGCATGCATGC` k=4 = 768/9; exact O/E for
`AAA` in `AAAAAAAAAA` k=3 = 64.0; observed Count + positions of a repeated k-mer (deterministic
window enumeration). SHOULD: `minCount` filter excludes below-threshold k-mers; null → 
`ArgumentNullException`, `k < 1` → `ArgumentOutOfRangeException`. COULD: `k > N` → no motifs.
No contradictions among sources — Compeau & Pevzner and the monaLisa/PeerJ corroboration agree
that expected = length/4^k and enrichment = O/E under the uniform model.
