---
type: concept
title: "k-mer Euclidean distance (alignment-free word-frequency distance)"
tags: [analysis, algorithm]
sources:
  - docs/Evidence/KMER-DIST-001-Evidence.md
  - docs/algorithms/K-mer/K-mer_Euclidean_Distance.md
source_commit: e59aee7693f6331696abe9d68f59844f570fcf17
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: kmer-dist-001-evidence
      evidence: "Test Unit ID: KMER-DIST-001; Algorithm: K-mer Euclidean Distance (alignment-free word-frequency distance)"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:kmer-jaccard-similarity
      source: kmer-dist-001-evidence
      evidence: "Both are alignment-free pairwise k-mer measures over the same length-k words, but Euclidean distance operates on normalized frequency VECTORS (Lau 2022 count/(L-k+1)) capturing abundance, whereas Jaccard operates on presence/absence k-mer SETS; Vinga & Almeida 2003 list Euclidean among the vector metrics over the 4^k word-frequency vector"
      confidence: high
      status: current
---

# k-mer Euclidean distance (alignment-free word-frequency distance)

The **Analysis** family's alignment-free pairwise **dissimilarity** measure:
`KmerAnalyzer.KmerDistance` scores how *different* two sequences are by the **Euclidean (L2)
distance between their normalized k-mer frequency vectors**. Each sequence is summarized by the
frequencies of its length-`k` substrings; the distance is the root-sum-of-squares of the per-word
frequency differences over the union of observed k-mers. It is exact, deterministic, needs no
alignment, and runs in `O(n + m)` — used for fast whole-genome phylogeny, clustering, and database
screening where alignment is impractical. Identical sequences give distance 0; more dissimilar
sequences give larger values. Validated under test unit **KMER-DIST-001** (record
[[kmer-dist-001-evidence]]); [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the evidence-artifact pattern.

## Core model (word-frequency vector + L2)

For a sequence `s` of length `L`, the number of overlapping k-mer windows is `L − k + 1`. The
**frequency** vector normalizes each k-mer count by that window total (Lau et al. 2022):

```
f_s(w) = count_s(w) / (L_s − k + 1)
```

Over the union `W` of k-mers occurring in either sequence `x` or `y`, the distance is

```
d(x, y) = sqrt( Σ_{w ∈ W}  ( f_x(w) − f_y(w) )² ),   d ≥ 0
```

where a word absent from a sequence contributes a `0` component (Zielezinski et al. 2017). The
union is a sparse representation of the full `4^k`-dimensional word-frequency vector (Vinga &
Almeida 2003) — unobserved words are `0` on both sides and contribute nothing. Boden et al. 2014
confirm the standard alignment-free Euclidean distance is taken over **relative (normalized)**
frequency vectors, not raw counts.

## Properties and invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `d(x, x) = 0` | Equal frequency vectors → zero sum of squares; "identical sequences yield a distance of 0" |
| INV-02 | `d(x, y) = d(y, x)` | `(f_x − f_y)² = (f_y − f_x)²`; Euclidean distance is a metric |
| INV-03 | `d(x, y) ≥ 0` | Square root of a sum of squares |
| INV-04 | Two disjoint single-k-mer sequences (each frequency 1) → `√2` | Vectors `(1,0)` and `(0,1)` ⇒ `√((1−0)²+(0−1)²)` |

## Contract and edge cases

- Inputs: two sequences (`string`) and `k` (`> 0`); no alphabet restriction — any character may
  form a k-mer.
- `k ≤ 0` → `ArgumentOutOfRangeException` (inherited from `CountKmers`).
- **Case-insensitive** (ASM-01) — inputs are upper-cased before counting, so mixed-case inputs
  produce identical k-mers.
- **Empty / too-short input (ASM-02)** — a null/empty sequence, or one shorter than `k`, yields an
  empty frequency vector treated as the **zero vector**: the distance then equals the L2 norm of the
  other sequence's frequency vector, and `0` when both are empty. (Sources define frequencies only
  for `L ≥ k`; this is the natural extension.)
- Entry points: `KmerDistance(seq1, seq2, k)`; `GetKmerFrequencies(string, int)` (count ÷ sum of
  counts = window total); `CountKmers(string, int)` (the underlying counter).

## Worked oracles (hand-derived)

| seq1 | seq2 | k | distance |
|------|------|---|----------|
| `ATGTGTG` | `CATGTG` | 3 | **√0.11 ≈ 0.33166247903553997** (Zielezinski 2017 Fig. 1: f_x=(0.2,0,0.4,0.4), f_y=(0.25,0.25,0.25,0.25)) |
| `AAAA` | `AAAT` | 1 | **√0.125 ≈ 0.3535533906** (f₁=(1.0,0), f₂=(0.75,0.25)) |
| any `s` | same `s` | any | **0** (INV-01) |
| disjoint single-k-mer | — | — | **√2** (INV-04) |

Note the count-based (un-normalized) form gives `√3 ≈ 1.732…` on the Fig. 1 example; `KmerDistance`
uses the **frequency** form (√0.11) per Lau et al. 2022.

## Deviations and assumptions

Two source-backed **assumptions**, neither a correctness gap: case folding (ASM-01) and empty →
zero-vector (ASM-02), both above. **Deviations: none** — the implementation follows Zielezinski's
word-vector model, Lau's `L − k + 1` normalization, and Boden's relative-frequency Euclidean
distance exactly.

**Not implemented:** count-based Euclidean and the other vector metrics (Manhattan, Canberra,
Chebyshev, cosine, D2 — Vinga & Almeida 2003 list several); spaced-word / spaced-k-mer frequencies
(Boden et al. 2014). Only the contiguous-k-mer frequency Euclidean variant is provided.

## Relation to other alignment-free measures

- **`alternative_to` [[kmer-jaccard-similarity]]** — the other alignment-free pairwise k-mer measure
  in this library. Both reduce sequences to their length-`k` words, but Jaccard is a **set
  resemblance** (presence/absence of distinct k-mers, `|A∩B|/|A∪B|`, → *similarity*), while this is
  an **L2 distance over frequency vectors** (captures k-mer *abundance*, → *dissimilarity*). Jaccard
  is insensitive to how often a k-mer occurs; the Euclidean distance is not. Use Jaccard for
  composition overlap; use the frequency distance when relative k-mer abundance matters.
- Both build on the same single-linear-pass k-mer counting primitive as
  [[both-strand-kmer-counting]] and [[asynchronous-kmer-counting]] (`CountKmers`); a suffix tree is
  evaluated-not-used (counting + a vector difference, not substring search).
- **Limitations** — frequency normalization removes overall length, so two sequences with identical
  relative composition but different lengths can have distance 0; contiguous k-mers only (no spaced
  words, no reverse-complement strand merge, no D2*/D2S background correction). The raw value is a
  dissimilarity, not a calibrated phylogenetic distance.

## Reference tools

Definitions trace to **Zielezinski et al. (2017)** (*Genome Biology* 18:186 — word-vector model,
Fig. 1 worked example, "very commonly computed by the Euclidean distance", `d(x,x)=0`),
**Lau et al. (2022)** (*NAR Genomics and Bioinformatics* — frequency = count ÷ `(L − k + 1)`,
Euclidean among the usable metrics), **Vinga & Almeida (2003)** (*Bioinformatics* 19(4):513 — the
`4^k` word-frequency vector and its family of dissimilarity scores), and **Boden et al. (2014)**
(*Bioinformatics* 30(14):1991 — Euclidean over relative-frequency vectors). No source
contradictions.
</content>
