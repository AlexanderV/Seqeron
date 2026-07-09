---
type: concept
title: "k-mer Jaccard similarity (alignment-free sequence resemblance)"
tags: [analysis, algorithm]
sources:
  - docs/Evidence/GENOMIC-SIMILARITY-001-Evidence.md
  - docs/algorithms/Analysis/Sequence_Similarity.md
source_commit: f2b9ce29b93a0977bf8cc2d4d003a59711a6534b
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: genomic-similarity-001-evidence
      evidence: "Test Unit ID: GENOMIC-SIMILARITY-001; Algorithm: Sequence Similarity (k-mer Jaccard index)"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:alignment-statistics
      source: genomic-similarity-001-evidence
      evidence: "Sequence_Similarity.md §2.5 contrasts k-mer Jaccard (set resemblance of k-mer composition, O(n+m), order-insensitive) with alignment identity (residue-by-residue over an alignment, O(n·m), positional) as two ways to measure sequence similarity"
      confidence: high
      status: current
---

# k-mer Jaccard similarity (alignment-free sequence resemblance)

The **Analysis** family's alignment-free pairwise-similarity measure: `GenomicAnalyzer.CalculateSimilarity`
scores how alike two DNA sequences are by the **Jaccard index of their distinct k-mer sets**, reported
as a percentage in `[0, 100]`. Each sequence is reduced to the set of its distinct length-`k` substrings;
the Jaccard index is the fraction of k-mers shared out of all distinct k-mers in either set. The result
is **exact** — full k-mer sets, no MinHash sketching — and computed in `O(n+m)`, making it a cheap
resemblance estimate that needs no alignment. Validated under test unit **GENOMIC-SIMILARITY-001**
(record [[genomic-similarity-001-evidence]]); [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the evidence-artifact pattern.

## Core model (Jaccard 1901)

For two finite sets A and B:

```
J(A,B) = |A ∩ B| / |A ∪ B| = |A ∩ B| / (|A| + |B| − |A ∩ B|),   0 ≤ J ≤ 1
```

Here **A and B are the sets of distinct k-mers** of the two sequences. Applied to k-mer sets J is
"the fraction of shared k-mers out of all distinct k-mers" (Ondov et al. 2016, *Mash*) — the same
alignment-free k-mer-Jaccard resemblance Mash popularized for genome/metagenome distance (Mash adds
MinHash *sketching*; this unit computes the exact index over full sets). The method reports `J × 100`.

## Properties and invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `0 ≤ result ≤ 100` | `0 ≤ J ≤ 1`, scaled ×100 |
| INV-02 | Identical non-empty sequences → 100.0 | A = B ⇒ A∩B = A∪B ⇒ J = 1 |
| INV-03 | Disjoint k-mer sets → 0.0 | A∩B = ∅ over a non-empty union |
| INV-04 | Symmetric: `result(a,b,k) = result(b,a,k)` | ∩ and ∪ are commutative |
| INV-05 | k-mers compared **as sets** — within-sequence repeats counted once | Jaccard over *distinct* k-mers (a `HashSet<string>`, not a multiset) |

## Contract and edge cases

- Inputs: two `DnaSequence` (non-null, A/C/G/T normalized/validated at the `DnaSequence` boundary → comparison is case-insensitive) and `kmerSize` (default **5**, `≥ 1`).
- `null` sequence → `ArgumentNullException`; `kmerSize < 1` → `ArgumentOutOfRangeException`.
- **Empty-union convention returns `0.0`** — when both sequences are empty (or both shorter than k) both k-mer sets are empty and the union is empty. Jaccard is mathematically *undefined* for an empty union ("not well-defined when μ(A ∪ B) = 0"); the implementation returns `0.0` as "no shared content". This is a documented implementation convention (ASM-1), not a literature-mandated value — either 0 or 1 appears in practice.
- One sequence empty, the other non-empty → `0.0` (empty intersection over a non-empty union).

## Worked oracles (hand-derived, k=3)

| seq1 | seq2 | \|A∩B\| | \|A∪B\| | result |
|------|------|---------|---------|--------|
| `ACGTACGT` | `ACGTACGA` | 4 | 5 | **80.0** (partial overlap) |
| `ACGT` | `ACGA` | 1 | 3 | **100/3 ≈ 33.33** (non-integer fraction) |
| `ACGTACGT` | `ACGTACGT` | 4 | 4 | **100.0** (identical) |
| `AAAAA` | `CCCCC` | 0 | 2 | **0.0** (disjoint) |

Distinct-set semantics also give `AAAAAA` vs `AAAA` at k=3 → **100** (both sets = `{AAA}`).

## Deviations and assumptions

Two source-backed **assumptions**, neither a correctness gap:

- **Empty-union → 0.0** (ASM-1) — the undefined-case convention above; documented and tested as the implementation contract, not asserted as the literature value.
- **×100 percentage scaling** — the formal index is in `[0,1]`; multiplying by 100 is a presentation convention that does not change relative ordering. The **default k = 5** is a project resolution choice (Mash uses k=21 for whole genomes); k only sets comparison resolution, never the formula.

**Not implemented:** MinHash sketching / approximate Jaccard (`j(A_s,B_s)=|A_s∩B_s|/s`) — exact sets suffice at the supported sizes; and Jaccard distance `1 − J` — derive as `1 − result/100`. **Deviations: none** beyond these.

## Relation to other similarity measures

- **`alternative_to` [[alignment-statistics]]** — the positional counterpart. Alignment percent identity/similarity is residue-by-residue over an alignment (`O(n·m)`, order-sensitive); k-mer Jaccard is set resemblance of composition (`O(n+m)`, order-*insensitive*: two sequences with identical k-mer composition but different arrangement score identically). Use alignment when position/edit distance matters.
- The same **k-mer Jaccard metric** is the alignment-free similarity that ranks best hits in [[ortholog-detection-reciprocal-best-hits]] (there a 5-mer Jaccard replaces BLAST bit-score), and it is the exact-set basis Mash sketches for [[average-nucleotide-identity|genome distance]].
- **Limitations** — exact-match k-mers only (no mismatch tolerance; one substitution perturbs up to k k-mers); composition-only (insensitive to order/position).

## Reference tools

Definitions trace to **Jaccard, Paul (1901)** (the index `|A∩B|/|A∪B|`, range `[0,1]`, defined for non-empty sets; *Bulletin de la Société vaudoise des sciences naturelles* 37(142):547–579) and **Ondov et al. (2016)** (*Mash*, k-mer-set Jaccard = fraction of shared k-mers, MinHash estimate `|A_s∩B_s|/s`; *Genome Biology* 17:132), with the Mash distance-estimation documentation corroborating the sketch form. No source contradictions — Jaccard's set definition and Mash's k-mer-set application are mutually consistent; the only implementation choices (empty-union → 0, ×100 scaling, default k=5) are documented conventions.
