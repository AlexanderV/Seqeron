---
type: concept
title: "K-mer positions (occurrence / position index)"
tags: [analysis, algorithm]
mcp_tools:
  - kmer_positions
sources:
  - docs/Evidence/KMER-POSITIONS-001-Evidence.md
  - docs/algorithms/K-mer/K-mer_Positions.md
source_commit: f43e0daf5dd06eba936e6e14939946e4cd980b67
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: kmer-positions-001-evidence
      evidence: "Test Unit ID: KMER-POSITIONS-001; Algorithm: K-mer Positions (find all start positions where a k-mer occurs in a sequence) (docs/algorithms/K-mer/K-mer_Positions.md)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:both-strand-kmer-counting
      source: kmer-positions-001-evidence
      evidence: "FindKmerPositions is a sibling KmerAnalyzer method; its ASSUMPTION 2 adopts the case-insensitive upper-casing convention of the sibling CountKmers methods (KMER-POSITIONS-001 Evidence)"
      confidence: high
      status: current
---

# K-mer positions (occurrence / position index)

The **K-mer** family's *where*, not *how many*: `KmerAnalyzer.FindKmerPositions(sequence, kmer)`
returns the ascending, **0-based** list of every starting index at which a given k-mer (a fixed
pattern) occurs in a sequence — a **position / occurrence index**, distinct from counting (which
reports *how many times* each k-mer occurs). It solves the classical **Pattern Matching Problem** —
"find all occurrences of a pattern in a string" — exactly (not heuristically). Validated under test
unit **KMER-POSITIONS-001** (record [[kmer-positions-001-evidence]]); [[test-unit-registry]] tracks
the unit and [[algorithm-validation-evidence]] describes the evidence-artifact pattern.

## Core model

Given a text `T` (the sequence, length `L`) and a pattern `P` (the k-mer, length `k`):

> **Occ(P, T) = { i ∈ [0, L − k] : T[i .. i+k) = P }**

reported in **ascending order** with **0-based** indexing (Rosalind BA1D). There are at most
`L − k + 1` candidate start positions (Wikipedia, "k-mer"). Occurrences may **overlap**, and every
overlapping start is included — e.g. `AA` in `AAAA` yields `0, 1, 2`.

## Positions vs counting (why this is its own unit)

This is the *inverse index* to the counting operations, not a special case of them. Counting
([[both-strand-kmer-counting]], [[asynchronous-kmer-counting]] and the sync `CountKmers` primitive)
answers *how many times* every observed k-mer occurs (a `Dictionary<string,int>`); k-mer positions
answers *at which offsets* one specified k-mer occurs (an ordered `IEnumerable<int>`). The count of
returned positions equals that k-mer's overlapping occurrence count (INV-03), so positions strictly
subsume the single-key count while also localizing it — the primitive behind motif scanning, repeat
localization, primer/probe placement, and read-mapping.

It is the **single-pattern K-mer-family sibling** of the multi-pattern
[[known-motif-search]] (`GenomicAnalyzer.FindMotif`): both report **all overlapping** 0-based
ascending occurrences by exact matching; known-motif-search generalizes to a *set* of query motifs
returning per-motif position lists, while `FindKmerPositions` locates one k-mer and returns one list.

## Properties and invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | every returned position `p` satisfies `T[p .. p+k) = P` (case-folded) | direct from the matching predicate |
| INV-02 | positions are 0-based and strictly ascending | single left-to-right scan; 0-based per spec |
| INV-03 | count of returned positions = overlapping occurrence count of `P` in `T` | every candidate index tested; overlaps reported |
| INV-04 | all positions lie in `[0, L − k]`; empty when `k > L` | loop bound `i ≤ L − k`; `L − k + 1` candidates |

## Contract and edge cases

- **Inputs:** `sequence` and `kmer` (both `string`, **case-insensitive** — upper-cased internally via
  `ToUpperInvariant`, mirroring the sibling `CountKmers`). Alphabet is **unrestricted text** — no
  DNA/RNA validation; any characters are matched literally after case-folding.
- **Output:** `IEnumerable<int>` — ascending 0-based start positions of every (possibly overlapping)
  occurrence, produced lazily by `yield return`.
- **Edge cases:** null/empty `sequence` or `kmer` → empty; `|kmer| > |sequence|` (`L − k + 1 ≤ 0`) →
  empty; k-mer absent → empty; k-mer equals the whole sequence → `[0]`. **No exception** is thrown for
  any of these (repository convention, no spec mandate).

## Worked oracles

- `FindKmerPositions("GATATATGCATATACTT", "ATAT")` → `[1, 3, 9]` — Rosalind BA1D sample; `ATAT`
  matches at i=1, 3, 9.
- `FindKmerPositions("AAAA", "AA")` → `[0, 1, 2]` — self-overlap; every overlapping start reported
  (`L − k + 1 = 3`).
- `AGAT` 2-mers / `ana` in `banana` → ascending starts (`AG`@0, `GA`@1, `AT`@2; `ana`→`1, 3`).

## Complexity and implementation

`O(L·k)` time (naive forward scan over `L − k + 1` windows, each compared in `O(k)` via
`ReadOnlySpan<char>.SequenceEqual` to avoid per-window substring allocation), `O(L)` space (the
upper-cased input copies). A **suffix tree was evaluated and rejected**: repository
`SuffixTree.FindAllOccurrences` counts overlapping occurrences correctly but returns positions in
**unordered** leaf-collection order (this unit requires ascending) and only amortizes its `O(n)`
construction when *many* patterns are queried against one preprocessed text — not repaid for a single
k-mer query. Correctness (overlapping, 0-based) is unaffected by the choice.

## Deviations and assumptions

- **No source contradictions; deviations = None.** Rosalind BA1D and Wikipedia agree on the
  occurrence-set definition, the overlapping rule, and the `L − k + 1` bound.
- **Three ASSUMPTIONS**, all repository-interoperability / API-shape (non-correctness-affecting on the
  all-uppercase evidence examples): **0-based** indexing (per the machine-checked BA1D exercise, over
  the textbook's 1-based prose), **case-insensitive** matching (sibling-`CountKmers` upper-casing), and
  **null/empty → empty** (no exception).

## Relation to other units

- Shares the case-insensitive upper-casing convention of the sibling K-mer counting methods
  ([[both-strand-kmer-counting]], [[asynchronous-kmer-counting]], and the pending sync `CountKmers`
  primitive) but is a distinct operation (index of offsets, not a count table).
- Sibling of the abundance-vector K-mer units [[k-mer-euclidean-distance]] and
  [[kmer-jaccard-similarity]], and of the space-enumeration [[k-mer-generation]] (which produces the
  full `n^k` universe rather than locating a specific k-mer in a sequence).
- The single-pattern counterpart of the multi-pattern exact-matcher [[known-motif-search]], sharing
  the exact all-overlapping-occurrences semantics also seen in [[longest-common-substring]] and
  [[dot-plot-word-match]].
