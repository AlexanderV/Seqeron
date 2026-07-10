---
type: concept
title: "Approximate pattern matching with mismatches (Hamming ball / d-neighborhood)"
tags: [motif, algorithm]
sources:
  - docs/Evidence/PAT-APPROX-003-Evidence.md
source_commit: dd94a7819ba9764ca0165e710710b83844931da9
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: pat-approx-003-evidence
      evidence: "Test Unit ID: PAT-APPROX-003; Algorithm: Best Match and Frequency Analysis (Approximate Pattern Matching, Frequent Words with Mismatches)"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:k-mer-positions
      source: pat-approx-003-evidence
      evidence: "Both locate a pattern's occurrences in a text by a forward window scan; k-mer positions requires EXACT equality (Hamming distance 0), while approximate matching accepts any window with HammingDistance(window, pattern) <= d (ROSALIND BA1H). d=0 makes them coincide (Count_0 = exact occurrence counting)."
      confidence: high
      status: current
---

# Approximate pattern matching with mismatches (Hamming ball / d-neighborhood)

**Approximate pattern matching** relaxes exact substring matching to allow up to `d`
**substitution mismatches** (no indels — a fixed-length, **Hamming-distance** model). A
length-`k` pattern *"appears with at most d mismatches"* at a window of the text when the
window and the pattern differ in at most `d` positions. This is the *k-mismatch* end of
pattern matching, distinct from the exact matchers ([[k-mer-positions]],
[[known-motif-search]]) and from the edit-distance / alignment family (which also allows
insertions and deletions). It is the Compeau & Pevzner *Bioinformatics Algorithms* ch. 1
motif-hunting toolkit — the way biology finds a degenerate signal (e.g. the *DnaA* box)
that never occurs as an exact literal. Validated under test unit **PAT-APPROX-003**;
the validation record is [[pat-approx-003-evidence]], [[test-unit-registry]] tracks the
unit, and [[algorithm-validation-evidence]] describes the artifact pattern.

## The core definition (Hamming ball, not edit distance)

Let `HammingDistance(u, v)` be the number of positions at which equal-length strings `u`
and `v` differ (unit **PAT-APPROX-001**). A length-`k` pattern `P` **appears as an
approximate occurrence** at start `i` of text `T` when

> `HammingDistance(P, T[i .. i+k)) <= d`   (ROSALIND BA1H)

Everything in the family is built from this predicate:

| Operation | Repository method | Returns |
|-----------|-------------------|---------|
| Approximate occurrences | `FindApproximateOccurrences` | all 0-based starts `i` with `HammingDistance(P, T[i..i+k)) <= d` |
| Count with mismatches (`Count_d`) | `CountApproximateOccurrences` | the number of such windows |
| Frequent words with mismatches | `FindFrequentKmersWithMismatches` | the most-frequent `k`-mer(s) counted over their `d`-neighborhood |
| d-neighborhood | `Neighbors` | the set of all `k`-mers within Hamming distance `d` of `P` |
| Best approximate match | `FindBestMatch` | the minimum-Hamming-distance equal-length window (leftmost) |

## Two counting semantics — occurrences vs neighborhood tally

There are **two different "counts"** in this family and conflating them is the classic
defect:

- **`Count_d(Text, Pattern)`** counts *windows of Text* within Hamming distance `d` of a
  **fixed** pattern (BA1H/BA1I). It is a scan: for each `k`-window, test the predicate.
  O(n·m).
- **Frequent words with mismatches** (BA1I) asks *which `k`-mer (over the whole 4^k
  space, not just substrings) has the largest `Count_d`*. The efficient implementation
  tallies **over each window's `d`-neighborhood**: for every window, enumerate
  `Neighbors(window, d)` and increment each neighbor's count (the reference
  `go-rosalind` "for kmer and neighbors" tally). The winner may be a `k`-mer that
  **never occurs exactly** in the text.

**Key consequence — the pattern need not be a substring.** `AAAAA` is the most frequent
5-mer with 1 mismatch in `AACAAGCTGATAAACATTTAAAGAG` even though `AAAAA` never appears
exactly (`Count_1 = 4`, over windows `AACAA, ATAAA, AAACA, AAAGA`).

## The d-neighborhood

`Neighbors(Pattern, d)` = the set of all `k`-mers whose Hamming distance from `Pattern`
does not exceed `d` (BA1N) — the **Hamming ball** of radius `d`. It is built by
recursive position-choice enumeration and **always contains `Pattern` itself** (the
identity, distance 0), so a window's own `k`-mer is always counted in the frequent-words
tally. For a 4-letter alphabet the size is a sum of binomial terms; for `k=3, d=1` it is
`1 + 3·(4−1) = 10` (e.g. `Neighbors(ACG, 1)` = {ACG, CCG, TCG, GCG, AAG, ATG, AGG, ACA,
ACC, ACT}). This neighbor enumeration is why the approach is **practical only for small
k and d** — the textbook bound is `k <= 12, d <= 3`; it grows combinatorially beyond
that.

## d = 0 degenerates to exact matching

When `d = 0` the whole family collapses to its exact counterparts: `Count_0` is exact
occurrence counting, `FindApproximateOccurrences` reduces to exact
[[k-mer-positions]], and `Neighbors(Pattern, 0) = {Pattern}`. This is the `alternative_to`
link — approximate and exact matching are the same operation at different mismatch
tolerances.

## Worked oracles

- **Frequent words with mismatches** (BA1I): `Text = ACGTTGCATGTCGCATGATGCATGAGAGCT`,
  `k = 4`, `d = 1` → **{GATG, ATGC, ATGT}**, each with max count **5**. *All* ties are
  returned.
- **Approximate occurrences** (BA1H): `Pattern = ATTCTGGA`,
  `Text = CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC`,
  `d = 3` → positions **{6, 7, 26, 27, 78}**, so `Count_3 = 5`.
- **Count_1 worked example**: `CountApproximateOccurrences(AACAAGCTGATAAACATTTAAAGAG,
  AAAAA, 1) = 4`.
- **d-neighborhood** (BA1N): `Neighbors(ACG, 1)` → 10 members (listed above).

## Corner cases and contract

- **Pattern absent as exact substring** — counting is over the Hamming ball, not exact
  occurrences (BA1H/BA1I).
- **Multiple ties (BA1I)** — when several `k`-mers share the max count, **all** are
  returned.
- **Neighborhood includes identity** — `Neighbors(P, d)` always contains `P`.
- **`FindBestMatch` tie-break (documented API convention).** No textbook problem defines
  a single "best approximate match". `FindBestMatch(sequence, pattern)` returns the
  equal-length window with **minimum** Hamming distance, scanning left-to-right and
  keeping the **first** window that achieves a strictly smaller distance (leftmost
  minimum-distance window; an exact match short-circuits with distance 0). The returned
  *distance value* is fully evidence-defined by `HammingDistance`; only the leftmost
  tie-break among equal-minimum windows is a repository convention, not a
  correctness-affecting parameter.

## Scope and relation to other units

This is the **substitution-only (Hamming)** approximate matcher. It does **not** model
insertions/deletions — for indel-tolerant matching use the edit-distance / alignment
family ([[global-alignment-needleman-wunsch]], [[semi-global-alignment-fitting]]). It
builds directly on the `HammingDistance` primitive (**PAT-APPROX-001**, BA1H
definition). It is the mismatch-tolerant `alternative_to` the exact single-pattern
locator [[k-mer-positions]] and shares the "surface over-represented `k`-mers" goal with
the exact de-novo discovery unit [[overrepresented-kmer-discovery]] — but frequent-words-
with-mismatches ranks `k`-mers by their **Hamming-ball** count rather than exact
occurrences, which is what lets it recover degenerate biological motifs.

## Reference sources

**ROSALIND BA1I** (frequent words with mismatches, `Count_d` definition + worked
example), **BA1H** (approximate-occurrence definition + positions sample), **BA1N**
(d-neighborhood definition + `ACG`/1 sample) — all Compeau & Pevzner, *Bioinformatics
Algorithms*, ch. 1. Reference implementations **charlesreid1/go-rosalind**
(`rosalind_ba1.go`, the neighbor-tally frequent-words approach + O(n·m) approximate
match) and **zonghui0228/Rosalind-Solutions** (`rosalind_ba1h.py`, O(n·m) confirmation).
**No deviations**; the only assumption is the `FindBestMatch` leftmost-minimum tie-break,
a documented API convention that does not change the returned minimum distance.
