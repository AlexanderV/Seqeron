---
type: concept
title: "Unique & min-count k-mers (frequency-1 singletons + Count ≥ t frequency filtering)"
tags: [analysis, algorithm]
mcp_tools:
  - kmers_with_min_count
  - unique_kmers
sources:
  - docs/Evidence/KMER-UNIQUE-001-Evidence.md
  - docs/algorithms/K-mer/Unique_And_MinCount_Kmers.md
source_commit: 8107afde5814b7194214921227878bd429f3fc04
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: kmer-unique-001-evidence
      evidence: "Test Unit ID: KMER-UNIQUE-001; Algorithm: Unique K-mers / K-mers with Minimum Count (docs/algorithms/K-mer/Unique_And_MinCount_Kmers.md)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:k-mer-statistics
      source: kmer-unique-001-evidence
      evidence: "BioInfoLogics distinguishes 'unique' (frequency exactly 1) from 'distinct' (each different k-mer counted once); FindUniqueKmers returns the count==1 subset, whereas KmerStatistics.UniqueKmers is the distinct count — the same ATCGATCAC k=3 example gives 5 unique vs 6 distinct"
      confidence: high
      status: current
---

# Unique & min-count k-mers (frequency-1 singletons + Count ≥ t frequency filtering)

The **K-mer** family's *frequency-filtering* operations: `KmerAnalyzer.FindUniqueKmers` and
`KmerAnalyzer.FindKmersWithMinCount`. Both operate on the per-k-mer **Count** (the number of
overlapping occurrences of a k-mer as a substring of the sequence — Compeau & Pevzner
`Count(Text, Pattern)`), and select k-mers at opposite ends of the frequency distribution:

- **`FindUniqueKmers(sequence, k)`** — the k-mers occurring **exactly once** (`Count == 1`), i.e.
  the *singletons*.
- **`FindKmersWithMinCount(sequence, k, minCount)`** — the k-mers with `Count ≥ minCount`
  (recurrent / over-represented k-mers), returned with their counts **ordered by Count descending**.

Validated under test unit **KMER-UNIQUE-001** (record [[kmer-unique-001-evidence]]);
[[test-unit-registry]] tracks the unit and [[algorithm-validation-evidence]] describes the
evidence-artifact pattern.

## The total / distinct / unique terminology

For a sequence `S` of length `L` and window length `k`, the three standard k-mer counts
(BioInfoLogics; Wikipedia) are:

| Term | Definition | For `ATCGATCAC`, k=3 |
|------|------------|----------------------|
| **Total** | all `L − k + 1` overlapping windows, duplicates included | 7 |
| **Distinct** | the number of *different* k-mer strings present | 6 |
| **Unique** | the k-mers whose frequency is **exactly 1** (singletons) | 5 → {TCG, CGA, GAT, TCA, CAC} |

`ATC` occurs twice (positions 0 and 4), so it is distinct but **not** unique. This is the canonical
definition `FindUniqueKmers` implements.

## Gotcha: "unique" here = singletons, ≠ `KmerStatistics.UniqueKmers`

The word "unique" is overloaded in the K-mer family. **This** unit's `FindUniqueKmers` returns the
`Count == 1` subset (the 5 singletons above). By contrast, [[k-mer-statistics]]'s
`KmerStatistics.UniqueKmers` field holds the **distinct** count `D` (each different k-mer counted
once — 6 for the same input). They coincide only when every distinct k-mer is a singleton
(no repeats). Keep the two straight: *singleton set* here vs *distinct count* there.

## Min-count filtering (`FindKmersWithMinCount`)

Selecting k-mers whose `Count ≥ t` is the standard way to isolate **recurrent / over-represented**
k-mers (Compeau & Pevzner: a *most-frequent* k-mer maximises `Count`; thresholding at `t` generalises
this). The result carries each qualifying k-mer with its count, **sorted by Count descending**.

## Worked oracles

| Sequence | k | Call | Result |
|----------|---|------|--------|
| `ATCGATCAC` | 3 | `FindUniqueKmers` | {TCG, CGA, GAT, TCA, CAC} (5; ATC=2 excluded) |
| `AGAT` | 2 | `FindUniqueKmers` | {AG, GA, AT} (all 3; every 2-mer Count 1) |
| `ACGTACGT` | 4 | `FindUniqueKmers` | {CGTA, GTAC, TACG} (ACGT=2 excluded) |
| `ACGTACGT` | 4 | `FindKmersWithMinCount(…, 2)` | {(ACGT, 2)} |
| `ACGTACGT` | 4 | `FindKmersWithMinCount(…, 1)` | all 4 distinct, Count-desc, ACGT (2) first |
| `AAAAA` | 3 | `FindUniqueKmers` | ∅ (single distinct 3-mer AAA has Count 3) |

For `ACGTACGT` k=4 the 5 total 4-mers are ACGT, CGTA, GTAC, TACG, ACGT → counts ACGT=2, others=1.

## Contract and edge cases

- **Inputs:** `sequence` (string; upper-cased internally, so matching is **case-insensitive**),
  `k` (`int`, `> 0`), and for the min-count method `minCount` (`int`).
- **`k ≤ 0`** → `ArgumentOutOfRangeException` (shared with the sibling counting methods).
- **`k > L` or empty/null sequence** → empty result for both methods (`L − k + 1 ≤ 0` ⇒ no k-mers).
- **Output:** `FindUniqueKmers` → a set/collection of the singleton k-mers; `FindKmersWithMinCount`
  → k-mer→count pairs ordered by Count descending.

## Deviations and assumptions

- **No source contradictions; deviations = None.** The Wikipedia `L − k + 1` total, BioInfoLogics
  distinct/unique distinction, and Compeau & Pevzner `Count(Text, Pattern)` are implemented verbatim.
- **ASSUMPTION (`minCount ≤ 1`):** the predicate `Count ≥ minCount` is then satisfied by every
  observed k-mer, so the method returns **all distinct** k-mers ordered by count — the mathematically
  consistent extension of `Count ≥ t`, treated as defined behaviour (not an invented value).
- **ASSUMPTION (case normalisation):** input is upper-cased so case variants count as the same
  k-mer (consistent with the sibling `KmerAnalyzer` methods); no source contradicts this.

## Relation to other units

- **Frequency filters over the shared count.** Both methods consume the same synchronous
  `KmerAnalyzer.CountKmers` sliding-window multiset that underlies the rest of the family
  ([[asynchronous-kmer-counting]], [[both-strand-kmer-counting]]), then filter it by `Count`.
- **Opposite ends of the same distribution:** `FindUniqueKmers` selects the frequency-1 tail;
  `FindKmersWithMinCount` selects the `Count ≥ t` head (recurrent k-mers) — the counterpart to the
  "most-frequent k-mers" query owned by [[k-mer-search]]. That same `FindUniqueKmers` singleton method
  is also documented by the K-mer Search spec (KMER-FIND-001) alongside the most-frequent and `(L, t)`
  clump operations synthesized on [[k-mer-search]].
- Distinct from [[k-mer-statistics]] (summarises the whole count profile, and whose `UniqueKmers`
  field is the *distinct* count — see the gotcha above), [[k-mer-generation]] (enumerates the full
  `n^k` universe, sequence-independent), [[k-mer-positions]] (indexes **where** one k-mer occurs),
  and [[k-mer-euclidean-distance]] / [[kmer-jaccard-similarity]] (compare **two** sequences' profiles).
