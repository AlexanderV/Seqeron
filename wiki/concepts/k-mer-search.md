---
type: concept
title: "K-mer search (special-interest k-mers: most-frequent words, singletons, (L,t) clumps)"
tags: [analysis, algorithm]
sources:
  - docs/algorithms/K-mer/K-mer_Search.md
source_commit: 1a7f80318ba001f378cbe1931974115a95515c13
created: 2026-07-13
updated: 2026-07-15
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: kmer-find-001-report
      evidence: "Test Unit ID: KMER-FIND-001; Algorithm: K-mer Search (docs/algorithms/K-mer/K-mer_Search.md), KmerAnalyzer.FindMostFrequentKmers / FindUniqueKmers / FindClumps"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:k-mer-counting
      source: kmer-find-001-report
      evidence: "FindMostFrequentKmers/FindUniqueKmers are filters over the exact CountKmers count map (§2.2, §5.2); FindClumps maintains its own mutable per-window count dictionary (docs/algorithms/K-mer/K-mer_Search.md)"
      confidence: high
      status: current
---

# K-mer search (special-interest k-mers: most-frequent words, singletons, (L,t) clumps)

The **K-mer** family's *pick-out-the-interesting-k-mers* surface: three `KmerAnalyzer` methods that
return k-mers of **special interest** in a sequence rather than the full count map —
`FindMostFrequentKmers` (the **frequent words**, Rosalind BA1B), `FindUniqueKmers` (the **singletons**,
count == 1), and `FindClumps` (patterns forming an **`(L, t)` clump**, Rosalind BA1E). All three are
built on the exact `L − k + 1` k-mer count; the frequent/unique searches are pure **filters** over the
shared [[k-mer-counting]] `CountKmers` map, while `FindClumps` runs its own **sliding-window** count
update with a deduplicating result set. Validated under test unit
[[kmer-find-001-report|KMER-FIND-001]];
[[test-unit-registry]] tracks the unit and [[algorithm-validation-evidence]] describes the
evidence-artifact pattern.

## Core model

For a sequence `S` of length `L` and window length `k`, over the `L − k + 1` overlapping windows with
per-k-mer count `Count(w)`:

| Operation | Method | Definition |
|-----------|--------|------------|
| **Most-frequent** | `FindMostFrequentKmers(string, int)` → `IEnumerable<string>` | `{ w : Count(w) = max_j Count(j) }` — every k-mer **tied** at the maximum observed count (Rosalind BA1B; ties all returned) |
| **Unique / singleton** | `FindUniqueKmers(string, int)` → `IEnumerable<string>` | `{ w : Count(w) = 1 }` — the k-mers occurring exactly once |
| **`(L, t)` clump** | `FindClumps(string, int, int, int)` → `IEnumerable<string>` | a pattern forms a clump if **some** window of length `L` contains **≥ `t`** of its occurrences; returns the deduplicated set of such patterns (Rosalind BA1E) |

A **most-frequent** k-mer maximises `Count`, and multiple k-mers may **tie** for the maximum. A
**unique** k-mer occurs exactly once. A pattern forms an **`(L, t)` clump** if some length-`L` window
contains at least `t` occurrences of it — biologically a signature of motif-rich regions such as
**origins of replication**. Sources: Rosalind BA1B, Rosalind BA1E, Wikipedia (K-mer).

## The clump algorithm (sliding-window count update)

`FindClumps(sequence, k, windowSize L, minOccurrences t)` initialises a mutable `Dictionary<string,int>`
count map over the **first** length-`L` window, records any k-mer whose count reaches `≥ t`, then
**slides** the window one position at a time, updating counts incrementally (drop the k-mer leaving the
left edge, add the one entering the right), and adds newly qualifying patterns to a `HashSet<string>`.
The set is what makes each clump-forming k-mer appear **at most once** in the output (INV-03), regardless
of how many windows support it.

## Properties and invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | every k-mer from `FindUniqueKmers` has `Count == 1` | filters the count dictionary by `value == 1` |
| INV-02 | every k-mer from `FindMostFrequentKmers` has the maximum observed count | filters by `counts.Values.Max()` |
| INV-03 | `FindClumps` returns each qualifying k-mer at most once | results stored in a `HashSet<string>` |
| INV-04 | all-equal-frequency input ⇒ `FindMostFrequentKmers` returns *all* observed k-mers | every k-mer is tied at the maximum |

## Contract and edge cases — the two validation regimes

The three methods split into **two validation regimes**, and the difference is the sharp edge of this
unit:

- **`FindMostFrequentKmers` / `FindUniqueKmers`** reuse [[k-mer-counting]]'s `CountKmers`, so they
  **inherit its validation**: `k ≤ 0` on a non-empty sequence **throws** `ArgumentOutOfRangeException`;
  null/empty sequence or `k > L` → empty result.
- **`FindClumps`** instead treats every invalid parameter as an **empty-result** condition — it **never
  throws**. It returns empty for null/empty sequence, `k ≤ 0`, `windowSize < k`,
  `windowSize > sequence.Length`, or `minOccurrences ≤ 0`.

All three **upper-case** the sequence internally, so matching is **case-insensitive** (shared with the
rest of the K-mer family).

| Case | Behaviour |
|------|-----------|
| Empty / null sequence | empty result (all three) |
| `k ≤ 0` | **throws** for most-frequent / unique; **empty** for `FindClumps` |
| `k > L` | empty result |
| `windowSize > L` | `FindClumps` → empty (no full window) |
| `windowSize < k` | `FindClumps` → empty (a window can't hold a full k-mer) |
| `minOccurrences ≤ 0` | `FindClumps` → empty |
| all k-mers equally frequent | `FindMostFrequentKmers` → all observed k-mers (INV-04) |

## Complexity

| Operation | Time | Space |
|-----------|------|-------|
| `FindMostFrequentKmers` | `O(n)` | `O(u)` |
| `FindUniqueKmers` | `O(n)` | `O(u)` |
| `FindClumps` | `O(n × (L − k + 1))` worst case (near-linear with efficient structures) | `O(u)` |

(`n = L − k + 1`; `u` = distinct k-mers held in the dictionaries/set.) Memory tracks the number of
distinct k-mers maintained — the same dictionary-pressure caveat as [[k-mer-counting]].

## Worked oracles

- `FindMostFrequentKmers("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4)` → `{CATG, GCAT}`, each occurring **3**
  times (Rosalind BA1B frequent-words example).
- `FindClumps("gatcagcataagggtcccTGCAATGCATGACAAGCCTGCAgttgttttac", k=4, L=25, t=3)` includes `TGCA` as a
  `(25, 3)` clump.

## Deviations and assumptions

- **No source contradictions; deviations = None.** The Rosalind BA1B most-frequent-words definition, the
  singleton definition, and the Rosalind BA1E `(L, t)` clump model are implemented verbatim.
- **Intentionally simplified.** `FindClumps` returns **only the set of qualifying k-mers**, not the
  windows in which each qualified — callers learn *which* patterns clump but not *where* each supporting
  window is. Reporting the full set of clump-supporting windows / multiplicity traces is **not
  implemented**; downstream custom analysis is the route if those locations are needed.

## Relation to other units

- **Filters over the shared count.** `FindMostFrequentKmers` and `FindUniqueKmers` filter the exact
  `KmerAnalyzer.CountKmers` multiset owned by [[k-mer-counting]]; `FindClumps` maintains its own
  per-window count map but over the same sliding-window count model.
- **Overlap with [[unique-and-mincount-kmers]] (deliberate).** `FindUniqueKmers` is the *same* singleton
  method synthesized there as part of KMER-UNIQUE-001 (with `FindKmersWithMinCount`); that page owns the
  total/distinct/**unique** terminology and the "unique = singletons ≠ `KmerStatistics.UniqueKmers`
  distinct count" gotcha. This unit adds the **most-frequent** (the `Count ≥ t` head's limiting case:
  the maxima) and the **`(L, t)` clump** operations, which are unrepresented elsewhere.
- **Distinct from [[k-mer-positions]].** That unit locates *where* one specified k-mer occurs (an
  ordered position index for a query pattern); this unit selects *which* k-mers are interesting by
  frequency/clumping and returns the k-mer strings, not their offsets.
- **Distinct from the distribution/summary siblings** [[k-mer-statistics]] (whose `MaxCount` is the
  scalar this unit's most-frequent set attains) and [[k-mer-frequency-analysis]] — those summarise the
  whole profile, whereas this unit extracts named special-interest k-mers.
- **Applications** (from the spec): unique k-mers for marker discovery / genomic fingerprinting; clump
  finding for motif-rich regions such as origins of replication.
