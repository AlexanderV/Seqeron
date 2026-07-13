---
type: concept
title: "K-mer counting (the canonical synchronous sliding-window count primitive)"
tags: [analysis, algorithm]
mcp_tools:
  - count_kmers
sources:
  - docs/algorithms/K-mer/K-mer_Counting.md
source_commit: a0600dbbba62a14760fc3c8398e37ab9405ff37c
created: 2026-07-13
updated: 2026-07-13
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: k-mer-counting
      evidence: "Test Unit ID: KMER-COUNT-001; Algorithm: K-mer Counting (docs/algorithms/K-mer/K-mer_Counting.md), KmerAnalyzer.CountKmers"
      confidence: high
      status: current
---

# K-mer counting (the canonical synchronous sliding-window count primitive)

`KmerAnalyzer.CountKmers` is the **foundational counting primitive** of the K-mer family: it slides a
window of length `k` across a sequence and returns a `Dictionary<string,int>` mapping each observed
length-`k` substring to the number of start positions that produced it. Every other K-mer counting
operation in the repository is a wrapper, variant surface, filter, or summary **over this method** —
the cancelable/progress [[asynchronous-kmer-counting]] delegates to it verbatim, the strand-aware
[[both-strand-kmer-counting]] calls it once on `S` and once on `RC(S)`, [[k-mer-statistics]] aggregates
its multiset, and [[unique-and-mincount-kmers]] filters it. This page owns the shared **`L − k + 1`
count definition** those pages defer to. Validated under test unit **KMER-COUNT-001**;
[[test-unit-registry]] tracks the unit and [[algorithm-validation-evidence]] describes the
evidence-artifact pattern.

## Core model

For a sequence `S` of length `L` and window length `k` with `1 ≤ k ≤ L`, the k-mers are the multiset
of the `L − k + 1` overlapping length-`k` windows `{ S[i..i+k) : 0 ≤ i ≤ L − k }`. The count of a
k-mer `w` is the number of start positions yielding `w`:

> **Count(w) = Σ_{i=0}^{L−k} 1( S[i..i+k) = w )**

stored in a dictionary keyed by the **uppercased** k-mer string. For a DNA alphabet of size 4 there are
`4^k` possible distinct k-mers, so the observed distinct count is at most `min(4^k, L − k + 1)` — but
that combinatorial bound applies **only** when the input is restricted to `A/C/G/T` (see the
alphabet note below). Sources: Wikipedia (K-mer: "a sequence of length L will have L − k + 1 k-mers
and there exist nᵏ total possible k-mers"), Rosalind (K-mer Composition), Compeau et al. 2011,
Marçais & Kingsford 2011.

## Entry-point surfaces

`CountKmers` is exposed through several overloads that all compute the identical count; they differ in
input type and execution concern, not arithmetic:

| Surface | Signature (class) | Distinguishing concern |
|---------|-------------------|------------------------|
| Canonical string | `CountKmers(string, int)` (`KmerAnalyzer`) | the reference implementation |
| `DnaSequence` wrapper | `CountKmers(DnaSequence, int)` (`KmerAnalyzer`) | delegates to the string path via `.Sequence` |
| Span-based | `CountKmersSpan(ReadOnlySpan<char>, int)` (`KmerAnalyzer` / `SequenceExtensions`) | allocation-lean span input; **validates `k` before length** (see gotcha) |
| Cancellation-aware | `CountKmers(string, int, CancellationToken, IProgress<double>?)` | cooperative cancellation + progress — see [[asynchronous-kmer-counting]] |
| Async | `CountKmersAsync(...)` | thread-pool offload — owned by [[asynchronous-kmer-counting]] |
| Both-strand | `CountKmersBothStrands(DnaSequence, int)` | forward + reverse-complement sum — owned by [[both-strand-kmer-counting]] |

## Properties and invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `Σ_w Count(w) = L − k + 1` whenever `1 ≤ k ≤ L` | the sliding window emits exactly one k-mer per valid start position |
| INV-02 | for DNA-alphabet input, distinct k-mers ≤ `min(4^k, L − k + 1)` | at most `4^k` DNA words and at most one emitted k-mer per window |
| INV-03 | counting is case-insensitive | string-based input is upper-cased before keys are emitted |
| INV-04 | empty/null `S` or `k > L` ⇒ empty dictionary | `L − k + 1 ≤ 0` ⇒ no valid windows |

## Contract and edge cases

- **Inputs:** `sequence` (`string`, `DnaSequence`, or `ReadOnlySpan<char>`; null/empty → empty
  dictionary; string/cancellation-aware paths upper-case internally so counting is
  **case-insensitive**), `k` (`int`; `k > L` → empty dictionary). Cancellation-aware overloads also
  take `CancellationToken` and optional `IProgress<double>`.
- **Output:** `Dictionary<string,int>` mapping each observed k-mer to its occurrence count.
- **Alphabet-agnostic:** the raw-string and span surfaces do **not** restrict the alphabet — ambiguous
  or non-ACGT symbols (e.g. IUPAC `N`) are counted literally as k-mer characters after upper-casing.
  The `4^k` bound therefore holds only for DNA-restricted input.
- **Edge cases:** empty/null → empty dictionary; `k > L` → empty dictionary; homopolymer (`AAAA`, k=2)
  → one key with count `L − k + 1`; lowercase/mixed case → same counts as uppercase.

## Gotcha: the `k ≤ 0` validation order differs by surface

The invalid-`k` guard fires at a **different point** depending on the surface, producing a genuinely
observable difference for empty input:

- **String / cancellation-aware overloads** short-circuit on null/empty input **before** validating
  `k`. So an **empty string with `k ≤ 0` returns an empty dictionary** (no throw), while a **non-empty
  string with `k ≤ 0` throws `ArgumentOutOfRangeException`**.
- **The span-based path validates `k` first**, before checking the sequence length — so span input with
  `k ≤ 0` throws `ArgumentOutOfRangeException` regardless of whether the span is empty.

The `DnaSequence` overloads delegate to the string path, so an empty `DnaSequence.Sequence` inherits
the string short-circuit (empty result before `k` is checked).

## Complexity

`O(n · k)` **effective** time (`n = L − k + 1` windows; each window **materializes and hashes a
length-`k` string key**, so the per-window cost is `O(k)`, not `O(1)`), `O(u · k)` space where `u` is
the number of distinct k-mers stored. `CountKmersBothStrands` runs the same scan twice (over `S` and
`RC(S)`). Large genomes or large `k` create substantial dictionary pressure precisely because of the
string-key materialization — the primitive is exact and general-purpose rather than memory-optimized
(no 2-bit packing or minimizer bucketing). A **suffix tree was evaluated and deliberately not used**
across the family: full-spectrum counting is a single linear count-all-windows pass, not a
fixed-pattern occurrence query, so a suffix tree would add `O(n)` construction with no benefit (the
same reuse-policy note the sibling counting pages record).

## Worked oracles

- `CountKmers("ATGG", 3)` → `{ ATG:1, TGG:1 }` (Wikipedia: "ATGG has two 3-mers"; `4 − 3 + 1 = 2`).
- `GTAGAGCTGT` (`L = 10`), total / distinct: **k=2 → 9 / 7** (GT and AG each occur twice), **k=3 →
  8 / 8** (all windows distinct), **k=4 → 7 / 7** (all distinct).
- `ATCGATCAC` (`L = 9`) k=3 → 7 total, 6 distinct (`ATC` occurs twice; the other 5 are singletons).
- `AAAA` (`L = 4`) k=2 → `{ AA:3 }` — homopolymer, every window identical, one key with count `L − k + 1`.

## Deviations and assumptions

- **No source contradictions; deviations = None** for the core sliding-window count — the Wikipedia /
  Rosalind `L − k + 1` model is implemented verbatim.
- **Intentionally simplified.** (1) The raw-string and span surfaces do **not** enforce a DNA alphabet,
  so the `4^k` combinatorial bound applies only to DNA-restricted input; (2) both-strand counting
  **sums** forward and reverse-complement counts rather than collapsing each k-mer to a single
  canonical key (owned by [[both-strand-kmer-counting]]).
- **Not implemented.** Canonical reverse-complement collapsing as the default representation (Jellyfish
  `-C` / Mash), and memory-optimized encodings (2-bit packing, minimizers) — downstream
  post-processing is the route if canonical keys are required.

## Relation to other units

- **Foundation of the counting family.** [[asynchronous-kmer-counting]] wraps the cancellation-aware
  overload for non-blocking execution (identical numeric result); [[both-strand-kmer-counting]] calls
  it on `S` and `RC(S)`; [[k-mer-statistics]] reduces its multiset to a `KmerStatistics` summary
  (including the Shannon k-entropy); [[unique-and-mincount-kmers]] filters it by per-k-mer count.
- **Distinct siblings that do *not* build on this count.** [[k-mer-generation]] enumerates the full
  `n^k` **universe** (sequence-independent, produces the address space rather than counting observed
  windows); [[k-mer-positions]] indexes **where** one k-mer occurs (the inverse index) rather than
  tallying how many; [[k-mer-euclidean-distance]] and [[kmer-jaccard-similarity]] compare the
  frequency/presence profiles of **two** sequences.
- **Downstream consumers.** The count profile feeds [[de-bruijn-graph-assembly]] and
  [[kmer-spectrum-error-correction]].
