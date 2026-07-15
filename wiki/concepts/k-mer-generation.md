---
type: concept
title: "K-mer generation (n^k universe enumeration)"
tags: [analysis, algorithm]
mcp_tools:
  - generate_all_kmers
sources:
  - docs/Evidence/KMER-GENERATE-001-Evidence.md
  - docs/algorithms/K-mer/K-mer_Generation.md
source_commit: 07cd59444e1d3f85403bf468c40c6de97216c385
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: kmer-generate-001-evidence
      evidence: "Test Unit ID: KMER-GENERATE-001; Algorithm: K-mer Generation (docs/algorithms/K-mer/K-mer_Generation.md)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:both-strand-kmer-counting
      source: kmer-generate-001-evidence
      evidence: "Sibling KmerAnalyzer K-mer-family operation; the default alphabet ACGT and k<=0 -> ArgumentOutOfRangeException contract match the sibling counting methods, but generation enumerates the full n^k universe rather than counting observed substrings"
      confidence: high
      status: current
---

# K-mer generation (n^k universe enumeration)

The **K-mer** family's *space-enumeration* operation: `KmerAnalyzer.GenerateAllKmers(int k, string
alphabet = "ACGT")` returns **every possible** k-mer of length `k` over the alphabet `Σ` — the
complete k-mer universe, `Σ^k`, of cardinality **`n^k`** (`n = |Σ|`; **4^k** for DNA). It is an
exact, deterministic, combinatorial enumeration that depends only on `k` and the alphabet, **not on
any sequence**. Validated under test unit **KMER-GENERATE-001** (record
[[kmer-generate-001-evidence]]); [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the evidence-artifact pattern.

## Generation vs counting — the core distinction

Generation is genuinely **distinct** from the counting operations of the family
([[both-strand-kmer-counting]], [[asynchronous-kmer-counting]], and the frequency-vector basis of
[[k-mer-euclidean-distance]] / [[kmer-jaccard-similarity]]):

| Aspect | K-mer **generation** (this) | K-mer **counting** (siblings) |
|--------|------------------------------|-------------------------------|
| Input | `k` + alphabet only | a sequence `S` (+ `k`) |
| Produces | **all** `n^k` possible words `Σ^k` | only the k-mers **observed** in `S` |
| Output size | exactly `n^k` (fixed by `k`, `Σ`) | ≤ `L − k + 1` distinct, sequence-dependent |
| Values | each word once (a **set**) | occurrence **counts** (a multiset/dictionary) |
| Use | frequency-array indices, background models, motif-space iteration | composition profiling, spectra, similarity |

A generated universe supplies the **address space** (e.g. the `n^k` slots of a frequency array) that
counting then fills; the two are complementary, not interchangeable.

## Core model

The set of all possible k-mers over `Σ` (`|Σ| = n`) is the **k-fold Cartesian product** `Σ^k` —
each of the `k` positions is chosen independently from the `n` symbols. Its size is `n^k`
("there exist n^k total possible k-mers, where n is number of possible monomers", Wikipedia; "the
possible combinations of k positions are computed as 4^k", BioInfoLogics). Enumeration realises
`itertools.product(Σ, repeat=k)`: extend a prefix one character at a time, iterating `Σ` in order at
each position, leftmost position outermost so the **rightmost position varies fastest** (odometer
ordering) — for a **sorted** alphabet this emits k-mers in **lexicographic** order.

## Properties and invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | output count = `n^k` (`4^k` for default DNA `"ACGT"`) | cardinality of the k-fold Cartesian product `Σ^k` |
| INV-02 | all emitted k-mers are distinct (output is a **set** of size `n^k`) | each k-mer is a unique length-`k` tuple over `Σ` |
| INV-03 | every k-mer has length exactly `k` and uses only alphabet characters | direct from the definition of a k-mer over `Σ` |
| INV-04 | for a **sorted** alphabet, emission order is lexicographic | odometer ordering — sorted input ⇒ sorted product tuples |

## Contract and edge cases

- **Inputs:** `k` (`int`, required, `> 0`); `alphabet` (`string`, default `"ACGT"`, non-empty,
  characters used **verbatim** — case-sensitive, no normalisation or de-duplication).
- **Output:** `IEnumerable<string>` of all `alphabet.Length^k` k-mers, streamed **lazily**
  (`yield return`) so callers can iterate the universe without materialising all `n^k` strings.
  Lexicographic when the alphabet is sorted (default `"ACGT"` is).
- **Edge cases:** `k = 1` → one k-mer per symbol (DNA → `A,C,G,T`); single-letter alphabet, any `k`
  → exactly one homopolymer (`1^k = 1`); `k ≤ 0` → `ArgumentOutOfRangeException` (a k-mer length must
  be positive); null/empty alphabet → `ArgumentException` (no symbols ⇒ no k-mers); **unsorted**
  alphabet → all `n^k` k-mers, but in the alphabet's positional order (INV-04 lapses); a **repeated**
  alphabet character yields repeated k-mers (no de-duplication).

## Worked oracles

- `GenerateAllKmers(1)` → `{A, C, G, T}` (4^1 = 4).
- `GenerateAllKmers(2)` → 16 (4^2) 2-mers, lexicographic: `AA, AC, AG, AT, CA, CC, CG, CT, GA, GC,
  GG, GT, TA, TC, TG, TT` — first `AA`, last `TT`.
- `GenerateAllKmers(3)` → 64 (4^3) distinct k-mers, first `AAA`, second `AAC`, last `TTT`.
- `GenerateAllKmers(2, "ACDEFGHIKLMNPQRSTVWY")` → 400 (20^2) protein 2-mers.
- `GenerateAllKmers(4, "A")` → `{AAAA}` (1^4 = 1).

## Complexity, limitations, and implementation

`O(k · n^k)` time (output is `n^k` strings of length `k`), `O(k)` working space via lazy recursion
(`GenerateKmersRecursive` private prefix-extension). Output size grows as `n^k` and is impractical to
**materialise** for large `k` (DNA `k ≥ 14` exceeds 2.6×10⁸ k-mers) — stream the lazy enumerable
instead. A **suffix tree does not apply**: generation *produces* the universe, performing no substring
search against a text.

## Deviations and assumptions

- **No source contradictions; deviations = None.** The `n^k` cardinality (Wikipedia + BioInfoLogics)
  and the Cartesian-product / odometer enumeration (`itertools.product`) are implemented verbatim
  (INV-01..04); nothing is intentionally simplified or left unimplemented.
- **ASSUMPTION (documented property):** the default `alphabet` is `"ACGT"` upper-case, matching the
  sibling `KmerAnalyzer` methods; the alphabet is used verbatim, so the lexicographic guarantee holds
  only when it is itself sorted (the default is).

## Relation to other units

- Same `KmerAnalyzer` K-mer family as [[both-strand-kmer-counting]], [[asynchronous-kmer-counting]],
  and [[k-mer-euclidean-distance]]; it shares their `k ≤ 0` → `ArgumentOutOfRangeException` and
  `"ACGT"` default-alphabet **API conventions**, but is built on an independent recursive
  Cartesian-product enumerator (`GenerateAllKmers`/`GenerateKmersRecursive`), **not** the sync
  `CountKmers` sliding-window primitive the counting siblings share.
- The generated `n^k` universe is the natural **index space** for the frequency vectors that
  [[k-mer-euclidean-distance]] and [[kmer-jaccard-similarity]] compare, and for the k-mer tables used
  in [[de-bruijn-graph-assembly]] and [[kmer-spectrum-error-correction]].
