---
type: concept
title: "Exact pattern search (suffix-tree Contains / Count / FindAll primitives)"
tags: [pattern-matching, algorithm, sequence-comparison]
sources:
  - docs/algorithms/Pattern_Matching/Exact_Pattern_Search.md
source_commit: c37254f8aec5dc7d726bdacffa11dd40809372b6
created: 2026-07-15
updated: 2026-07-15
---

# Exact pattern search (suffix-tree Contains / Count / FindAll primitives)

**Exact pattern search** answers *"where does pattern `P` occur in text `T`?"* — the
classical exact-matching problem, *find every start position `i` where
`T[i..i+m-1] = P`*. This is the **engine layer** (test unit **PAT-EXACT-001**): the three
core `SuffixTree` search primitives plus the two DNA wrappers built on them. It sits
*beneath* the caller-facing motif finders — the multi-pattern
[[known-motif-search]] (`GenomicAnalyzer.FindMotif`) and the single-string repeat
family [[longest-repeated-substring]] both drive this same suffix-tree exact-match
traversal. The artifact/validation pattern for the unit is described by
[[algorithm-validation-evidence]].

## Suffix-tree characterisation

In a suffix tree of `T` (Ukkonen-style `O(n)` construction), every suffix of `T` is a
root-to-leaf path, so `P` occurs at `i` **iff** `P` is a prefix of the suffix beginning
at `i`. Exact search therefore reduces to a single traversal:

1. Start at the root; match `P` character-by-character against edge labels.
2. On the **first mismatch** before `P` is exhausted → `P` is absent.
3. When `P` is exhausted, **every leaf below the matched locus** is one occurrence.

Because each leaf is a distinct suffix = one start position, counting occurrences is
just counting the leaves under the matched node — which the tree precomputes.

## The three core primitives (`SuffixTree`)

| Operation | Returns | Time | Space | How |
|-----------|---------|------|-------|-----|
| `Contains(P)` | `bool` — ≥ 1 occurrence | `O(m)` | `O(1)` aux | matched/not-matched specialization of the traversal (INV-03) |
| `CountOccurrences(P)` | `int` occurrence count | `O(m)` | `O(1)` aux | reads the matched node's precomputed `LeafCount` (no leaf walk) |
| `FindAllOccurrences(P)` | occurrence positions | `O(m + z)` | `O(z)` | collects the `z` leaves under the matched node |

`m` = pattern length, `z` = number of matches. The build itself is `O(n)` time and
space (amortized once per indexed text; the tree is cached per `DnaSequence`).

**Invariants (from the spec):**
- **INV-01** — every reported `i` satisfies `T[i..i+m-1] = P` (leaves collected only after a full edge-by-edge match).
- **INV-02** — `CountOccurrences(P) = FindAllOccurrences(P).Count` on the core API.
- **INV-03** — `Contains(P) ⇔ at least one exact occurrence`.

## The empty-pattern split (the sharp edge)

The **core tree** and the **DNA wrappers** disagree on the degenerate empty pattern — a
deliberate contract difference:

| Input | Core `SuffixTree` | Genomics wrappers |
|-------|-------------------|-------------------|
| `null` pattern | **throws** `ArgumentNullException` | wrappers return **empty** for null/empty motif |
| empty `""` pattern | `Contains` → **true**; `CountOccurrences` → **text length**; `FindAllOccurrences` → **all** positions `[0..n-1]` | returns an **empty** result (wrapper-level guard) |
| pattern longer than text | no matches | no matches |

So the empty pattern is "occurs everywhere" at the primitive level but "no meaningful
occurrence" at the biology level.

## Two DNA wrappers — same engine, different ordering

Both uppercase the motif first (case-insensitive DNA matching) and return empty for a
null/empty motif, but they differ in output ordering — a documented, non-correctness
divergence:

- **`MotifFinder.FindExactMotif(DnaSequence, motif)`** — uppercases, then **sorts** the
  positions before yielding them.
- **`GenomicAnalyzer.FindMotif(DnaSequence, motif)`** — uppercases, then returns the
  underlying suffix-tree position list **directly** (order as produced by leaf
  collection). This is the entry point synthesized by [[known-motif-search]] for the
  multi-motif set case.

## Implementation notes

The repository `SuffixTree` optimizes the hot path: **hybrid SIMD** comparisons for edge
fragments ≥ 8 chars (scalar for shorter), **precomputed `LeafCount`** so
`CountOccurrences` is `O(m)` with no leaf enumeration, and **thread-static buffers** to
cut allocations during traversal and leaf collection. Entry points:
`SuffixTree.Search.cs`, `MotifFinder.cs`, `GenomicAnalyzer.cs`.

## Scope — exact only

This surface is **exact-equality only**. Approximate, mismatch-tolerant, and ambiguity
matching are separate units the spec explicitly redirects to:

- substitution-tolerant (Hamming) → [[approximate-pattern-matching-mismatches]];
- indel-tolerant (Levenshtein) → [[edit-distance]];
- IUPAC ambiguity codes → [[iupac-degenerate-consensus]] (and IUPAC degenerate matching).

The same suffix-tree exact-match engine also backs the two-string
[[longest-common-substring]] and the [[dot-plot-word-match]] word engine — all members
of the deepest-node / exact-occurrence suffix-tree family.

## Worked oracles

| Text | Pattern | Occurrences (0-based) |
|------|---------|-----------------------|
| `banana` | `ana` | `[1, 3]` (overlapping — both reported) |
| `mississippi` | `issi` | `[1, 4]` |
| `GATATATGCATATACTT` | `ATAT` | `[1, 3, 9]` |

Overlapping matches are all reported (the `ana`/`banana` and `ATAT`/`GATATAT...` cases),
consistent with the all-occurrences correctness rule shared with [[known-motif-search]]
and [[k-mer-positions]].

## Sources and deviations

Gusfield (1997) *Algorithms on Strings, Trees and Sequences*, Ukkonen (1995) on-line
suffix-tree construction, Wikipedia "Suffix tree", and Rosalind SUBS (find-a-motif). No
algorithm deviations; the only design choices are the empty-pattern core/wrapper split
and the `MotifFinder`-sorts-vs-`GenomicAnalyzer`-unsorted ordering, neither of which any
source defines as a correctness parameter.
