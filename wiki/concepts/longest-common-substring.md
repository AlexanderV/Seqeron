---
type: concept
title: "Longest common substring / common-region detection (generalized suffix tree)"
tags: [sequence-comparison, algorithm]
sources:
  - docs/Evidence/GENOMIC-COMMON-001-Evidence.md
  - docs/algorithms/Sequence_Comparison/Common_Region_Detection.md
source_commit: 60f2b4f40b2211ce94edbea7a5a5928b42b90ce9
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: genomic-common-001-evidence
      evidence: "Test Unit ID: GENOMIC-COMMON-001 ... Algorithm: Longest Common Substring / Common Region Detection (generalized suffix tree)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:dot-plot-word-match
      source: genomic-common-001-evidence
      evidence: "Both are exact-match sequence-comparison ops built on the same generalized suffix tree — the LCS is the deepest two-string internal node, the dot plot's exact-word engine finds occurrences via the same suffix tree (docs/algorithms/Comparative_Genomics/Dot_Plot_Generation.md)"
      confidence: medium
      status: current
---

# Longest common substring / common-region detection (generalized suffix tree)

The **longest common substring (LCS)** of two strings S (length m) and T (length n)
is *"a longest string which is substring of both S and T"* (Wikipedia). The defining
property is **contiguity**: unlike the longest common *subsequence*, which allows
insertions/deletions within the shared text, the LCS problem seeks a **contiguous**
run shared by both inputs. Seqeron exposes this as **common-region detection** —
`FindLongestCommonRegion` (the single maximal shared region) and
`FindCommonRegions(minLength)` (every distinct shared region ≥ `minLength`), both
built on `SuffixTree.LongestCommonSubstringInfo`. Validated under test unit
**GENOMIC-COMMON-001**; the validation record is [[genomic-common-001-evidence]],
[[test-unit-registry]] tracks the unit, and [[algorithm-validation-evidence]]
describes the artifact pattern.

## Substring, not subsequence

The one distinction a correct implementation must honour: the returned region is a
**contiguous** substring of *both* inputs, **not** a gapped subsequence. `FindCommonRegions`
returns exact shared runs — nothing is stitched across a gap. (The gapped/edit view
of similarity lives elsewhere: edit distance, alignment; this unit is the
exact-contiguous-match view.)

## The generalized-suffix-tree algorithm

The classic solution (Gusfield 1997, via Wikipedia) builds a **generalized suffix
tree (GST)** over both strings and finds the LCS in **Θ(n + m)** time:

1. Build a GST containing every suffix of S and of T, each leaf tagged with which
   string it came from.
2. The LCS is the **path label from the root to the deepest internal node whose
   subtree contains leaves from *both* strings** (GeeksforGeeks). The deepest such
   node's path label is the longest run present in both inputs.
3. Finding it is a single DFS over the tree.

Both stages are linear: GST construction is O(M + N) and the DFS is O(M + N), so the
whole operation is **linear in time and space** (GeeksforGeeks). This is the
exact-match engine Seqeron also reuses for the [[dot-plot-word-match|dot plot]]'s
word occurrences. The **single-string** counterpart — the longest substring recurring
within *one* sequence — is [[longest-repeated-substring]] (deepest internal node with
≥ 2 leaves, rather than leaves from both strings).

## API contract and invariants

| Operation | Returns |
|-----------|---------|
| `FindLongestCommonRegion(a, b)` | the maximal contiguous shared region: substring, length, and **0-based** start positions in *both* sequences |
| `FindCommonRegions(a, b, minLength)` | *all* distinct contiguous shared substrings of length ≥ `minLength`, each with its positions |

- **No common substring** → `CommonRegion.None` (empty string, length 0, positions
  **−1**). The empty string is a substring of both, so length 0 is the correct floor
  (Wikipedia definition).
- **Identical sequences** → the whole sequence is the LCS at positions **0 / 0** (a
  string is a substring of itself).
- **Empty input** → `None` / empty enumeration.
- **Position invariant**: the returned substring actually occurs at the reported
  positions in both sequences (a correctness cross-check).

## Ties — the deterministic tie-break

Several distinct substrings can share the maximal length. Wikipedia's example
`BADANAT` / `CANADAS` shares **two** maximal substrings, `ADA` and `ANA`; no
authoritative source mandates *which* one to return (Wikipedia reports **all**,
GeeksforGeeks returns **one**). Seqeron makes a documented deterministic choice:
**"the first one found in `other` is returned"** (`SuffixTree.LongestCommonSubstringInfo`
XML doc). This returns *a* correct LCS — it fixes only the representative, never which
lengths are maximal.

- **Oracle** (DNA analogue of the Wikipedia tie): `CACAGAG` vs `TACATAGAT` share two
  distinct maximal length-3 substrings `ACA` and `AGA`; the first-found-in-`other`
  rule selects **`ACA`** (it ends earlier in `TACATAGAT`).

## Reference sources

**Wikipedia — "Longest common substring"** (the formal definition, the
substring-vs-subsequence contiguity distinction, the Θ(n+m) GST claim, and the tie /
3-string worked examples) and **GeeksforGeeks — "Suffix Tree Application 5"** (the
`xabxac`/`abcabxabcd` → `abxa` length-4 example, the "deepest two-string internal
node" GST mechanism, and the O(M+N) bound), with the primary attribution to
**Gusfield 1997** (*Algorithms on Strings, Trees and Sequences*) for the linear GST
result. **No deviations**; the only assumption is the deterministic first-in-`other`
tie-break, which the sources leave unspecified.
