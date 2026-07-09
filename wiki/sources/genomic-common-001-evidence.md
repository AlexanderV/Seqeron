---
type: source
title: "Evidence: GENOMIC-COMMON-001 (Longest common substring / common-region detection)"
tags: [validation, sequence-comparison]
doc_path: docs/Evidence/GENOMIC-COMMON-001-Evidence.md
sources:
  - docs/Evidence/GENOMIC-COMMON-001-Evidence.md
source_commit: 60f2b4f40b2211ce94edbea7a5a5928b42b90ce9
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: GENOMIC-COMMON-001

The validation-evidence artifact for test unit **GENOMIC-COMMON-001** — **Longest
Common Substring / Common Region Detection** via a generalized suffix tree
(`FindLongestCommonRegion` + `FindCommonRegions`, over
`SuffixTree.LongestCommonSubstringInfo`). It is one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the
algorithm, its parameters, invariants, worked oracles, and corner cases are
synthesized in [[longest-common-substring]]. See [[test-unit-registry]] for how
units are tracked.

## What this file records

- **Online sources:**
  - **Wikipedia — "Longest common substring"** (rank 4, citing Gusfield 1997) — the
    formal definition (*"Given two strings, S of length m and T of length n, find a
    longest string which is substring of both S and T"*), the **contiguity**
    distinction from the *longest common **subsequence*** ("seeks a contiguous
    substring", not a gapped subsequence), the **Θ(n + m)** generalized-suffix-tree
    complexity claim, and two worked examples — the **tie** case (`BADANAT`/`CANADAS`
    share *two* distinct maximal substrings `ADA` and `ANA`) and the 3-string unique
    case (`ABABC`/`BABCA`/`ABCBA` → only `ABC`).
  - **GeeksforGeeks — "Suffix Tree Application 5 — Longest Common Substring"**
    (rank 3, reference-implementation description) — worked output `xabxac` /
    `abcabxabcd` → **`abxa`, length 4**; the **generalized-suffix-tree (GST)
    mechanism** (the LCS is the path label root → *the deepest internal node whose
    subtree contains leaves from **both** strings*); and the **O(M+N)** linear
    time-and-space bound (GST build + a DFS).
- **Algorithm behaviour (from the artifact):** `FindLongestCommonRegion` returns the
  maximal **contiguous** shared substring with its length and **0-based** start
  positions in both sequences; no common character → `CommonRegion.None` (empty,
  length 0, positions −1); identical sequences → the whole sequence at positions 0/0.
  `FindCommonRegions(minLength)` enumerates *all* distinct contiguous shared
  substrings of length ≥ `minLength` with their positions.
- **Datasets (documented oracles):**
  - *Tie / determinism* (Wikipedia tie property, DNA analogue): `CACAGAG` vs
    `TACATAGAT` share two distinct maximal length-3 substrings `ACA` and `AGA`; the
    documented **first-found-in-`other`** tie-break selects `ACA` (it ends earlier in
    `TACATAGAT`). The GeeksforGeeks literal strings are non-DNA, so only their
    *length + contiguity* property transfers; DNA cases were cross-checked by an
    independent **O(n³) brute-force** `s[i:j] in t` enumeration (not the repo impl).

## Deviations and assumptions

**Deviations: none.** One documented **assumption** — the **tie-break rule**: when
several distinct substrings share the maximal length, no authoritative source
mandates *which* to return (Wikipedia reports all; GeeksforGeeks returns one). The
repository's `SuffixTree.LongestCommonSubstringInfo` documents and implements
"**the first one found in `other` is returned**"
(`src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Algorithms.cs`, XML doc). This is a
deterministic, documented choice that returns *a* correct LCS — it changes only the
representative, never which lengths are maximal.

Recommended coverage (MUST): maximal contiguous region with length + 0-based
positions; deterministic tie (`CACAGAG`/`TACATAGAT` → `ACA`); no common substring →
`None`; identical → whole sequence at 0/0; `FindCommonRegions(minLength)` enumerates
all distinct shared substrings ≥ minLength. SHOULD: empty input → `None`. COULD:
invariant that the returned substring actually occurs at the reported positions in
both sequences. No contradictions among sources — Wikipedia and GeeksforGeeks agree
on the contiguity definition and the GST mechanism; they differ only on *reporting*
ties (all vs one), which the repo resolves deterministically.
