---
type: concept
title: "Longest repeated substring + repeat enumeration (single-string suffix tree)"
tags: [sequence-comparison, algorithm]
sources:
  - docs/Evidence/GENOMIC-REPEAT-001-Evidence.md
  - docs/algorithms/Repeat_Analysis/Repeat_Detection.md
source_commit: c52b50ebde38808027f1f8c3dadf32592547a738
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: genomic-repeat-001-evidence
      evidence: "Test Unit ID: GENOMIC-REPEAT-001 ... Algorithm: Repeat Detection — Longest Repeated Substring (LRS) and all repeated substrings via suffix tree"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:longest-common-substring
      source: genomic-repeat-001-evidence
      evidence: "Both are exact suffix-tree deepest-internal-node problems on GenomicAnalyzer: the LRS is the deepest internal node with ≥2 leaves in a SINGLE string's suffix tree; the LCS is the deepest internal node with leaves from BOTH strings in a generalized suffix tree (docs/algorithms/Sequence_Comparison/Common_Region_Detection.md)."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:repetitive-element-detection
      source: genomic-repeat-001-evidence
      evidence: "docs/algorithms/Repeat_Analysis/Repeat_Detection.md §2.5 contrasts LRS/repeat-enumeration (any recurring substring via the suffix tree, occurrences anywhere, may overlap) with FindTandemRepeats (consecutive repeating unit via a windowed scan) — sibling repeat sub-problems in the repeats family."
      confidence: medium
      status: current
---

# Longest repeated substring + repeat enumeration (single-string suffix tree)

Finding the substrings that **recur within one DNA sequence** — the exact-match
counterpart to the two-string [[longest-common-substring]]. Seqeron exposes two
operations on `GenomicAnalyzer`, validated together as test unit **GENOMIC-REPEAT-001**:

- **`FindLongestRepeat(DnaSequence)`** — the **Longest Repeated Substring (LRS)**: the
  longest substring occurring **at least twice**.
- **`FindRepeats(DnaSequence, minLength)`** — enumerates **every** distinct substring
  occurring ≥ 2× with length ≥ `minLength`.

The validation record is [[genomic-repeat-001-evidence]]; [[test-unit-registry]] tracks
the unit and [[algorithm-validation-evidence]] describes the artifact pattern. This is a
*distinct* operation from the annotation-family [[repetitive-element-detection]] anchor
(tandem / inverted / RepeatMasker-class assignment) — it finds **any** recurring
substring, dispersed or overlapping, not a consecutive repeat unit.

## The suffix-tree characterisation

In a suffix tree of the text (with a terminal sentinel), **every leaf is a distinct
suffix = one position**, so:

- a substring occurring **k times** is the path label of an **internal node with k
  leaves** below it;
- a substring occurring **once** ends inside a **leaf edge**.

Repeat finding is therefore a traversal problem (the classical Gusfield ch. 5–7
application family).

**LRS** = the path label of the **deepest internal node** (measured by *string depth* =
number of characters from root, which equals the repeat length) that has **≥ 2 leaves**
below it (CMU §2.1: *"Find the deepest node that has ≥ 2 leaves under it"*; Wikipedia:
deepest internal node with more than one child). One O(n) build + one deepest-node query
→ **linear**.

**All repeats ≥ minLength**: every substring occurring ≥ 2× is the longest common prefix
(LCP) of two suffixes. Seqeron sorts the suffixes, takes each adjacent-pair LCP, and — for
completeness — emits **every prefix** (length `max(1, minLength)` … LCP-length) that
occurs ≥ 2×, not just the full LCP. (The single-echo-of-the-full-LCP bug that omitted
shorter repeated prefixes like `ACG`/`TAC` was the FINDINGS_REGISTER correction.)
`FindRepeats` is worst-case **O(n²)** time and space — unsuitable for whole-genome input.

## LRS vs LCS — same engine, one vs two strings

| | LRS (this unit, GENOMIC-REPEAT-001) | LCS ([[longest-common-substring]], GENOMIC-COMMON-001) |
|---|---|---|
| Input | one sequence | two sequences |
| Tree | suffix tree of `T` | **generalized** suffix tree over both |
| Winner | deepest internal node with **≥ 2 leaves** | deepest internal node with leaves from **both** strings |
| API | `FindLongestRepeat` / `FindRepeats` | `FindLongestCommonRegion` / `FindCommonRegions` |

Both are the exact-match deepest-internal-node family Seqeron builds on its repository
`SuffixTree` (also reused by the [[dot-plot-word-match]] word engine).

## API contract and invariants

`RepeatInfo` carries `Sequence`, ascending **0-based** `Positions`, `Length` (=
`Sequence.Length`), `Count` (= `Positions.Count`), and `IsEmpty`.

- **INV**: every returned repeat has `Count ≥ 2` (internal node, never a leaf); `Length`
  = root-to-node string depth; each position is a true 0-based occurrence.
- **No repeat** (e.g. `ACGT`) or **empty sequence** → `FindLongestRepeat` = `RepeatInfo.None`,
  `FindRepeats` = empty enumeration.
- **Overlaps counted**: `AAAAAAAAAA` → `AAAAAAAAA` at {0,1}; `ATATATA` → `ATATA` at {0,2}.
- **`minLength ≤ 0`** clamped to `max(1, minLength)` → only substrings occurring ≥ 2×, no
  zero-length "repeats".
- Restricted to the ACGT alphabet (`DnaSequence`); does **not** classify maximal vs
  supermaximal repeats and does **not** search the reverse complement (dispersed inverted
  repeats) — see [[repetitive-element-detection]] for the inverted/tandem sub-problems.

## Worked oracles

- `FindLongestRepeat("ATCGATCGA")` → `ATCGA`, length 5, count 2, positions {0, 4}
  (Wikipedia `ATCGATCGA$`→`ATCGA`).
- `FindLongestRepeat("AAAAAAAAAA")` → `AAAAAAAAA`, {0,1}; `FindLongestRepeat("ATATATA")`
  → `ATATA`, {0,2}; `FindLongestRepeat("ACGT")` / `""` → `None`.
- `FindRepeats("ACGTACGTTTTTACGT", 3)` → exactly 8 substrings: `ACG`@{0,4,12},
  `ACGT`@{0,4,12}, `CGT`@{1,5,13}, `TAC`@{3,11}, `TACG`@{3,11}, `TACGT`@{3,11},
  `TTT`@{7,8,9}, `TTTT`@{7,8}.

## Sources and deviations

CMU 15-451 Lecture #10 (§2.1 verbatim LRS definition) + Wikipedia "Longest repeated
substring problem" (deepest-internal-node rule, Θ(n), `ATCGATCGA`→`ATCGA`) + GeeksforGeeks
"Suffix Tree Application 3" (selection rule + worked examples) + JHU/Langmead (Gusfield 5.4
grounding, maximal-repeat note). **No algorithm deviations**; two assumptions — the
**tie-break** among equal-length longest repeats (spec requires only *a* longest; not
exercised by any MUST case) and the **ascending `Positions`** output convention.
