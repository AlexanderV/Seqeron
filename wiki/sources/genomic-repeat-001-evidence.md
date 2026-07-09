---
type: source
title: "Evidence: GENOMIC-REPEAT-001 (Longest repeated substring + all-repeats enumeration)"
tags: [validation, sequence-comparison]
doc_path: docs/Evidence/GENOMIC-REPEAT-001-Evidence.md
sources:
  - docs/Evidence/GENOMIC-REPEAT-001-Evidence.md
source_commit: c52b50ebde38808027f1f8c3dadf32592547a738
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: GENOMIC-REPEAT-001

The validation-evidence artifact for test unit **GENOMIC-REPEAT-001** — **Repeat
Detection: Longest Repeated Substring (LRS) + all-repeats enumeration** via a suffix
tree (`GenomicAnalyzer.FindLongestRepeat` + `GenomicAnalyzer.FindRepeats`, over
`SuffixTree.LongestRepeatedSubstring` / `FindAllOccurrences` / `GetAllSuffixes`). It is
one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence
artifact]] pattern; the algorithm, its parameters, invariants, worked oracles, and
corner cases are synthesized in [[longest-repeated-substring]]. See
[[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **CMU 15-451/651 Lecture #10 (Suffix Trees and Arrays)** (rank 1) — the verbatim
    LRS definition and algorithm (§2.1): *"Find the longest string `r` such that `r`
    occurs at least twice in `T`: Find the deepest node that has ≥ 2 leaves under it"*;
    the internal-node ⇔ repeated-substring mapping (each leaf = one suffix; a substring
    occurring ≥ 2× is an internal node with ≥ 2 leaves below; occurring once = a leaf);
    O(t)-space linear suffix-tree solution.
  - **Wikipedia — "Longest repeated substring problem"** (rank 4) — problem statement
    (*"longest substring… that occurs at least twice"*), the deepest-internal-node
    (> 1 child) selection rule with a terminal `$`, the Θ(n) linear-time/space bound,
    and the worked example `ATCGATCGA$` → `ATCGA`.
  - **GeeksforGeeks — "Suffix Tree Application 3"** (rank 3) — the verbatim selection
    rule (LRS ends at the internal node farthest from the root, path-label length =
    depth), Ukkonen O(N) build + O(N) deepest-node scan, and worked examples
    `GEEKSFORGEEKS`→`GEEKS`, `AAAAAAAAAA`→`AAAAAAAAA`, `ABCDEFG`→none,
    `ABABABA`→`ABABA`, `banana`→`ana`.
  - **JHU (Langmead) — Suffix Trees** (rank 1) — the Gusfield 5.4 grounding for the
    suffix-tree repeat-application family; internal-node ⇒ shared substring; the
    maximal-repeat note (≤ n maximal repeats, REPuter/Kurtz).
- **Algorithm behaviour (from the artifact):** `FindLongestRepeat(DnaSequence)` returns
  a single `RepeatInfo` (the deepest internal node) or `RepeatInfo.None`;
  `FindRepeats(DnaSequence, minLength)` enumerates *every* distinct substring occurring
  ≥ 2× with length ≥ `minLength`. `RepeatInfo` carries the substring, ascending 0-based
  `Positions`, `Length` (= substring length), `Count` (≥ 2 for any non-`None`), and
  `IsEmpty`. Occurrences **may overlap**.
- **Datasets (documented oracles):**
  - Wikipedia LRS: `ATCGATCGA` → `ATCGA`, length 5, count 2, positions {0, 4}.
  - GeeksforGeeks (DNA analogues — the structural property is alphabet-independent, so
    non-ACGT strings are remapped to A/C/G/T with identical repeat structure since
    `FindLongestRepeat` takes a `DnaSequence`): `AAAAAAAAAA`→`AAAAAAAAA` (overlapping,
    positions {0,1}); `ACGT`→none; `ATATATA` (analogue of `ABABABA`)→`ATATA` (positions
    {0,2}); empty string → `None`.
  - **All-repeats enumeration**, `FindRepeats("ACGTACGTTTTTACGT", 3)` → exactly 8
    substrings: `ACG`@{0,4,12}, `ACGT`@{0,4,12}, `CGT`@{1,5,13}, `TAC`@{3,11},
    `TACG`@{3,11}, `TACGT`@{3,11}, `TTT`@{7,8,9}, `TTTT`@{7,8}. Ground truth is an
    independent **brute-force** enumeration of every distinct substring occurring ≥ 2×
    (cross-checked against the sorted-suffix-LCP prefix set). *(2026-06-15 correction:
    an earlier 5-entry row listed only the maximal-length LCP per adjacent suffix pair
    and omitted the shorter repeated prefixes `ACG`, `TAC`, `TACG` — a code echo of a
    defective implementation; see FINDINGS_REGISTER.)*

## Corner cases and assumptions

Documented corner cases: **no repeat** (every substring unique → no qualifying internal
node → `None`); **overlapping repeats allowed** (the ≥ 2-leaves definition does not
require disjoint occurrences); **ties** — the sources require only *a* longest repeated
substring (Wikipedia), so any equal-length winner is correct.

**Deviations: none for the algorithm.** Two documented **assumptions**: (1) **tie-break**
among equal-length longest repeats — spec requires only "a" longest, repo returns
whichever the deepest-node bookkeeping records (not exercised by any MUST case, all cited
inputs have a unique longest repeat); (2) **occurrence ordering** — `Positions` is sorted
ascending (`OrderBy`), an output-shape convention; the position *set* is definition-fixed.
The matching algorithm doc (`docs/algorithms/Repeat_Analysis/Repeat_Detection.md`) also
records `minLength ≤ 0` clamped to `max(1, minLength)` in `FindRepeats` (a fixed
deviation: without it the empty string was emitted as a zero-length "repeat"). No source
contradictions — CMU, Wikipedia, GeeksforGeeks, and JHU/Gusfield agree on the
deepest-internal-node characterisation and overlap allowance.
