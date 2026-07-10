---
type: source
title: "PAT-APPROX-003 Evidence — approximate pattern matching / frequent words with mismatches"
tags: [validation, testing, motif]
doc_path: docs/Evidence/PAT-APPROX-003-Evidence.md
sources:
  - docs/Evidence/PAT-APPROX-003-Evidence.md
source_commit: dd94a7819ba9764ca0165e710710b83844931da9
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# PAT-APPROX-003 Evidence — approximate pattern matching, frequent words with mismatches

Per-unit [[algorithm-validation-evidence|evidence artifact]] for test unit
**PAT-APPROX-003** ("Best Match and Frequency Analysis"): the Hamming-distance
**approximate pattern-matching** family. Synthesized into the concept
[[approximate-pattern-matching-mismatches]]; see there for the full model. This page is
the concise source record. Full artifact: `docs/Evidence/PAT-APPROX-003-Evidence.md`.

## What the source establishes

- **Definition (BA1H).** A `k`-mer `Pattern` appears in `Text` with at most `d`
  mismatches at any window `Pattern'` with `HammingDistance(Pattern, Pattern') <= d`.
- **`Count_d`** = number of such windows; **frequent words with mismatches (BA1I)** =
  the `k`-mer(s) with the largest `Count_d`, tallied over each window's
  **d-neighborhood** (Hamming ball); the winner need not occur exactly.
- **`Neighbors(Pattern, d)` (BA1N)** = all `k`-mers within Hamming distance `d`, includes
  the identity; size `1 + 3k = 10` for `k=3, d=1`.
- **Practical bound** `k <= 12, d <= 3` (neighbor enumeration is combinatorial).

## Methods covered

`FindApproximateOccurrences`, `CountApproximateOccurrences` (`Count_d`),
`FindFrequentKmersWithMismatches`, `Neighbors`, `FindBestMatch`. Built on the
`HammingDistance` primitive (PAT-APPROX-001 / BA1H).

## Test datasets (oracles)

| Dataset | Params | Expected |
|---------|--------|----------|
| BA1I frequent words | `Text=ACGTTGCATGTCGCATGATGCATGAGAGCT`, k=4, d=1 | {GATG, ATGC, ATGT}, max count 5 |
| BA1H approx. occurrences | `Pattern=ATTCTGGA`, 100-nt Text, d=3 | positions {6,7,26,27,78}, `Count_3=5` |
| Count_1 worked example | `Text=AACAAGCTGATAAACATTTAAAGAG`, `AAAAA`, d=1 | 4 (AACAA, ATAAA, AAACA, AAAGA) |
| BA1N d-neighborhood | `ACG`, d=1 | 10 members incl. ACG itself |

## Deviations and assumptions

- **Deviations: none** — matches the ROSALIND/textbook BA1H/BA1I/BA1N definitions and the
  `go-rosalind` / Rosalind-Solutions reference implementations (O(n·m) scan, neighbor
  tally).
- **One assumption** — `FindBestMatch` returns the **leftmost** minimum-Hamming-distance
  window (exact match short-circuits). A documented API tie-break convention only; it does
  not change the returned minimum distance.

## Corner cases documented

Pattern absent as exact substring (count is over the Hamming ball); all ties returned in
BA1I; `d=0` degenerates to exact matching (`Count_0` exact, `Neighbors(P,0)={P}`);
neighborhood always contains the pattern.

## Sources cited by the artifact

ROSALIND BA1H / BA1I / BA1N (Compeau & Pevzner, *Bioinformatics Algorithms*, ch. 1);
reference implementations charlesreid1/go-rosalind (`rosalind_ba1.go`) and
zonghui0228/Rosalind-Solutions (`rosalind_ba1h.py`).
