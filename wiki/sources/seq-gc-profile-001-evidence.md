---
type: source
title: "Evidence: SEQ-GC-PROFILE-001 (sliding-window GC-content profile)"
tags: [validation, sequence-statistics, composition, chromosome]
doc_path: docs/Evidence/SEQ-GC-PROFILE-001-Evidence.md
sources:
  - docs/Evidence/SEQ-GC-PROFILE-001-Evidence.md
source_commit: 599fc94985c5e39969feee53560e6db69d7bb21f
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SEQ-GC-PROFILE-001

The validation-evidence artifact for test unit **SEQ-GC-PROFILE-001** — a **sliding-window
GC-content profile**: a fixed-width window slides along the sequence and each fully-contained
window emits its **GC% = `(G + C)/(A + T + G + C) × 100`** (the ×100 percentage form). This is
the **standalone GC%-only scan** that the composite `GcAnalysisResult`
(**SEQ-GC-ANALYSIS-001**, [[seq-gc-analysis-001-evidence]]) wraps as its windowed GC-content
channel — the **same measure and same window geometry**, just the GC% channel emitted on its
own (no GC-skew profile, no compositional variance). The method itself is written up on
[[windowed-gc-profile-and-variance]]; this page records the source trace and worked oracles.
One instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern; [[test-unit-registry]] tracks the unit, and [[seq-entropy-profile-001-evidence]] is the
information-theoretic sibling (the standalone windowed-entropy profile).

## What this file records

- **The windowed GC% value** — per fully-contained window, `GC% = (G+C)/(A+T+G+C) × 100`, a
  value bounded in `[0, 100]` (numerator ≤ denominator). Window count =
  `⌊(n − w)/step⌋ + 1` with offsets `0, step, 2·step, …`; a sequence **shorter than one window**
  (and null/empty) ⇒ **empty profile**.
- **Denominator excludes non-standard symbols (N etc.)** — only A/T/U/G/C count toward the
  length; N and other non-standard symbols are removed. This matches Biopython's default
  `ambiguous="remove"` and the Wikipedia A+T+G+C denominator (window `GGAN` → `2/3 × 100 =
  66.66…%`, N excluded).
- **RNA U is a non-GC base** (equivalent to T); a window with U scores the same as with T.
- **Online sources:**
  - **Wikipedia "GC-content"** (rank 4, citing primaries) — the closed-form
    `GC% = (G+C)/(A+T+G+C) × 100%`; denominator is the four standard bases (ambiguous symbols
    not part of the standard-base denominator); also records the `(A+T)/(G+C)` AT/GC ratio (not
    used by this unit). Undefined when `A+T+G+C = 0` (division by zero).
  - **Biopython `Bio.SeqUtils.gc_fraction`** (rank 3, reference implementation) — returns a
    `[0,1]` fraction (`×100` ⇒ the Wikipedia percentage); `S` counts as GC; default
    `ambiguous="remove"` drops non-`ACTGSWU` symbols from the denominator, `ambiguous="ignore"`
    keeps all symbols. Doctests (verbatim): `gc_fraction("ACTGN","ignore")` → `0.40` (len 5),
    `gc_fraction("ACTGN","remove")` → `0.50` (len 4 after removing N), `gc_fraction("ACTG")` →
    `0.50`, RNA `"GGAUCUUCGGAUCU"` → `0.50`. `GC123` confirms the ×100 percentage is the
    standard external presentation.
- **Datasets:**
  - Biopython `gc_fraction` doctests scaled ×100: `ACTG` remove → 50.0, `ACTGN` remove → 50.0,
    `ACTGN` ignore → 40.0, RNA `GGAUCUUCGGAUCU` remove → 50.0.
  - Hand-derived windows from the Wikipedia formula (×100): `GGGG` → 100.0, `AAAA` → 0.0,
    `ATGC` → 50.0, `GGGA` → 75.0, `GCAT` → 50.0, `GGAN` (N excluded) → 66.66666666666666.
- **Corner cases / failure modes:** window with **no standard base** (all-N) ⇒ `A+T+G+C = 0`
  division by zero → the repository returns **0** (convention, see below); `windowSize > length`
  / null / empty ⇒ empty profile; case-insensitive (case-folded before counting); every window
  value bounded `[0, 100]`.

## Deviations and assumptions

**Two documented assumptions, only one still open:**

1. **Empty-window convention — window with no standard base ⇒ GC% `0`** (open assumption).
   Wikipedia and Biopython leave GC content **undefined** when `A+T+G+C = 0` (division by zero);
   the repository returns `0`. This mirrors the sibling `GcSkewCalculator` (SEQ-GC-ANALYSIS-001),
   which returns `GcContent = 0` for the no-G/C / zero-division case, so it is **consistent
   within the repository**, but it is not dictated by the external sources. Applies only to the
   degenerate all-N window.
2. **Denominator excludes non-standard symbols (N etc.)** — **resolved by evidence**, not an
   open assumption: Biopython's default `ambiguous="remove"` and the Wikipedia A+T+G+C
   denominator both exclude N. Counting only A/T/U/G/C matches the `remove` convention.

Recommended coverage (from the artifact): MUST — each window's GC% `= (G+C)/(A+T+G+C)×100`
(exact 0/50/75/100), N excluded from the denominator (`GGAN` → 66.6…%), window count
`⌊(n−w)/step⌋+1` and offsets `0, step, 2·step, …`, RNA U counts as non-GC. SHOULD — all-N
window ⇒ 0; `windowSize > length` / null / empty ⇒ empty profile. COULD — case-insensitivity;
every window bounded `[0, 100]`. No source contradictions — Wikipedia GC-content and Biopython
`gc_fraction` agree on the numerator (G+C), the A+T+G+C denominator, the `remove` N-handling,
and the ×100 percentage presentation.
