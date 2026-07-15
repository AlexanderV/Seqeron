---
type: source
title: "Evidence: SEQ-GC-ANALYSIS-001 (Comprehensive GC analysis — GC content, GC/AT skew, windowed profiles, compositional variance)"
tags: [validation, sequence-statistics, composition, chromosome]
doc_path: docs/Evidence/SEQ-GC-ANALYSIS-001-Evidence.md
sources:
  - docs/Evidence/SEQ-GC-ANALYSIS-001-Evidence.md
source_commit: f6fc5f03fffb7fd2053db36d0ad79995b8affe3e
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SEQ-GC-ANALYSIS-001

The validation-evidence artifact for test unit **SEQ-GC-ANALYSIS-001** — **Comprehensive GC
analysis**, the composite `GcAnalysisResult` that bundles the whole-sequence GC/AT scalars
with a **sliding-window profile** and its **population variance**. It is one instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern;
[[test-unit-registry]] tracks the unit.

**No new concept for the overall scalars — they are already covered.** The overall GC content
is [[base-composition]]; the overall GC skew and AT skew are [[nucleotide-composition-skew]].
What this unit *adds* — the windowed GC-content / GC-skew profiles and the population-variance
summaries of those windows — is synthesized on the concept
[[windowed-gc-profile-and-variance]]. This page records only what the artifact itself pins.

## What this file records

- **Six outputs of `GcAnalysisResult`** (per the worked datasets):
  - **OverallGcContent (%)** = `(G + C) / (A + T + G + C) × 100` — the ×100 *percentage*
    convention (Brock textbook / Madigan & Martinko 2003), not Biopython's `[0,1]`
    `gc_fraction`; the two differ only by the factor 100. `GGGCCAT` → `5/7×100 = 71.428…`.
  - **OverallGcSkew** = `(G − C)/(G + C)` (Lobry 1996). `GGGCCAT` → `1/5 = 0.2`.
  - **OverallAtSkew** = `(A − T)/(A + T)` (Charneski/Lobry). `GGGCCAT` → `0/2 = 0.0`.
  - **windowed GC%/GC-skew lists** — per fully-contained window `(WindowStart, WindowEnd,
    Position)`; a sequence shorter than one window → **empty lists**.
  - **GcSkewVariance / GcContentVariance** — the **population variance** `σ² = Σ(xᵢ−μ)²/N`
    (÷N, *not* the Bessel-corrected ÷(N−1)) of the per-window values. Window set `GG`/`CC`
    (w=2, s=2 over `GGCC`): skews `+1`/`−1` → `GcSkewVariance = ((1−0)²+(−1−0)²)/2 = 1.0`;
    GC% `100`/`100` → `GcContentVariance = 0.0`. Empty window set → variances `0`.
- **Online sources:**
  - **Wikipedia "GC-content"** (rank 4, citing Madigan & Martinko 2003 *Brock Biology of
    Microorganisms* 10th ed.) — the `GC% = (G+C)/(A+T+G+C)×100%` formula; denominator is all
    four bases.
  - **Wikipedia "GC skew"** (rank 4, citing Lobry 1996 + Grigoriev 1998) — `GC skew =
    (G−C)/(G+C)`, spectrum `−1…+1` (`+1 ⇔ C=0`, `−1 ⇔ G=0`); a **sign switch locates the
    replication origin/terminus** (leading vs lagging strand), Grigoriev 1998's cumulative-skew
    diagram.
  - **Biopython `Bio.SeqUtils`** (rank 3) — `GC_skew` "for multiple windows along the
    sequence" (the windowing precedent), **zero-division ⇒ 0** when `G+C=0`, **ambiguous bases
    ignored** (only G and C count); `gc_fraction` returns a `[0,1]` fraction (= the percentage
    ÷100).
  - **Cuemath "Population Variance"** (rank 4) — `σ² = Σ(xᵢ−μ)²/n` (÷N), worked anchor
    `{12,13,12,14,19}` → μ=14, Σ(xᵢ−μ)²=34, variance `34/5 = 6.8`; pins the estimator choice
    independently of the implementation.
- **Corner cases / failure modes:** window with no G/C ⇒ skew `0` (not NaN); non-ACGT symbols
  ignored; pure-G window ⇒ GC skew `+1`, pure-C ⇒ `−1`; no-G/C sequence ⇒ GC% `0`; **sequence
  shorter than the window ⇒ empty windowed lists ⇒ window-derived variances `0`, while the
  overall scalars are still defined over the whole sequence**; null `DnaSequence` ⇒
  `ArgumentNullException`, null/empty string ⇒ zero/empty result.

## Deviations and assumptions

**Two documented assumptions, both labelling/estimator choices — neither correctness-affecting
at the formula level:**

1. **GC content reported as a percentage (×100), not Biopython's `[0,1]` fraction** — both are
   documented conventions (Brock ×100% vs `gc_fraction` `[0,1]`); the repository convention
   (matching the existing `GcAnalysisResult` and the sibling `CalculateGcContent`) is
   percentage, and tests pin the exact percentage so the units are locked.
2. **Window "variability" = population variance (÷N), not sample variance (÷N−1)** — the
   checklist names "variability" without an estimator; population variance is the natural
   choice because the windows **are** the entire population of windows (not a sample), matching
   the Cuemath definition. Sample variance would change the value.

Recommended coverage (from the artifact): MUST — OverallGcContent `(G+C)/total×100`,
OverallGcSkew `(G−C)/(G+C)`, OverallAtSkew `(A−T)/(A+T)`, GcSkewVariance/GcContentVariance =
population variance of the windowed values (6.8 anchor + `GGCC` windowed dataset), and the
windowed lists' count + `WindowStart`/`WindowEnd`/`Position`. SHOULD — sequence shorter than
window ⇒ empty lists / variances 0 / scalars still computed; null ⇒ `ArgumentNullException`,
null/empty string ⇒ zero/empty. COULD — string and `DnaSequence` overloads agree (delegation).
No source contradictions — GC-content (Brock), GC skew (Lobry/Grigoriev), Biopython, and the
Cuemath population-variance anchor are mutually consistent on formula, range, and estimator.
