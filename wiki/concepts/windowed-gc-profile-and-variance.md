---
type: concept
title: "Windowed GC profile & compositional variance"
tags: [sequence-statistics, composition, chromosome]
sources:
  - docs/Evidence/SEQ-GC-ANALYSIS-001-Evidence.md
  - docs/Evidence/SEQ-GC-PROFILE-001-Evidence.md
  - docs/Evidence/SEQ-REPLICATION-001-Evidence.md
source_commit: c094b65e4a89b3c3c146d655c12489e6d28e8564
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: seq-gc-analysis-001-evidence
      evidence: "Test Unit ID: SEQ-GC-ANALYSIS-001 ... Algorithm: Comprehensive GC Analysis (GC content, GC skew, AT skew, windowed profiles, compositional variance)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: seq-gc-profile-001-evidence
      evidence: "Test Unit ID: SEQ-GC-PROFILE-001; Algorithm: GC Content Profile (sliding-window GC content) — the same windowed GC%=(G+C)/(A+T+G+C)×100 and window geometry as the composite's GC-content channel, emitted standalone (GC% only, no skew profile, no variance)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:base-composition
      source: seq-gc-profile-001-evidence
      evidence: "SEQ-GC-PROFILE-001 windows the GC%=(G+C)/(A+T+G+C)×100 quantity base-composition defines (Wikipedia GC-content; Biopython gc_fraction ×100), N excluded from the denominator (remove mode)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:base-composition
      source: seq-gc-analysis-001-evidence
      evidence: "OverallGcContent (%) = (G+C)/(A+T+G+C)×100 — the same GC-content quantity base-composition covers, here bundled into GcAnalysisResult and windowed"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:nucleotide-composition-skew
      source: seq-gc-analysis-001-evidence
      evidence: "OverallGcSkew = (G−C)/(G+C) and OverallAtSkew = (A−T)/(A+T) — the skew family this GcAnalysisResult computes as scalars and as windowed profiles"
      confidence: high
      status: current
---

# Windowed GC profile & compositional variance

**Comprehensive GC analysis** (`GcAnalysisResult`, test unit **SEQ-GC-ANALYSIS-001**,
[[seq-gc-analysis-001-evidence]]) is the *composite* GC unit: it bundles the whole-sequence
GC/AT scalars with a **sliding-window profile** of GC content and GC skew, plus the
**population variance** of each windowed series. The scalars are not new — this concept is the
home for the two things that are: the **windowed GC profiles** and the **compositional
variance** summaries. The windowed GC%-content channel also ships as a **standalone unit**,
**SEQ-GC-PROFILE-001** ([[seq-gc-profile-001-evidence]]) — GC% only, no skew profile and no
variance; see the [standalone-unit section](#standalone-entry-point-the-gc-only-window-profile-seq-gc-profile-001) below.

## The six outputs

| Output | Formula | Home |
|--------|---------|------|
| `OverallGcContent` (%) | `(G+C)/(A+T+G+C) × 100` | [[base-composition]] |
| `OverallGcSkew` | `(G−C)/(G+C)` | [[nucleotide-composition-skew]] |
| `OverallAtSkew` | `(A−T)/(A+T)` | [[nucleotide-composition-skew]] |
| windowed GC% / GC-skew lists | per-window `(WindowStart, WindowEnd, Position)` | *this concept* |
| `GcContentVariance` | `σ² = Σ(xᵢ−μ)²/N` of the windowed GC% | *this concept* |
| `GcSkewVariance` | `σ² = Σ(xᵢ−μ)²/N` of the windowed GC skew | *this concept* |

The overall scalars are the same quantities [[base-composition]] (GC content) and
[[nucleotide-composition-skew]] (GC/AT skew) already define — this unit re-exposes them inside
one result object rather than re-deriving them. Note the one **units** difference: this unit
reports GC content as a **percentage (×100)** (Brock / Madigan & Martinko convention), whereas
[[base-composition]] carries the Biopython `gc_fraction` **`[0,1]`** value; the two differ only
by the factor 100.

## Windowed profile

A fixed-width window slides along the sequence, emitting a per-window GC% and GC skew tagged
with its `WindowStart`, `WindowEnd`, and center `Position`. Only **fully-contained** windows
are emitted; a sequence **shorter than one window** produces **empty** lists (and therefore
zero window-derived variance), while the overall scalars are still computed over the whole
sequence. This is the same "multiple windows along the sequence" idea Biopython's `GC_skew`
formalizes, and the same zero-denominator convention applies per window (`G+C = 0 ⇒ skew 0`,
non-ACGT symbols ignored — see [[nucleotide-composition-skew]]).

## Compositional variance — population variance of the windows

`GcContentVariance` and `GcSkewVariance` are each the **population variance**
`σ² = Σ(xᵢ−μ)²/N` (÷N, *not* the Bessel-corrected ÷(N−1)) of their windowed series. The ÷N
choice is a **documented assumption**: the windows *are* the entire population of windows over
this sequence, not a sample of a larger set, so population variance is the natural estimator
(Cuemath population-variance definition, anchored by the worked example `{12,13,12,14,19}` →
`6.8`). Worked GC datasets: `GGGCCAT` (whole-sequence window) → GC% `71.428…`, GC skew `0.2`,
AT skew `0.0`; windowing `GGCC` at w=2/s=2 → windows `GG` (skew `+1`, GC% 100) and `CC`
(skew `−1`, GC% 100) ⇒ `GcSkewVariance = ((1)²+(−1)²)/2 = 1.0`, `GcContentVariance = 0.0`.

## Standalone entry point: the GC%-only window profile (SEQ-GC-PROFILE-001)

The windowed **GC-content** channel above has a **standalone entry point** validated as test unit
**SEQ-GC-PROFILE-001** — a **sliding-window GC-content profile** that slides a fixed-width window
along the sequence and emits, per fully-contained window, its **`GC% = (G+C)/(A+T+G+C) × 100`**
*alone*, without the paired GC-skew profile or the population-variance summaries. It is the **same
measure and the same window geometry** as this composite's GC-content channel — not a different
statistic — just that channel emitted on its own. Same conventions apply: **GC% percentage (×100)**
units, **N and other non-standard symbols excluded from the denominator** (matching Biopython's
default `ambiguous="remove"` and the Wikipedia A+T+G+C denominator — window `GGAN` → `2/3×100 =
66.66…%`), **RNA U treated as a non-GC base** (= T), and the **empty-window convention** (a window
with no standard base ⇒ `A+T+G+C = 0` division by zero ⇒ **GC% `0`**, mirroring the sibling
`GcSkewCalculator`). Window count `⌊(n − w)/step⌋ + 1` with offsets `0, step, 2·step, …`; a sequence
**shorter than one window** (or null/empty) ⇒ **empty profile**; every window value is bounded
`[0, 100]`.

Worked window oracles (×100): `GGGG` → **100.0**, `AAAA` → **0.0**, `ATGC` → **50.0**, `GGGA` →
**75.0**, `GCAT` → **50.0**, `GGAN` (N excluded) → **66.66666666666666**; Biopython `gc_fraction`
doctests scaled ×100 — `ACTG` → 50.0, `ACTGN` (remove) → 50.0, `ACTGN` (ignore) → 40.0, RNA
`GGAUCUUCGGAUCU` → 50.0. Source trace and full dataset in [[seq-gc-profile-001-evidence]];
[[test-unit-registry]] tracks the unit. This is the composition analogue of the standalone
[[windowed-sequence-complexity-profile|windowed Shannon-entropy profile]] (SEQ-ENTROPY-PROFILE-001)
— a per-window scalar channel factored out of a richer composite scan.

## Why it matters

The variance summaries quantify **compositional heterogeneity along a sequence** — how much GC
content and strand skew fluctuate window-to-window. This is the scalar counterpart of the
isochore / GC-heterogeneity view of a genome, and the cumulative version of GC skew is what
locates the **replication origin/terminus** (its sign switch; Grigoriev 1998) — the distinct
argmin/argmax locating algorithm synthesized on [[replication-origin-cumulative-skew]]
(SEQ-REPLICATION-001), the biological motivation the skew family carries in
[[nucleotide-composition-skew]]. It is the same
GC-variability signal the centromere heuristic in [[centromere-analysis]] leans on, made
explicit as a windowed profile plus a variance. The per-window entropy/complexity profile of
[[windowed-sequence-complexity-profile]] is the information-theoretic sibling of this
composition-specific windowed scan.

## Scope / assumptions

Formulas are fully sourced (GC content — Brock/Wikipedia; GC/AT skew — Lobry 1996 /
Charneski 2011 / Wikipedia; population variance — Cuemath). Only two **labelling/estimator**
choices are assumed, neither correctness-affecting at the formula level: GC content reported as
a **percentage** (×100, matching the repository `GcAnalysisResult`/`CalculateGcContent`) rather
than the `[0,1]` fraction, and window "variability" taken as **population variance** (÷N)
rather than sample variance (÷N−1). See [[seq-gc-analysis-001-evidence]].
