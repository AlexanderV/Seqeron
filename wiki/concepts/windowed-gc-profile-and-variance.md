---
type: concept
title: "Windowed GC profile & compositional variance"
tags: [sequence-statistics, composition, chromosome]
sources:
  - docs/Evidence/SEQ-GC-ANALYSIS-001-Evidence.md
source_commit: f6fc5f03fffb7fd2053db36d0ad79995b8affe3e
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
variance** summaries.

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

## Why it matters

The variance summaries quantify **compositional heterogeneity along a sequence** — how much GC
content and strand skew fluctuate window-to-window. This is the scalar counterpart of the
isochore / GC-heterogeneity view of a genome, and the cumulative version of GC skew is what
locates the **replication origin/terminus** (its sign switch; Grigoriev 1998), the biological
motivation the skew family carries in [[nucleotide-composition-skew]]. It is the same
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
