---
type: concept
title: "GC skew — the GcSkewCalculator surface (scalar → windowed → cumulative → origin)"
tags: [sequence-statistics, composition, chromosome]
mcp_tools:
  - gc_skew
  - windowed_gc_skew
sources:
  - docs/algorithms/Sequence_Composition/GC_Skew.md
source_commit: 0a33b97be4e29c589fb37eb5beb5c1f629d5b937
created: 2026-07-17
updated: 2026-07-17
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: gc-skew
      evidence: "Test Unit ID: SEQ-GCSKEW-001; Algorithm Group: Sequence Composition; the standalone GC-skew unit flagged by the sibling skew pages, here ingested as the GcSkewCalculator class-level anchor"
      confidence: high
      status: current
    - predicate: specializes
      object: concept:nucleotide-composition-skew
      source: gc-skew
      evidence: "CalculateGcSkew = (G−C)/(G+C) is the GC member of the Lobry-1996 skew family; this page is the whole-sequence scalar entry point + the windowed/cumulative/origin surface of the same GcSkewCalculator, the scalar formula & range being defined on nucleotide-composition-skew"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:windowed-gc-profile-and-variance
      source: gc-skew
      evidence: "CalculateWindowedGcSkew emits the per-window GC-skew series that AnalyzeGcContent (SEQ-GC-ANALYSIS-001) bundles into GcAnalysisResult with its population variance; same recount-per-window cores"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:replication-origin-cumulative-skew
      source: gc-skew
      evidence: "CalculateCumulativeGcSkew + PredictReplicationOrigin are the cumulative-diagram and argmin/argmax origin-locator surface of the same class; kept distinct as the downstream consumer of GC skew"
      confidence: high
      status: current
---

# GC skew — the GcSkewCalculator surface

**GC skew** is the strand-composition statistic `(G − C) / (G + C)`, validated as the standalone
test unit **SEQ-GCSKEW-001** over the `GcSkewCalculator` class (`Seqeron.Genomics.Analysis`,
spec `docs/algorithms/Sequence_Composition/GC_Skew.md`, implementation status *Simplified*). This
is the `gc-skew` unit the sibling skew pages flagged as "separately flagged, not yet ingested" —
it is the **class-level anchor** for the four-tier GC-skew surface: the whole-sequence **scalar**,
the **windowed** profile, the **cumulative** running sum, and the heuristic **origin/terminus**
predictor. The scalar formula and range live on [[nucleotide-composition-skew]]; the windowed
profile + variance on [[windowed-gc-profile-and-variance]]; the cumulative diagram + argmin/argmax
locator on [[replication-origin-cumulative-skew]]. This page carries the **whole-sequence scalar
entry point** and the **cross-overload parameter-validation contract** that spans all four tiers.

## The four entry points

| Method | Tier | Returns | Home for the detail |
|--------|------|---------|---------------------|
| `CalculateGcSkew(sequence)` | scalar | `double` overall `(G−C)/(G+C)` | *this page* + [[nucleotide-composition-skew]] |
| `CalculateWindowedGcSkew(sequence, windowSize=1000, stepSize=100)` | windowed | `IEnumerable<GcSkewPoint>` | [[windowed-gc-profile-and-variance]] |
| `CalculateCumulativeGcSkew(sequence, windowSize=1000)` | cumulative | `IEnumerable<CumulativeGcSkewPoint>` (fixed, **non-overlapping** windows — the routine sets `stepSize = windowSize`) | [[replication-origin-cumulative-skew]] |
| `PredictReplicationOrigin(sequence, windowSize=1000)` | locator | `ReplicationOriginPrediction` (`PredictedOrigin`, `PredictedTerminus`, cumulative-skew values, significance flag) | [[replication-origin-cumulative-skew]] |

The same class also ships `CalculateAtSkew(...)` (the AT member, primary-documented via
`Extended_GC_Skew_Analysis/AT_Skew.md` on [[nucleotide-composition-skew]]) and the composite
`AnalyzeGcContent(...)` (SEQ-GC-ANALYSIS-001, [[windowed-gc-profile-and-variance]]).

## Whole-sequence scalar — `CalculateGcSkew`

A single **O(n) time / O(1) space** pass counts `G` and `C` and returns `(G − C) / (G + C)`:

- **INV-01 — bounded** `−1 ≤ skew ≤ 1` (the numerator is bounded by the denominator in absolute
  value). `+1 ⇔ C = 0`, `−1 ⇔ G = 0` among the G/C pair.
- **INV-02 — zero-denominator ⇒ `0`.** Empty input, or a region with no `G`/`C`, returns `0.0`
  (never NaN or an exception) — the Biopython `GC_skew` convention (see
  [[nucleotide-composition-skew]] for the case-insensitive / ignore-ambiguous counting rules).
- **INV-03 — windowed positions are reported at the window center** `WindowStart + WindowSize / 2`.

The scalar has two overloads: `CalculateGcSkew(string)` returns `0` for empty input, while
`CalculateGcSkew(DnaSequence)` throws `ArgumentNullException` on a null value object.

## The parameter-validation sharp edge (spans all overloads)

The distinctive, easy-to-miss contract of this class is that **window/step validation is
inconsistent across the overloads** — a genuine sharp edge not carried by the sibling pages:

- **Validated:** the *typed* `CalculateWindowedGcSkew(DnaSequence, …)` and
  `CalculateCumulativeGcSkew(DnaSequence, …)` overloads throw **`ArgumentOutOfRangeException`**
  when `windowSize < 1` or `stepSize < 1`.
- **Unvalidated:** the *raw-string* windowed/cumulative overloads only guard empty input and then
  delegate straight to the shared core loops. Nonpositive sizes are **unsupported rather than
  rejected** — in particular **`stepSize = 0` in `CalculateWindowedGcSkew(string, …)`** and
  **`windowSize = 0` in `CalculateCumulativeGcSkew(string, …)`** can produce **non-terminating
  enumerations**.
- **Bypassed:** the higher-level typed helpers `PredictReplicationOrigin(...)` and
  `AnalyzeGcContent(...)` also call the unvalidated cores directly, so they **assume** positive
  `windowSize`/`stepSize` without enforcing them — the validated typed windowed/cumulative entry
  points are *not* on their path. `PredictReplicationOrigin(...)` returns a **zeroed prediction**
  when no cumulative points are available.

Practical rule: pass window/step through the **typed `DnaSequence` overloads** if you want the
`ArgumentOutOfRangeException` guard; treat any nonpositive size on the string / helper paths as
undefined (and potentially hanging) behaviour.

## Scope / assumptions

The core skew formula and the sliding-window / cumulative variants are **fully sourced**
(Lobry 1996; Grigoriev 1998; Tillier & Collins 2000; Wikipedia "GC skew"). The *Simplified*
status attaches only to the **origin/terminus predictor** — it maps the cumulative-skew minimum
directly to the origin and the maximum to the terminus, assumes a single circular chromosome with
bidirectional replication, and does not correct for genome rearrangements or horizontal transfer;
that heuristic (and its significance predicate) is synthesized on
[[replication-origin-cumulative-skew]], kept **distinct from GC skew itself** as the downstream
consumer. See [[test-unit-registry]] for the unit record and [[algorithm-validation-evidence]]
for the artifact pattern.
