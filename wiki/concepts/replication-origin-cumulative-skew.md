---
type: concept
title: "Replication origin/terminus via cumulative GC-skew extrema"
tags: [sequence-statistics, composition, chromosome]
mcp_tools:
  - cumulative_gc_skew
  - predict_replication_origin
sources:
  - docs/Evidence/SEQ-REPLICATION-001-Evidence.md
  - docs/algorithms/Sequence_Composition/Replication_Origin_Prediction.md
source_commit: 52fb606581c4d8e490ad414659bcccf9e72c9d9a
created: 2026-07-10
updated: 2026-07-17
graph:
  relationships:
    - predicate: relates_to
      object: concept:nucleotide-composition-skew
      source: seq-replication-001-evidence
      evidence: "SEQ-REPLICATION-001 integrates the per-base GC skew into a cumulative diagram whose sign flip at origin/terminus is the strand asymmetry nucleotide-composition-skew defines (Lobry 1996)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: seq-replication-001-evidence
      evidence: "Test Unit ID: SEQ-REPLICATION-001; Algorithm: Replication Origin Prediction (cumulative GC-skew minimum)"
      confidence: high
      status: current
---

# Replication origin/terminus via cumulative GC-skew extrema

**Replication origin prediction** (test unit **SEQ-REPLICATION-001**,
[[seq-replication-001-evidence]]) is the classic ori-finding method: integrate the strand
skew along a chromosome into a **cumulative diagram**, then read the **origin off its global
minimum** and the **terminus off its global maximum**. It is the *locating* algorithm that the
scalar and windowed skew pages ([[nucleotide-composition-skew]],
[[windowed-gc-profile-and-variance]]) both flagged as the biological payoff of skew — realised
here as a distinct algorithm with its own outputs (`PredictedOrigin`, `PredictedTerminus`,
`IsSignificant`).

## The cumulative skew — an integer running sum

Unlike the ratio-normalised `(G−C)/(G+C)` scalar of [[nucleotide-composition-skew]], the
locator uses the **integer cumulative skew** of Rosalind's *Minimum Skew Problem* (BA1F):

- Walk the sequence left to right accumulating a running total: **G contributes +1, C
  contributes −1, A and T contribute 0**, starting at **`Skew_0 = 0`** before any base.
- This yields **`|Genome| + 1` prefix values** `Skew_0 … Skew_n`, indexed by prefix position
  `i ∈ [0, |Genome|]` (a *boundary* index between bases, not a base index).

Worked diagrams: `CCGGGG` → `0,−1,−2,−1,0,+1,+2`; `GGGCCC` → `0,+1,+2,+3,+2,+1,0`;
`AATT` → `0,0,0,0,0` (flat — A/T never move the diagram).

## Locating rule — argmin = origin, argmax = terminus

- **`PredictedOrigin` = position of the global minimum** of the cumulative diagram
  (Grigoriev 1998; Wikipedia GC skew: *"the minimum value corresponds to the origin"*).
- **`PredictedTerminus` = position of the global maximum** (*"the maximum value corresponds to
  the terminal"*). On a real bacterial chromosome the two extrema sit **≈ half the chromosome
  length apart** (Grigoriev 1998).
- **Tie-break:** when several positions share the extreme value, return the **first/smallest**
  minimizing index. The Rosalind BA1F sample (100 nt) has **two** minimizers, `53` and `97`
  (min value −4); this unit returns `PredictedOrigin = 53` — the in-session re-derivation
  reproduces the published `53 97` exactly.

The biological basis: the leading strand is **G-rich** and the lagging strand C-rich, so skew
holds one sign along each replichore and **flips at the origin and terminus** (Lobry 1996). The
running integral of that per-base sign therefore turns downward until the origin and upward
after it, making the origin a valley and the terminus a peak.

## `IsSignificant` — a threshold-free amplitude check

There is **no authoritative numeric significance cutoff** for an origin call. An earlier
implementation used an invented constant (`amplitude > count × 0.01`), now **removed** as
untraceable. `IsSignificant` is redefined as the weakest non-invented predicate **`max > min`**
— the diagram has non-zero amplitude, i.e. a detectable strand-composition asymmetry exists
(Lobry 1996 / Grigoriev 1998). A **flat diagram** (no G/C, or balanced) gives `max = min` ⇒
origin = terminus = 0 and `IsSignificant = false`. Callers wanting a quantitative confidence
should inspect the skew **amplitude** (`max − min`) directly.

## Implementation contract — entry points, output, invariants

The primary spec (SEQ-REPLICATION-001) fixes the surface in
`Seqeron.Genomics.Analysis/GcSkewCalculator.cs` (status **Production**):

- **`GcSkewCalculator.PredictReplicationOrigin(DnaSequence)`** — canonical method; throws
  `ArgumentNullException` on null.
- **`GcSkewCalculator.PredictReplicationOrigin(string)`** — thin overload; upper-cases and
  delegates to the same core, returning a **zero prediction** (`PredictedOrigin =
  PredictedTerminus = 0`, skews `0`, `IsSignificant = false`) for null/empty input.

There is **no window or threshold parameter** — the only constants are the per-base increments
G:+1 / C:−1 / A,T/other:0. The output is the record struct **`ReplicationOriginPrediction`**:

| Field | Type | Meaning |
|-------|------|---------|
| `PredictedOrigin` | `int` | first prefix index minimizing the cumulative skew (*ori*) |
| `PredictedTerminus` | `int` | first prefix index maximizing the cumulative skew (*ter*) |
| `OriginSkew` | `double` | skew value at the minimum (≤ 0) |
| `TerminusSkew` | `double` | skew value at the maximum (≥ 0) |
| `IsSignificant` | `bool` | `true` ⇔ `max > min` (non-zero amplitude) |

**Invariants** (INV-01…06): origin = first argmin, terminus = first argmax (strict `<`/`>` keep
the smallest tied index); `OriginSkew ≤ 0 ≤ TerminusSkew` because `Skew_0 = 0` is always a prefix
value; `0 ≤ PredictedOrigin, PredictedTerminus ≤ n`; `IsSignificant ⇔ max > min`; A/T never move
the diagram. Computation is a **single O(n)-time, O(1)-space pass** — the diagram is never
materialized (min/max tracked in scalars), and the repository suffix tree is not applicable
(this is a scalar fold, not substring search).

**Edge cases** (spec §6.1): no G/C (`AAAATTTT`) → origin = terminus = 0, not significant; single
`G` → origin 0 (skew 0), terminus 1 (skew +1); tied extremum → first index; null `DnaSequence` →
`ArgumentNullException`; null/empty `string` → zero prediction. The API deliberately returns a
**single** origin/terminus even when several positions tie for the extremum, whereas Rosalind
BA1F enumerates *all* minimizers.

## Why it matters

This is the first wiki entry where the skew family's headline application — *finding the
replication origin* — is an actual algorithm rather than a motivating footnote. It consumes the
same per-base G/C tally that [[base-composition]] counts and [[nucleotide-composition-skew]]
normalises, but its distinguishing move is **cumulative integration + extremum search**, the
Grigoriev-1998 diagram. The [[windowed-gc-profile-and-variance]] scan is the *variability* view
of the same skew signal; this is its *localisation* view.

## Scope / assumptions

The counting convention (G:+1/C:−1/A,T:0, `Skew_0 = 0`, prefix indices `[0, n]`) and the
min = origin / max = terminus locating rule are **fully sourced** (Rosalind BA1F; Grigoriev
1998; Lobry 1996; corroborated by Wikipedia "GC skew"). Only the **`IsSignificant` predicate**
(`max > min`) is a documented assumption — chosen precisely because no source defines a
significance threshold; it is the evidence-neutral floor, not an invented constant. See
[[seq-replication-001-evidence]].
