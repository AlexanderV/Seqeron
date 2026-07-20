---
type: concept
title: "Windowed sequence complexity profile (sliding-window Shannon + linguistic complexity)"
tags: [analysis, algorithm]
mcp_tools:
  - windowed_complexity
sources:
  - docs/algorithms/Complexity/Windowed_Complexity.md
  - docs/Evidence/SEQ-COMPLEX-WINDOW-001-Evidence.md
  - docs/Evidence/SEQ-ENTROPY-PROFILE-001-Evidence.md
  - docs/algorithms/Statistics/Entropy_Profile.md
  - docs/algorithms/Sequence_Composition/Linguistic_Complexity.md
source_commit: dd6f9c3fb21684add3c59107cae5ea989bbd3315
created: 2026-07-10
updated: 2026-07-17
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: seq-complex-window-001-evidence
      evidence: "Test Unit ID: SEQ-COMPLEX-WINDOW-001 ... Algorithm: Windowed Sequence Complexity (sliding-window complexity profile)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: seq-entropy-profile-001-evidence
      evidence: "Test Unit ID: SEQ-ENTROPY-PROFILE-001; Algorithm: Shannon Entropy Profile (sliding-window per-symbol Shannon entropy) — the same per-window Shannon H=−Σpᵢlog₂pᵢ over single-base composition, emitted as a standalone entropy profile (second entry point)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:entropy-profile
      source: seq-entropy-profile-001-evidence
      evidence: "The standalone Shannon entropy profile SEQ-ENTROPY-PROFILE-001 is homed on entropy-profile: SequenceStatistics.CalculateEntropyProfile (general-alphabet, all letters, W=50/step=1, IEnumerable<double>) — same per-window Shannon statistic, distinct method from this page's DNA-canonical CalculateWindowedComplexity ComplexityPoint Shannon channel (A/T/G/C, W=64/step=10)."
      confidence: high
      status: current
---

# Windowed sequence complexity profile (sliding-window Shannon + linguistic complexity)

A **sliding-window scan** that walks fixed-size windows across a sequence and, for each window,
emits a `ComplexityPoint` carrying two complexity values — a **Shannon entropy** and a **linguistic
complexity** — together with the window's coordinates. The output is a **complexity profile**:
complexity as a function of position along the sequence (Troyanskaya et al. 2002), the raw signal a
caller thresholds to locate low-complexity regions. Validated as test unit
**SEQ-COMPLEX-WINDOW-001**; the source trace and worked oracles are in
[[seq-complex-window-001-evidence]], [[test-unit-registry]] tracks the unit, and see
[[algorithm-validation-evidence]] for the artifact pattern. The **same per-window Shannon scan** is
also exposed standalone (entropy profile only, no linguistic complexity) as a second entry point,
**SEQ-ENTROPY-PROFILE-001** — see the [second-entry-point section](#second-entry-point-the-standalone-shannon-entropy-profile-seq-entropy-profile-001) below and [[seq-entropy-profile-001-evidence]].

## Where it sits in the complexity family — the *profiling* layer

This is a **distinct operation** from its `SEQ-COMPLEX-*` siblings, not a re-derivation of one:
every sibling reduces a **whole sequence to a single scalar**, whereas this unit is the
**sliding-window driver** that produces a **per-position profile** and is the family's first home
of **linguistic complexity**:

- vs. **compression complexity** — [[sequence-complexity-compression-lempel-ziv]] (Lempel–Ziv LZ76)
  emits one normalized scalar for the whole string. This unit re-scores *local windows* and reports
  each window's value along the sequence.
- vs. **DUST triplet score** — [[dust-low-complexity-score]] is the DNA low-complexity **masker**:
  one triplet-frequency scalar per window/sequence, a *high* score meaning *low* complexity. The
  windowed profile is a companion *scanning* view; a caller masks by thresholding the emitted
  profile, but the per-window statistics here are Shannon entropy + linguistic complexity, not the
  `∑ c(c−1)/2` DUST statistic.
- vs. **k-mer k-entropy** — [[k-mer-statistics]]'s entropy field / `CalculateKmerEntropy`
  (SEQ-COMPLEX-KMER-001) is a fixed-`k` Shannon entropy of the k-mer distribution over the whole
  sequence. The per-window Shannon here is over the **single-base composition** of a window (n = 4),
  and it is paired with linguistic complexity.
- vs. **protein SEG** — [[protein-low-complexity-seg]] is the protein-side sliding-window
  low-complexity detector (Shannon entropy of amino-acid composition). This is the DNA analogue's
  *profiling* form, adding linguistic complexity as a second per-window statistic.

## The two per-window measures

For a window `x` of length `w`:

**Shannon entropy of composition** — `H = −Σ pᵦ log₂ pᵦ` over the four base frequencies `pᵦ` in the
window, in **bits** (Shannon 1948; convention `0·log 0 = 0`). A **uniform** window (all four bases
equally frequent) gives the DNA maximum `log₂4 = 2.0`; a **homopolymer** window (one base,
deterministic) gives `0`. Invariant `0 ≤ H ≤ 2`.

**Linguistic complexity (summation variant)** — vocabulary usage aggregated over word lengths:

```
LC = ( Σᵢ Vᵢ ) / ( Σᵢ Vmax,i )      i = 1 … maxWordLength
Vmax,i = min(4^i, N − i + 1)
```

where `Vᵢ` = number of **distinct** length-`i` subwords observed in the window, `N` = window length,
`4^i` = the DNA-alphabet theoretical maximum vocabulary, and `N−i+1` = the positional maximum
(Gabrielian & Bolshoy 1999; Troyanskaya 2002). `maxWordLength = min(6, windowSize)` (the cap of 6 is
a repository efficiency parameter, ASSUMPTION). Range `0 < LC ≤ 1`; a **repetitive** window has few
distinct subwords relative to the maximum ⇒ low LC. This is the same LC as the standalone scalar unit
SEQ-COMPLEX-001 — [[linguistic-complexity]] (`SequenceComplexity.CalculateLinguisticComplexity`) —
applied per window via the shared `CalculateLinguisticComplexityCore`.

## Window geometry and the `ComplexityPoint` contract

The driver enumerates windows fully contained in the sequence and steps by `stepSize`:

| Quantity | Value |
|----------|-------|
| Number of windows | `floor((L − w)/s) + 1` for `L ≥ w`; **0** when `L < w` |
| WindowStart / WindowEnd | 0-based, **end-inclusive** (`[start, start+w−1]`) |
| Position (center) | `WindowStart + windowSize/2` (integer division; ASSUMPTION — label only) |

Only windows with `i + w ≤ L` are emitted — a **trailing partial fragment is never emitted**, so a
sequence shorter than one window yields an **empty profile**. Defaults `windowSize = 64`,
`stepSize = 10` (ASSUMPTION — common windowing practice, caller-overridable; do not affect the value
computed for an explicitly specified window).

## Implementation surface

The primary spec (`docs/algorithms/Complexity/Windowed_Complexity.md`) exposes the driver as
`SequenceComplexity.CalculateWindowedComplexity(DnaSequence, int windowSize = 64, int stepSize = 10)`
on `Seqeron.Genomics.Analysis` (`SequenceComplexity.cs`), returning a **lazily-evaluated**
`IEnumerable<ComplexityPoint>`. It **delegates per-window metrics** to the same
`CalculateShannonEntropyCore` / `CalculateLinguisticComplexityCore` helpers used by the standalone
scalar units, so window values match those metrics exactly. Cost is `O((L/s)·w²)` time — one pass over
≈`L/s` windows, each enumerating subwords of length `1…min(6,w)` — with per-window space bounded by
the distinct-subword count. A **suffix tree was evaluated and rejected**: this is a scoring scan over
each window rather than an exact-match occurrence search, and the small LC word-length cap (≤6) keeps
direct enumeration cheap; the linear-time suffix-tree profile of Troyanskaya et al. (2002) is *not*
implemented — the exact per-window enumeration is authoritative for the bounded word lengths used.

## Worked oracles

- **`ACGTACGT`** window (length 8, `A=C=G=T=2`): Shannon `2.0`; `Σ Vᵢ = 23`, `Σ Vmax,i = 29` ⇒
  **LC = 23/29 = 0.7931034482758621**.
- **`AAAAAAAA`** poly-A window: Shannon `0.0`; `Σ Vᵢ = 6`, `Σ Vmax,i = 29` ⇒
  **LC = 6/29 = 0.20689655172413793**.
- **Geometry** `L = 24, w = 8, s = 8` ⇒ **3** windows: [0,7] center 4, [8,15] center 12, [16,23]
  center 20.

## Sibling: the standalone Shannon entropy profile (SEQ-ENTROPY-PROFILE-001 → [[entropy-profile]])

A **sibling sliding-window profiler** validated as test unit **SEQ-ENTROPY-PROFILE-001** emits, per
window, the **entropy value alone** (`H = −Σ pᵢ log₂ pᵢ` in **bits**), *without* the paired
linguistic complexity — now homed on its own page [[entropy-profile]]. It is the **same statistic**
(per-window Shannon `H` over single-symbol, `k = 1` composition), but a **distinct method**: it is
`SequenceStatistics.CalculateEntropyProfile` (defaults `W = 50`, `step = 1`, returning
`IEnumerable<double>`) over the **general alphabet** — counting **all letters**, so `N`/degenerate
symbols count as themselves and protein windows can exceed `2` bits — whereas this page's
`ComplexityPoint` Shannon channel is the **DNA-canonical** `SequenceComplexity` path (A/T/G/C only,
defaults `W = 64`, `step = 10`). Convention `0·log 0 = 0`; DNA bounds `0 ≤ H ≤ log₂k ≤ 2` (uniform
4-base window ⇒ `2.0`, homopolymer ⇒ `0`). Full contract, invariants, and worked oracles on
[[entropy-profile]]; source trace in [[seq-entropy-profile-001-evidence]]; [[test-unit-registry]]
tracks the unit.

This is the **character-level (k = 1) counterpart** of the k-mer/block Shannon **k-entropy**
([[k-mer-statistics]] `AnalyzeKmers.Entropy` / SEQ-COMPLEX-KMER-001 `CalculateKmerEntropy`), which is
the same formula over the `L − k + 1` overlapping k-mer distribution of a **whole sequence**; this
unit instead scans **per-position windows** and uses the mono-nucleotide (n = 4) alphabet.

Worked per-window oracles (bits): `AAAA` → **0.0** (homopolymer); `AATT` → **1.0**; `ATGC` → **2.0**
(= log₂4, uniform); `AAAT` → **0.8112781244591328** (3:1 skew); `AATG`/`GCAA` → **1.5**; `AAATTC` →
**1.4591479170272448**. Sliding profiles: `AAATGC` w = 4 s = 1 → `[0.8112781244591328, 1.5, 2.0]`;
`AAATGCAA` w = 4 s = 2 → `[0.8112781244591328, 2.0, 1.5]`. Corner cases mirror the profile driver:
`windowSize > length` ⇒ **empty** profile (no partial trailing window), `windowSize == length` ⇒ a
single value; case-folded before counting (case-insensitive). One non-value-affecting ASSUMPTION —
the mono-symbol (k = 1) alphabet is the implementation's modelling choice; the cited sources define
`H` over any symbol distribution.

## Corner cases and contract

- **`L < w`** ⇒ empty profile (no partial trailing window).
- **Bounds:** `0 ≤ ShannonEntropy ≤ 2`, `0 ≤ LinguisticComplexity ≤ 1` across all windows.
- **Argument validation:** null `DnaSequence` ⇒ `ArgumentNullException`; `windowSize < 1` and
  `stepSize < 1` ⇒ `ArgumentOutOfRangeException`.

## References

Shannon C.E. (1948) *Bell System Technical Journal* 27(3):379–423; Troyanskaya O.G. et al. (2002)
*Bioinformatics* 18(5):679–688 (PMID 12050064); Gabrielian A. & Bolshoy A. (1999) *Computers &
Chemistry* 23(3–4):263–274; Trifonov E.N. (1990) *Making sense of the human genome*. Full citations
in [[seq-complex-window-001-evidence]] (do not duplicate here). A
[[research-grade-limitations|research-grade]] implementation combining standard Shannon entropy and
linguistic complexity over a sliding window.
