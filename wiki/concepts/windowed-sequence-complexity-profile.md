---
type: concept
title: "Windowed sequence complexity profile (sliding-window Shannon + linguistic complexity)"
tags: [analysis, algorithm]
sources:
  - docs/Evidence/SEQ-COMPLEX-WINDOW-001-Evidence.md
source_commit: 177b19d39600766e767f4acb04f34eca05bb79d7
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: seq-complex-window-001-evidence
      evidence: "Test Unit ID: SEQ-COMPLEX-WINDOW-001 ... Algorithm: Windowed Sequence Complexity (sliding-window complexity profile)"
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
[[algorithm-validation-evidence]] for the artifact pattern.

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
distinct subwords relative to the maximum ⇒ low LC. This is the same LC as the standalone unit
SEQ-COMPLEX-001 (its docs `Linguistic_Complexity.md §2.2`; not yet separately ingested), applied
per window.

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

## Worked oracles

- **`ACGTACGT`** window (length 8, `A=C=G=T=2`): Shannon `2.0`; `Σ Vᵢ = 23`, `Σ Vmax,i = 29` ⇒
  **LC = 23/29 = 0.7931034482758621**.
- **`AAAAAAAA`** poly-A window: Shannon `0.0`; `Σ Vᵢ = 6`, `Σ Vmax,i = 29` ⇒
  **LC = 6/29 = 0.20689655172413793**.
- **Geometry** `L = 24, w = 8, s = 8` ⇒ **3** windows: [0,7] center 4, [8,15] center 12, [16,23]
  center 20.

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
