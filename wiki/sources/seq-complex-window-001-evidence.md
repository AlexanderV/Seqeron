---
type: source
title: "Evidence: SEQ-COMPLEX-WINDOW-001 (windowed sequence complexity profile)"
tags: [validation, analysis]
doc_path: docs/Evidence/SEQ-COMPLEX-WINDOW-001-Evidence.md
sources:
  - docs/Evidence/SEQ-COMPLEX-WINDOW-001-Evidence.md
source_commit: 177b19d39600766e767f4acb04f34eca05bb79d7
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SEQ-COMPLEX-WINDOW-001

The validation-evidence artifact for test unit **SEQ-COMPLEX-WINDOW-001** — the **windowed
sequence complexity** operation: a sliding-window scan that emits, per window, a `ComplexityPoint`
carrying both a **Shannon entropy** and a **linguistic complexity** value, producing a per-position
*complexity profile*. It is the **profiling / scanning layer** of the complexity/entropy family
(scalar-per-sequence siblings are Lempel–Ziv, DUST, k-mer entropy). One instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the method itself is
written up on the concept page [[windowed-sequence-complexity-profile]]. This file records the
source trace and worked oracles. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **Troyanskaya, Arbell, Koren, Landau & Bolshoy (2002)** *Bioinformatics* 18(5):679–688
    (primary, rank 1; PMID 12050064) — "Sequence complexity profiles of prokaryotic genomic
    sequences: a fast algorithm for calculating linguistic complexity." Computes **complexity
    profiles** (complexity as a function of position along a genome — precisely this unit's
    sliding-window output); linguistic complexity is driven by counts of distinct subwords (the
    summation form), suffix-tree computable in time linear in genome size.
  - **Linguistic sequence complexity** (Wikipedia, citing Trifonov 1990; Gabrielian & Bolshoy 1999;
    Troyanskaya 2002; rank 4) — vocabulary usage `Uᵢ` = ratio of the actual vocabulary size to the
    maximal possible vocabulary; **maximum possible vocabulary at level i** `Vmax,i = min(4^i, N−i+1)`
    (`4^i` = DNA-alphabet theoretical max, `N−i+1` = positional max for length `N`); value range
    `0 < C < 1`. The repo's per-window LC uses the **summation variant**
    `LC = (Σ Vᵢ) / (Σ Vmax,i)` (its standalone unit SEQ-COMPLEX-001; not yet separately ingested).
  - **Shannon entropy** (Wikipedia, citing Shannon 1948; rank 4) — `H(X) = −Σ p(x) log_b p(x)`,
    base 2 ⇒ **bits**; maximum `log_b(n)` for a uniform distribution over `n` outcomes (DNA n=4 ⇒
    max **log₂4 = 2**); deterministic ⇒ **0**; convention `0·log 0 = 0`. Primary: Shannon 1948,
    *Bell System Technical Journal*.
- **Datasets (hand-derived oracles):**

  | Window | Shannon (bits) | Σ Vᵢ | Σ Vmax,i | LinguisticComplexity |
  |--------|----------------|------|----------|----------------------|
  | `ACGTACGT` (uniform, A=C=G=T=2) | **2.0** (log₂4) | 23 | 29 | **23/29 = 0.7931034482758621** |
  | `AAAAAAAA` (poly-A homopolymer) | **0.0** (deterministic) | 6 | 29 | **6/29 = 0.20689655172413793** |

  LC per-length tables (maxWordLength = min(6, windowSize) = 6) given in the artifact; `Vmax,i =
  min(4^i, 8−i+1)`.
- **Window-enumeration geometry** (L = 24, w = 8, s = 8): number of windows = `floor((24−8)/8)+1 =
  **3**`; window 0 = [0,7] center 4, window 1 = [8,15] center 12, window 2 = [16,23] center 20
  (Position = `WindowStart + windowSize/2`, integer division; 0-based, end-inclusive).
- **Documented corner cases / failure modes:** uniform window ⇒ Shannon 2.0; homopolymer window ⇒
  Shannon 0; highly repetitive window ⇒ low LC (few distinct subwords vs the maximum); the sliding
  driver only emits windows **fully contained** in the sequence (`i + w ≤ L`) — a trailing partial
  fragment is **not** emitted, so `L < w` ⇒ empty profile; null DnaSequence ⇒ ArgumentNullException;
  windowSize < 1 and stepSize < 1 ⇒ ArgumentOutOfRangeException.
- **Recommended coverage:** window count `floor((L−w)/s)+1` for `L ≥ w` and 0 when `L < w`; Shannon
  = 2.0 (uniform) / 0.0 (homopolymer); LC = 23/29 (`ACGTACGT`) / 6/29 (poly-A); per-window
  WindowStart/WindowEnd/Position coordinates; argument-validation throws; invariants
  `0 ≤ Shannon ≤ 2` and `0 ≤ LC ≤ 1`.

## Deviations and assumptions

Three **ASSUMPTIONs**, all flagged in the artifact and all non-correctness-affecting for the
complexity *values*:

1. **Center-position convention** — reported `Position = WindowStart + windowSize/2` (integer
   division). No cited source mandates the center label (sources specify the profile *values*, not
   the position label); it is the documented `ComplexityPoint` contract.
2. **Default windowSize = 64, stepSize = 10** — follow common low-complexity windowing practice but
   are not fixed by a single cited source; caller-overridable, and do not affect the value computed
   for an explicitly specified window.
3. **Per-window LC uses maxWordLength = min(6, windowSize)** — the cap of 6 mirrors Gabrielian &
   Bolshoy's "limit vocabulary assessment to W rather than all N−1 values" efficiency choice; the
   specific cap is a repository parameter consistent with the existing LC unit.

## Contradictions

No source contradictions — Troyanskaya 2002 (complexity profile + linguistic complexity), the
Wikipedia linguistic-complexity `Uᵢ = Vᵢ/min(4^i, N−i+1)` definition, and Shannon 1948 (uniform max
log₂4 = 2, deterministic ⇒ 0) are mutually consistent; the hand-derived `ACGTACGT` and poly-A
oracles follow directly from those definitions.
