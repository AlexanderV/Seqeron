---
type: source
title: "Evidence: EPIGEN-CHROM-001 (Chromatin state prediction from histone marks)"
tags: [validation, epigenetics]
doc_path: docs/Evidence/EPIGEN-CHROM-001-Evidence.md
sources:
  - docs/Evidence/EPIGEN-CHROM-001-Evidence.md
source_commit: ed84544abdf23b804d6718e7dfbccaf97bffbe4e
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: EPIGEN-CHROM-001

The validation-evidence artifact for test unit **EPIGEN-CHROM-001** — **chromatin state prediction from
histone modification marks** (ChromHMM-style, `PredictChromatinState` +
`AnnotateHistoneModifications` + `FindAccessibleRegions`). This is the **third ingested unit of the
Epigenetics family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The algorithm is synthesized in its own
concept, [[chromatin-state-prediction]]; [[test-unit-registry]] tracks the unit. Siblings:
[[epigenetic-age-horvath-clock]] and [[bisulfite-methylation-calling]] (both DNA-methylation units —
this one instead consumes histone-mark ChIP-seq signals).

## What this file records

- **Online sources (mutually consistent, no contradictions):**
  - **Ernst & Kellis (2012)** "ChromHMM: automating chromatin-state discovery and characterization",
    *Nature Methods* 9(3):215–216 (authority rank 1, primary) — states are the "major re-occurring
    combinatorial and spatial patterns of marks" capturing promoters/enhancers/transcribed/repressed/
    repetitive regions; the multivariate HMM "explicitly models the **presence or absence** of each
    chromatin mark."
  - **ChromHMM software/manual** (compbio.mit.edu/ChromHMM, rank 3) — the **binary mark model**: raw
    signal is converted to present/absent (1/0) calls per mark via `BinarizeBed`/`BinarizeBam`, then
    `LearnModel` operates on the binary calls (justifies treating state as a function of binarized
    marks, not magnitudes).
  - **Roadmap Epigenomics** chromatin state learning (rank 2, consortium standard) — **core 15-state**
    marks = {H3K4me3, H3K4me1, H3K36me3, H3K27me3, H3K9me3}; **expanded 18-state** adds **H3K27ac**
    (exactly the six marks `PredictChromatinState` takes); the state → characteristic-mark mapping (TssA
    → H3K4me3, Tx/TxWk → H3K36me3, Enh/EnhG → H3K4me1, Het → H3K9me3, TssBiv → H3K4me3+H3K27me3, EnhBiv →
    H3K4me1+H3K27me3, ReprPC → H3K27me3, Quies → none); H3K27ac subdivides active vs weak enhancers/TSS.
  - **Per-mark primaries** (Wikipedia rank 4, each citing a primary): **H3K4me3** active promoters/TSS
    (Liang 2004, PNAS), **H3K4me1** active/primed enhancers (Rada-Iglesias 2018, Nat Genet), **H3K27ac**
    active-enhancer mark separating active from poised enhancers (Creyghton 2010, PNAS), **H3K27me3**
    Polycomb/PRC2 repression (Ferrari 2014, Mol Cell), **H3K9me3** heterochromatin (Nicetto 2019,
    Science), **H3K36me3** transcribed gene bodies (Kimura 2013, J Hum Genet review).

- **Documented corner cases / failure modes:** no mark present → **LowSignal** (Roadmap Quies);
  co-occurring active + repressive marks (H3K4me3 + H3K27me3) = **bivalent** signature (TssBiv), not a
  contradiction; **combinatorial precedence** — a promoter signature (H3K4me3) takes precedence over an
  enhancer signature (H3K4me1) at the same locus; **magnitude is not ordinal beyond the threshold** —
  once a mark is present, more magnitude does not change the state (state = function of the set of
  present marks); negative/zero signals treated as absent.

- **Dataset (canonical Roadmap present/absent → state oracles):** H3K4me3 (±H3K27ac) →
  ActivePromoter (TssA); H3K4me1+H3K27ac → ActiveEnhancer; H3K4me1 no H3K27ac → WeakEnhancer; H3K36me3
  → Transcribed (Tx); H3K27me3 alone → Repressed (ReprPC); H3K9me3 alone → Heterochromatin (Het);
  H3K4me3+H3K27me3 → BivalentPromoter (TssBiv); H3K4me1+H3K27me3 → BivalentEnhancer (EnhBiv); none →
  LowSignal (Quies).

- **Test-coverage recommendations:** MUST — each canonical signature maps to its Roadmap state; bivalent
  (H3K4me3+H3K27me3) → BivalentPromoter (not ActivePromoter/Repressed); no mark → LowSignal; binary
  invariance (same pattern, different magnitudes → same state). SHOULD — `AnnotateHistoneModifications`
  labels each region by its single mark's state; `FindAccessibleRegions` merges contiguous
  above-threshold positions and excludes sub-`minWidth` regions. COULD — negative/zero signals treated
  as absent.

## Deviations and assumptions

- **ASSUMPTION (presence-call threshold value):** ChromHMM binarizes with a Poisson background model
  from raw read counts; a single fixed numeric threshold on an already-normalized [0, 1] signal is not
  specified by the sources. The implementation exposes the presence threshold as a caller-supplied
  parameter (**default 0.5**) and documents it; the state-assignment logic given the present/absent
  pattern is fully source-backed and is what the tests verify. Tests choose magnitudes unambiguously
  above/below the call so results do not depend on the exact default.
- **ASSUMPTION (promoter-over-enhancer precedence):** when H3K4me3 and H3K4me1 co-occur without
  repressive marks the locus is classified as **promoter** (H3K4me3 dominates), consistent with TSS
  states ranking above enhancer states in the Roadmap mnemonic ordering. Marked an assumption because
  Roadmap derives this from spatial HMM context, not a single-locus rule.

No source contradictions.
