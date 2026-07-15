---
type: concept
title: "Chromatin state prediction (ChromHMM-style histone-mark annotation)"
tags: [epigenetics, algorithm]
mcp_tools:
  - annotate_histone_modifications
  - find_accessible_regions
  - predict_chromatin_state
sources:
  - docs/Evidence/EPIGEN-CHROM-001-Evidence.md
  - docs/algorithms/Epigenetics/Chromatin_State_Prediction.md
source_commit: ed84544abdf23b804d6718e7dfbccaf97bffbe4e
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: epigen-chrom-001-evidence
      evidence: "Test Unit ID: EPIGEN-CHROM-001 ... Algorithm: Chromatin State Prediction from histone modification marks (Epigenetics family)"
      confidence: high
      status: current
---

# Chromatin state prediction (ChromHMM-style histone-mark annotation)

Annotating a genomic locus with its **chromatin state** from the **combination of histone
modification marks** present there, following the **ChromHMM** model (Ernst & Kellis 2012) and the
**Roadmap Epigenomics** 15-/18-state definitions. This is the **third ingested unit of the Epigenetics
family** and a **genuinely distinct algorithm** from its siblings [[epigenetic-age-horvath-clock]] (age
from β-values), [[bisulfite-methylation-calling]] (methylation calling), and
[[cpg-island-detection]] (sequence-only CpG site / island detection) — it operates on **histone
ChIP-seq mark signals**, not on DNA methylation. Validated under test unit **EPIGEN-CHROM-001**; the
record is [[epigen-chrom-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern.

## Binary (present/absent) mark model

ChromHMM "is based on a multivariate Hidden Markov Model that explicitly models the **presence or
absence** of each chromatin mark." Before learning states, raw signal is **binarized** to present/absent
(1/0) calls per mark (ChromHMM `BinarizeBed`/`BinarizeBam`); `LearnModel` then operates on those binary
calls. Chromatin state is therefore a function of the **set of present marks**, not their magnitudes.
Seqeron's `EpigeneticsAnalyzer.PredictChromatinState` takes the **six marks of the Roadmap 18-state
model** — H3K4me3, H3K4me1, **H3K27ac**, H3K36me3, H3K27me3, H3K9me3 — as normalized [0, 1] enrichment
values, calls each present when it exceeds a **presence threshold** (default 0.5, caller-supplied), and
maps the present/absent pattern to a state.

**Binary invariance** (a load-bearing invariant): once a mark is above the call, increasing its
magnitude does **not** change the state — two inputs with the **same present/absent pattern yield the
same state** (ChromHMM binary model).

## State assignment (Roadmap combinatorial signatures)

Chromatin states are the "major re-occurring combinatorial and spatial patterns of marks" and "capture
known classes of genomic elements such as promoters, enhancers, transcribed, repressed, and repetitive
regions." The present-mark pattern maps to a state (Roadmap mnemonic in parentheses):

| Present marks (above call) | State | Roadmap | Characteristic mark |
|----------------------------|-------|---------|---------------------|
| H3K4me3 (± H3K27ac) | **ActivePromoter** | TssA | H3K4me3 (active-promoter/TSS mark) |
| H3K4me1 **+** H3K27ac | **ActiveEnhancer** | active Enh | H3K4me1 + H3K27ac |
| H3K4me1, **no** H3K27ac | **WeakEnhancer** | poised/weak Enh | H3K4me1 |
| H3K36me3 | **Transcribed** | Tx / TxWk | H3K36me3 (transcribed gene body) |
| H3K27me3 (alone) | **Repressed** | ReprPC | H3K27me3 (Polycomb/PRC2) |
| H3K9me3 (alone) | **Heterochromatin** | Het | H3K9me3 |
| H3K4me3 **+** H3K27me3 | **BivalentPromoter** | TssBiv | H3K4me3 + H3K27me3 |
| H3K4me1 **+** H3K27me3 | **BivalentEnhancer** | EnhBiv | H3K4me1 + H3K27me3 |
| none | **LowSignal** | Quies | no enrichment |

Per-mark biology (each a Wikipedia article citing a primary): H3K4me3 = active promoters near TSS
(Liang 2004); H3K4me1 = active/primed enhancers (Rada-Iglesias 2018); H3K27ac = **active** enhancer mark
that "separates active from poised enhancers" (Creyghton 2010) — present on an enhancer marks it active,
absent leaves it weak/poised; H3K27me3 = Polycomb repression (Ferrari 2014); H3K9me3 = heterochromatin
(Nicetto 2019); H3K36me3 = actively transcribed gene bodies (Kimura 2013).

## Two combinatorial rules

- **Bivalency, not contradiction:** co-occurring active + repressive marks are a *state of their own*,
  not an error. H3K4me3 **+** H3K27me3 is the canonical bivalent/poised TSS signature (Roadmap TssBiv) →
  **BivalentPromoter**, classified as bivalent rather than as either active or repressed alone;
  likewise H3K4me1 + H3K27me3 → **BivalentEnhancer** (EnhBiv).
- **Promoter dominates enhancer:** when both active marks H3K4me3 (promoter) and H3K4me1 (enhancer)
  co-occur at one locus without repressive marks, the locus is classified as **promoter** (H3K4me3
  wins), matching the Roadmap TSS-vs-enhancer distinction and the mnemonic ordering (TSS states rank
  above enhancer states). Marked an **assumption** because Roadmap derives this from spatial HMM
  context, not a single-locus rule.

## Companion operations

- **`AnnotateHistoneModifications`** — labels each supplied region by the state of its **single**
  characteristic mark (per-mark delegation / identity mapping), rather than a full combinatorial call.
- **`FindAccessibleRegions`** — ATAC-seq-like peak calling: merges **contiguous above-threshold**
  positions into one accessible region and **excludes** regions narrower than `minWidth`.

## Invariants and edge cases

- **INV (binary invariance):** identical present/absent pattern ⇒ identical state, regardless of the
  magnitudes above the call.
- **INV:** no mark present ⇒ **LowSignal** (Quies).
- **INV:** bivalent signature (H3K4me3 + H3K27me3) ⇒ **BivalentPromoter** — never ActivePromoter or
  Repressed.
- Negative or zero mark signals are treated as **absent** (below any positive call).
- `FindAccessibleRegions` merges adjacent above-threshold positions and drops sub-`minWidth` peaks.

## Scope and limitations

A [[research-grade-limitations|research-grade]] correctness reference for the **state-assignment logic
given a present/absent mark pattern**, which is fully source-backed, **not** a trained multivariate HMM.
The implementation does **not** learn state emission/transition probabilities from data (no `LearnModel`
step), does not perform ChromHMM's Poisson-background binarization from raw read counts, and does not
model spatial (neighbouring-bin) context. The **presence-call threshold** is exposed as a caller
parameter (default 0.5 on an already-normalized [0, 1] signal) rather than derived from a background
model — an assumption; tests choose mark magnitudes unambiguously above/below the call so results do not
depend on the exact default. For production chromatin-state segmentation use ChromHMM or Segway. No
source contradictions — ChromHMM (Ernst & Kellis 2012), the ChromHMM manual, Roadmap Epigenomics, and
the per-mark primaries are mutually consistent.
