---
type: source
title: "Evidence: PROTMOTIF-LC-001 (Protein low-complexity region detection — SEG)"
tags: [validation, protein]
doc_path: docs/Evidence/PROTMOTIF-LC-001-Evidence.md
sources:
  - docs/Evidence/PROTMOTIF-LC-001-Evidence.md
source_commit: dda0efbabe0d4288e2a5f4e50964fb84a531cf23
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PROTMOTIF-LC-001

The validation-evidence artifact for test unit **PROTMOTIF-LC-001** — **low-complexity region
detection** by the **SEG algorithm** (Wootton & Federhen 1993). This is **the same algorithm**
already synthesized under [[protein-low-complexity-seg]] (test unit DISORDER-LC-001): identical
defaults (W=12, K1=2.2, K2=2.5), identical Shannon-entropy-in-bits/residue complexity measure,
and the same two-stage trigger/extend scan. PROTMOTIF-LC-001 is the **ProteinMotif-family
registration** of the SEG unit — a second Evidence file tracing the same method — so this page
records only what it adds and defers the model, oracle table and deviations to the concept.
It is one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence
artifact]] pattern; see [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources** (all consistent with the DISORDER-LC-001 set, some additional):
  - **NCBI `ncbi-seg` man page** (rank 2, official program spec) — trigger window `W` [default
    12], trigger complexity `K1` [default 2.2 bits], extension complexity `K2` [default 2.5
    bits]; complexity "in units of bits", range 0…4.322 = log₂20; "only values greater than K1
    are effective in extending"; complexity "defined by equation (3) of Wootton & Federhen".
  - **NCBI BLAST `blast_seg.c`** (rank 3, reference impl) — verbatim constants `kSegWindow=12`,
    `kSegLocut=2.2`, `kSegHicut=2.5`; `s_Entropy()` Shannon entropy over residue counts
    normalized by log alphabet size; `s_LnPerm()` = "K2 entropy per Wootton & Federhen eq. 3";
    precomputed `lnfact[]` log-factorial table for the exact permutation form; two-stage
    algorithm (approximate raw segments, then local optimization).
  - **SeqComplex `SeqComplex.pm`** (rank 3) — `sub ce` Shannon form `K = −Σ pᵢ·log₂ pᵢ`;
    `sub cwf` Wootton–Federhen Stirling per-residue form `(Σlog_N(W) − Σlog_N(nᵦ))/tot`;
    `log_k(base,num)=log(num)/log(base)`.
  - **universalmotif `sequence_complexity`** (rank 3, Bioconductor) — WF score "reflects the
    numbers of each unique letter"; cites Wootton & Federhen 1993, Trifonov 1990, Orlov &
    Potapov 2004.
  - **Pei & Grishin 2005**, *Bioinformatics* 21(2):160 (rank 1) — SEG = "information measure of
    the complexity state vector … residue composition on a sliding window, no regard to pattern
    or periodicity"; two-pass, defaults W=12, K1=2.2, K2=2.5.
  - **Mier et al. / Shashidhara**, *Bioinformatics* 22(24):2980 (rank 1) — Shannon entropy
    complexity `−Σᵢ₌₁²⁰ pᵢ log pᵢ`, sum over the 20 amino acids.
- **Datasets (documented oracles):**
  - SEG defaults: W=12, K1=2.2, K2=2.5 bits/residue, max complexity log₂20 = 4.321928.
  - Worked window complexities (L=12, Shannon `K=−Σpᵢ·log₂pᵢ`, computed independently):
    12×A (homopolymer) → 0.000000; 11A/1B → 0.413817; 10A/2B → 0.650022; 8A/4B → 0.918296;
    6A/6B → 1.000000; 12 distinct → log₂12 = 3.584963.
- **Corner cases / failure modes:** sequence shorter than `W` has no complete trigger window →
  no low-complexity segment; extension `K2` must exceed `K1` to be effective; cutoffs must lie
  in [0, log₂N] (≤ 4.322 for 20 aa); a homopolymer window has complexity 0; entropy is defined
  only when total count > 1 (single-symbol window contributes 0).

## Deviations and assumptions

Two documented assumptions, both matching the DISORDER-LC-001 record:

1. **Complexity measure = Shannon entropy in bits/residue** — Wootton & Federhen equation (3)
   has two interconvertible operational forms: the exact compositional `K=(1/L)·log_N(L!/Πnᵢ!)`
   (NCBI `lnfact[]` permutation form) and the Shannon `−Σpᵢ·log₂pᵢ` form (NCBI `s_Entropy`,
   SeqComplex `ce`, Mier et al.). This implementation uses the Shannon bits/residue form —
   the one whose range (0…log₂20) and "bits" units are stated verbatim in the official man page.
2. **Empty / short-sequence behavior** — the spec defines complexity only for full windows of
   length `W`; sources do not prescribe a value for inputs shorter than `W`. The implementation
   returns no regions (empty result), consistent with "no complete trigger window exists."

Recommended coverage (MUST): homopolymer window → complexity 0, reported as a low-complexity
region; maximally diverse 12-distinct window → log₂12 ≈ 3.585 > K2, not low-complexity;
per-window complexity equals the exact Shannon value for biased windows; SEG defaults
W=12/K1=2.2/K2=2.5; poly-Q tract in a diverse protein → single region with inclusive
boundaries. SHOULD: sequence shorter than the window → no regions; two separated biased tracts
→ two regions. COULD: region complexity ∈ [0, log₂20]. No contradictions among sources — this
artifact fully agrees with the DISORDER-LC-001 evidence on the SEG method.
