---
type: source
title: "Evidence: POP-SELECT-001 (Selection signature detection — iHS / EHH scan)"
tags: [validation, population-genetics]
doc_path: docs/Evidence/POP-SELECT-001-Evidence.md
sources:
  - docs/Evidence/POP-SELECT-001-Evidence.md
source_commit: 0d5c33fdcb6b264e682c9d593a37867121e26d99
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: POP-SELECT-001

The validation-evidence artifact for test unit **POP-SELECT-001** — **recent-positive-selection
signature detection** via the **integrated Haplotype Score (iHS)** and its underlying **Extended
Haplotype Homozygosity (EHH)** (`CalculateEhh`, `CalculateIHS`, `StandardizeIHS`,
`ScanForSelection`). It is one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the formulae, worked oracles, invariants,
sign convention, and corner cases are synthesized in the dedicated concept
[[selection-scan-ihs-ehh]]. This is a population-genetics `POP-*` unit that **consumes** derived/
ancestral allele frequencies from [[allele-genotype-frequencies]] (MAF filter, frequency-bin
standardization) and quantifies the same haplotype-structure phenomenon measured by
[[linkage-disequilibrium]]. It sits in the family anchored by [[ancestry-estimation-admixture]]; see
[[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources (authority-ranked, all cross-verified):**
  - **Voight, Kudaravalli, Wen & Pritchard (2006)**, *PLoS Biology* 4(3):e72 — the primary iHS
    reference. Defines iHH as the trapezoidal area under the EHH curve summed both directions from the
    core SNP (`iHH_A` ancestral, `iHH_D` derived); **unstandardized iHS = ln(iHH_A / iHH_D)**
    (ancestral numerator); standardization to an approx. standard normal within
    derived-allele-frequency bins; integration truncated where EHH drops below **0.05**; genome scan
    quantified by the **proportion of SNPs with |iHS| > 2** in ~50-SNP windows; **MAF > 5%** filter.
  - **Sabeti et al. (2002)**, *Nature* 419:832–837 — originator of EHH: the probability that two
    randomly chosen core-allele chromosomes are homozygous over a surrounding region.
  - **Szpiech & Hernandez (2014)** *selscan*, *MBE* 31(10):2824–2827 — reference implementation.
    EHH Eq. 3 `EHH_c = Σ_{h} C(n_h,2) / C(n_c,2)`; iHH by trapezoidal quadrature weighted by genetic
    distance; **explicit note that selscan uses `ln(iHH_1/iHH_0)` (derived/ancestral) — the opposite
    sign to Voight** (ancestral/derived). Confirms the two conventions differ by sign.
  - **rehh vignette** (Gautier, Klassmann & Vitalis), CRAN — EHH formula algebraically identical
    (`1/(n_a(n_a−1)) · Σ n_k(n_k−1)`); default cutoff `limehh = 0.05`; default `freqbin = 0.05`.
    Worked `calc_ehh()` output for SNP F1205400: `IHH_A = 284429.9`, `IHH_D = 2057107.4` ⇒ Voight
    unstandardized iHS = ln(284429.9/2057107.4) = **−1.978569274** (verified this session).

- **Datasets / oracles:**
  - **rehh worked ratio:** `iHH_A = 284429.9`, `iHH_D = 2057107.4` ⇒ ln(A/D) = −1.978569274.
  - **Constructed haplotype panel** (core index 2, positions 0,10,20,30,40): 3× identical derived
    `AA1GG` flanks + three all-distinct ancestral haplotypes ⇒ `EHH_D = 1.0` at each flank,
    `EHH_A = 0.0` (truncated at first flank), `iHH_D = 40.0`, `iHH_A = 10.0`, unstandardized
    iHS = ln(10/40) = ln(0.25) = **−1.386294361**, derived freq 0.5 (long derived haplotype ⇒ negative
    iHS under Voight).
  - **EHH unit values (selscan Eq. 3):** `11,11,11,10` (n_c=4) ⇒ (C(3,2)+C(1,2))/C(4,2) = 3/6 = 0.5;
    `00,00,01,01` ⇒ 2/6 = 0.3333; single haplotype ⇒ 1.0; three all-distinct ⇒ 0.0.

- **Documented corner cases / failure modes:** monomorphic core or no ancestral state ⇒ no iHS
  reported (throws for monomorphic core); balanced EHH decay ⇒ `iHH_A/iHH_D ≈ 1` ⇒ unstandardized
  iHS ≈ 0; integration truncates at the first marker where EHH < 0.05; single-element frequency bin ⇒
  standardized iHS = 0; edge cases throw on null inputs, empty/inconsistent-length haplotypes, invalid
  core allele, out-of-range `coreIndex`. Long-gap / chromosome-end masking (>200 kb) is a
  data-curation detail, not part of the in-memory core score.

## Deviations and assumptions

No algorithm **deviation** — the core `EHH`, `iHH` (trapezoid + 0.05 cutoff), and unstandardized
`iHS = ln(iHH_A/iHH_D)` are exact matches to Voight (2006) / selscan / rehh, and the implementation
adopts **Voight's ancestral/derived sign convention** (the opposite of selscan). Two documented
**assumptions** affect only the standardization magnitude, not the sign, ordering, or the canonical
unstandardized iHS: (1) the standardization SD uses the **sample estimator (N−1)** — Voight does not
state N vs N−1; (2) frequency binning defaults to **20 equal-width bins (0.05)** matching rehh
`freqbin`, but the bin count is an explicit, caller-overridable parameter. No source contradictions
(the Voight/selscan sign difference is a documented convention, not a conflict). Open Questions: none.
