---
type: source
title: "Evidence: ONCO-HRD-001 (HRD composite genomic-scar score — LOH + TAI + LST)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-HRD-001-Evidence.md
sources:
  - docs/Evidence/ONCO-HRD-001-Evidence.md
source_commit: ea6bdcb6f4ff447762681f38ff91dfacc4853d66
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-HRD-001

The validation-evidence artifact for test unit **ONCO-HRD-001** — the **Homologous Recombination
Deficiency (HRD) composite genomic-scar score** (`HRD = LOH + TAI + LST`). The **nineteenth ingested
unit of the Oncology family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is synthesized in its
own concept, [[homologous-recombination-deficiency-score]]; [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (mutually consistent, no contradictions):**
  - **Telli et al. (2016), Clin Cancer Res 22(15):3764–3773** (rank 1) — composite HRD score is an
    **unweighted sum of LOH, TAI, and LST**; HRD-high defined as **score ≥ 42** (or BRCA1/2 mutation).
  - **Stewart et al. (2022) review, Oncologist 27(3):167–174, PMC8914493** (rank 1) — component
    definitions: gLOH/HRD-LOH = LOH regions **> 15 Mb and < whole chromosome** (Abkevich); TAI =
    regions of allelic imbalance extending to a sub-telomere but not crossing the centromere (Birkbak);
    LST = chromosome breaks (Popova); myChoice CDx `gLOH + TAI + LST`; cutoff **42** (independent
    corroboration of Telli).
  - **Birkbak et al. (2012), Cancer Discov, PMC3806629** (rank 1) — allelic imbalance = the two alleles'
    CN not equal with ≥ 1 allele present; **telomeric** AI/CNA = extends to a sub-telomere but does not
    cross the centromere.
  - **Popova et al. (2012), Cancer Res 72(21):5454** (rank 1) — LST = a chromosomal breakpoint between
    **adjacent ≥ 10 Mb regions** after smoothing/filtering `< 3 Mb` small-scale CN variation.
  - **oncoscanR `score_loh`** (rank 3) — LOH segments `> 15 Mb`, excluding any chromosome with a global
    LOH alteration (confirms 15 Mb minimum + whole-chromosome exclusion, Abkevich 2012).
  - **scarHRD reference R (Sztupinszki 2018)** (rank 3) — `calc.hrd.R` (LOH: minor CN 0 & major CN ≠ 0,
    size `> 15 Mb`, drop whole-chromosome-LOH; **no centromere table needed** — exactly what `DetectLOH`
    does); `calc.ai_new.R` (TAI: pre-filter `< min.size` default 1 Mb; AI = 2 iff major ≠ minor;
    first-segment-end < centromere start → p-telomeric `AI=1`, last-segment-start > centromere end →
    q-telomeric `AI=1`; single-segment imbalance → whole-chr `AI=3`; TAI = count of `AI==1`);
    `calc.lst.R` (LST: autosomes only, p/q arm split at centromere, iterative 3 Mb shrink, adjacent
    pair both ≥ 10 Mb AND gap < 3 Mb → 1 break); `scar_score.R` (`sum_HRD0 <- res_lst + res_hrd +
    res_ai[1]` — the unweighted three-way sum).
  - **UCSC cytoBand `acen`** (rank 5) + **NCBI GRC modeled centromeres** (rank 2) — per-chromosome
    centromere `[start, end]` for **GRCh38 and GRCh37**, embedded in the Evidence as a published table
    and cross-verified between the two databases (agree to cytoband resolution).

- **Documented corner cases / failure modes:** boundary at 42 inclusive (42 → HRD-high, 41 → negative);
  near-diploid low-signal tumours → low sum → HRD-negative (`0+0+0=0`); components are non-negative
  integer event counts (negative input invalid). TAI: interstitial and single-segment whole-chromosome
  imbalance not counted, sub-1 Mb segments dropped, first segment crossing the centromere not telomeric.
  LST: `< 2` segments skipped, sex chromosomes ignored, iterative 3 Mb smoothing can expose a break.

- **Datasets (deterministic worked oracles):**
  - **Composite (from the sum + cutoff):** (14,14,14) → 42 HRD-high (boundary) · (14,13,14) → 41
    HRD-negative · (20,15,12) → 47 HRD-high · (5,4,3) → 12 HRD-negative · (0,0,0) → 0 HRD-negative.
  - **Centromere table:** embedded GRCh38 / GRCh37 per-chromosome `acen` `[start,end]` matches UCSC
    cytoBand values.

- **Coverage recommendations:** MUST test composite `LOH+TAI+LST` for representative triples; MUST test
  classification at the 42 boundary (42 → high, 41 → negative); MUST test end-to-end `DetectHRD`; MUST
  test the TAI telomeric-vs-interstitial + sub-1 Mb rules and the centromere-crossing negative; MUST
  test the LST adjacent-≥10 Mb / <3 Mb-gap break rule + 3 Mb smoothing + skip rules; MUST verify the
  embedded centromere table vs UCSC; SHOULD test near-diploid all-zero → 0 and negative-component
  rejection; COULD test commutativity of the sum.

## Deviations and assumptions

- **RESOLVED — TAI/LST now derived from segments.** The prior blocker (centromere `chrominfo` table
  unretrievable) is resolved: the per-chromosome centromere boundaries are the UCSC cytoBand `acen`
  regions, cross-verified against the NCBI GRC modeled-centromere table and embedded as a published
  GRCh38/GRCh37 table. LOH, TAI and LST are all derived from `AlleleSpecificSegment` in
  `DetectHRD(segments)`.
- **ASSUMPTION — even-ploidy / standard allelic-imbalance path for TAI.** `AlleleSpecificSegment`
  carries major/minor allele CN but not scarHRD's ASCAT per-sample ploidy / aberrant-cell-fraction
  columns, so the implemented AI rule is scarHRD's default even/diploid path (AI present ⟺ major ≠
  minor, the literal `seg[,7]==seg[,8]` test). The odd-ploidy ploidy-normalised branch is not
  reproduced (it needs the absent ploidy column); this is the dominant path and matches Birkbak's
  "regions of allelic imbalance" (documented as an intentional simplification in §5.3 of the algorithm
  doc).

No source contradictions — the primary papers (Telli/Stewart/Abkevich/Birkbak/Popova), the scarHRD
reference implementation, and the UCSC/NCBI centromere coordinate databases each cover a disjoint part
and agree on the unweighted sum, the 42 cutoff, and the LOH/TAI/LST component rules.
