---
type: source
title: "Evidence: ONCO-LOH-001 (genome-wide LOH detection + HRD-LOH score + per-chromosome LOH fraction)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-LOH-001-Evidence.md
sources:
  - docs/Evidence/ONCO-LOH-001-Evidence.md
source_commit: af81049fe33cf6b62fe3bf6864944552392fdeed
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-LOH-001

The validation-evidence artifact for test unit **ONCO-LOH-001** — **Loss of Heterozygosity (LOH)
detection and the HRD-LOH genomic-scar count** (`DetectLOH` + `CalculateLOHFraction`). The
**twenty-first ingested unit of the Oncology family** and one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is
synthesized in its own concept, [[loss-of-heterozygosity-detection]]; [[test-unit-registry]]
tracks the unit.

## What this file records

- **Online sources (mutually consistent, no contradictions):**
  - **Abkevich et al. (2012), Br J Cancer 107(10):1776–1782, PMID 23047548** (rank 1, primary) —
    the **HRD score = the number of intermediate-size LOH regions** in a tumour; intermediate-size
    LOH (not whole-chromosome LOH) correlates with defective BRCA1/BRCA2 (P = 10⁻¹¹).
  - **scarHRD reference R (Sztupinszki), `scarHRD.md`** (rank 3) — HRD-LOH = "the number of 15 Mb
    exceeding LOH regions which do not cover the whole chromosome"; the **15 Mb cut-off was
    arbitrarily selected** (middle of the studied interval).
  - **scarHRD `R/calc.hrd.R`** (rank 3, exact counting algorithm, verbatim R) — whole-chromosome
    exclusion (`chrDel`: a chromosome where **all** segments have minor CN `nB == 0` is a global-LOH
    chromosome, dropped); LOH segment = **minor CN 0 AND major CN ≠ 0** (`segSamp[,nB]==0 &
    segSamp[,nA]!=0` — allelic loss, not homozygous deletion); size filter **strictly `>`**
    (`segLOH[,4]-segLOH[,3] > sizelimit1`); score = count of remaining LOH segments; major CN `> 1`
    **capped to 1** before adjacent equal-state segments are merged (LOH state, not absolute CN,
    drives merging).
  - **scarHRD `R/scar_score.R`** (rank 3) — input column order (sample, chr, start, end, total CN,
    A-allele CN, B-allele CN); `sizelimitLOH = 15e6` default; within `calc.hrd` column 7 = A (major),
    column 8 = B (minor).
  - **oncoscanR `score_loh` (Christinat)** (rank 3, independent corroboration) — "All LOH segments
    larger than 15Mb but excluding chromosome with a global LOH alteration"; **merges overlapping or
    adjacent LOH segments (separated by 1bp)** before the size filter; "based on Abkevich et al.,
    Br J Cancer 2012".

- **Documented corner cases / failure modes:** homozygous deletion (minor = 0 **and** major = 0) is
  NOT LOH (fails `nA != 0`); whole-chromosome LOH (all segments minor = 0) excluded; strict size
  boundary (a segment of length **exactly** 15,000,000 bp is NOT counted — only `> 15 Mb`);
  **length = end − start** (columns 4 − 3, not `end − start + 1`); adjacent/overlapping same-state
  LOH segments merged first, so two adjacent < 15 Mb pieces forming a > 15 Mb region count as one.

- **Datasets (deterministic worked oracles):**
  - **HRD-LOH count on synthetic allele-specific segments → score = 1.** Seven segments straddling
    the 15 Mb boundary: chr1 0–20 Mb (major 1 / minor 0) is the **only** counted region; chr1
    20–60 Mb minor ≠ 1 → not LOH; chr2 0–10 Mb LOH but ≤ 15 Mb; chr3 0–16 Mb whole-chr LOH →
    excluded; chr4 0–30 Mb major 0 → homozygous deletion; chr5 0–15 Mb length exactly 15 Mb → not
    `> 15 Mb`; chr5 15–50 Mb minor ≠ 0.
  - **Per-chromosome LOH fraction (invariant 0 ≤ f ≤ 1):** chr1 (0–20 Mb minor 0, 20–60 Mb minor 1)
    → 20/60 = 0.333…; chr2 all het → 0.0; chr3 all LOH → 1.0.

- **Coverage recommendations:** MUST test HRD-LOH count on the synthetic dataset = 1; MUST test the
  exactly-15 Mb strict-`>` boundary (not counted); MUST test homozygous deletion (minor 0, major 0)
  not counted; MUST test whole-chromosome-LOH exclusion; MUST test a heterozygous-retained
  (minor ≠ 0) segment is not LOH; MUST test `CalculateLOHFraction` ∈ [0,1] on representative
  chromosomes; SHOULD test adjacent-LOH-merge → one region and null/empty → score 0, fraction 0;
  COULD test ordering invariance (shuffling input segments does not change the count).

## Deviations and assumptions

- **ASSUMPTION — LOH-fraction definition.** Abkevich 2012 / scarHRD define only the **count** of
  qualifying LOH regions (the HRD-LOH score), not a per-chromosome "LOH fraction". The Registry lists
  `CalculateLOHFraction(chromosome)` with invariant `0 ≤ fraction ≤ 1`. LOH fraction is defined as the
  **length-weighted fraction** of a chromosome's covered length lying under LOH segments (minor 0 &
  major ≠ 0) — the natural quantity satisfying the invariant and consistent with the segment-based LOH
  definition. Only this aggregation into a 0–1 fraction is a definitional/API choice; the
  correctness-affecting LOH-segment criterion itself is fully source-backed.
- **ASSUMPTION — input shape.** `DetectLOH` consumes allele-specific copy-number **segments** (chr,
  start, end, major CN, minor CN) — the scarHRD `seg`-table shape — not raw `tumorVcf/normalVcf`.
  Raw-VCF parsing and the upstream segmentation / BAF model are out of scope (handled upstream, cf.
  ONCO-CNA-001 / ONCO-ASCAT-001). API-shape decision; the counting logic is unchanged.

No source contradictions — the primary paper (Abkevich), the scarHRD reference implementation, and
the independent oncoscanR implementation each cover a disjoint part and agree on the LOH criterion,
the 15 Mb strict cut-off, and the whole-chromosome exclusion.
</content>
</invoke>
