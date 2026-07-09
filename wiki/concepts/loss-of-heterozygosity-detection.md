---
type: concept
title: "Genome-wide LOH detection (HRD-LOH count + per-chromosome LOH fraction)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-LOH-001-Evidence.md
source_commit: af81049fe33cf6b62fe3bf6864944552392fdeed
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-loh-001-evidence
      evidence: "Test Unit ID: ONCO-LOH-001 ... Algorithm: Loss of Heterozygosity (LOH) detection and HRD-LOH genomic-scar score"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:allele-specific-copy-number-ascat
      source: onco-loh-001-evidence
      evidence: "DetectLOH takes allele-specific copy-number segments (chromosome, start, end, major CN, minor CN), i.e. the scarHRD seg-table shape ... the upstream segmentation/BAF model are out of scope for this unit (handled upstream)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:homologous-recombination-deficiency-score
      source: onco-loh-001-evidence
      evidence: "The definition of a sample's HRD-LOH score is the number of 15 Mb exceeding LOH regions which do not cover the whole chromosome — the LOH component that DetectHRD sums with TAI and LST."
      confidence: high
      status: current
---

# Genome-wide LOH detection (HRD-LOH count + per-chromosome LOH fraction)

The Oncology family's standalone **loss-of-heterozygosity (LOH) caller** (**ONCO-LOH-001**):
`DetectLOH` counts the **HRD-LOH genomic-scar regions** (intermediate-size LOH segments) over a
set of allele-specific copy-number segments, and `CalculateLOHFraction` reports the length-weighted
LOH fraction of a single chromosome. The literature-traced record is [[onco-loh-001-evidence]];
[[test-unit-registry]] tracks the unit and [[algorithm-validation-evidence]] describes the
evidence-artifact pattern.

This is exactly the **HRD-LOH component** that the composite genomic-scar score
[[homologous-recombination-deficiency-score]] (ONCO-HRD-001) sums with TAI and LST — that unit
delegates its LOH count to this `DetectLOH`. Here the same count is the standalone deliverable, plus
the per-chromosome LOH-fraction API that HRD does not expose.

## 1. The LOH-segment criterion (Abkevich 2012 / scarHRD `calc.hrd.R`)

A segment is a counted **LOH region** when **all** of the following hold (reference `calc.hrd.R`,
verbatim R):

```
minor CN (nB) == 0   AND   major CN (nA) != 0      # allelic loss, one allele retained
length = end - start  >  15 Mb (15,000,000 bp)      # strictly greater; length = end − start
chromosome is not a whole-chromosome ("global") LOH chromosome
```

- **Minor 0 & major ≠ 0** — loss of *one* allele (the heterozygosity), the other retained. A segment
  with **minor 0 AND major 0** is a **homozygous deletion**, not LOH (both alleles lost), and fails
  the `nA != 0` clause.
- **Strict `> 15 Mb`** — the scarHRD comparison is `segLOH[,4]-segLOH[,3] > sizelimit1` with
  `sizelimitLOH = 15e6`. A segment of length **exactly** 15,000,000 bp is **not** counted. The 15 Mb
  cut-off was "arbitrarily selected … approximately in the middle of the interval" (scarHRD) — it
  captures *intermediate-size* LOH, which Abkevich 2012 showed correlates with BRCA1/2 deficiency
  (P = 10⁻¹¹), whereas whole-chromosome LOH does not.
- **Length = end − start** (columns 4 − 3), **not** `end − start + 1`.

**Whole-chromosome ("global LOH") exclusion (`chrDel`).** For each chromosome, if **every** segment
has minor CN 0, the whole chromosome is a global-LOH chromosome and **all** of its LOH segments are
dropped before counting (`if(all(segSamp[…,nB]==0)) chrDel <- c(chrDel, j)` then
`segLOH <- segLOH[!segLOH[,2] %in% chrDel,]`).

**HRD-LOH score** = `nrow(segLOH)` after those filters — the count of qualifying LOH regions.

## 2. Merging before counting

Major CN `> 1` is **capped to 1** before adjacent equal-state segments are merged
(`segSamp[segSamp[,nA]>1,nA] <- 1`), so the **LOH/non-LOH state** — not the absolute copy number —
drives which neighbours merge. Independently, oncoscanR `score_loh` **merges overlapping or adjacent
LOH segments (separated by ≤ 1 bp)** before the size filter: two adjacent < 15 Mb LOH pieces that
together exceed 15 Mb count as **one** region.

## 3. Per-chromosome LOH fraction (`CalculateLOHFraction`)

`CalculateLOHFraction(chromosome)` returns the **length-weighted fraction** of the chromosome's
covered length that lies under LOH segments (minor 0 & major ≠ 0), with invariant `0 ≤ fraction ≤ 1`
(LOH lengths are a subset of covered length). **This aggregation is a definitional / API choice** —
Abkevich 2012 and scarHRD define only the *count* of qualifying regions, not a fraction; the choice is
the natural quantity satisfying the stated invariant and is consistent with (not contradicted by) the
source LOH-segment criterion.

## Worked oracles

**HRD-LOH count → 1** (only the first segment qualifies):

| Chr | Start–End (bp) | Major | Minor | Length | Counted? | Why |
|-----|----------------|-------|-------|--------|----------|-----|
| 1 | 0–20 M | 1 | 0 | 20 M | **yes** | minor 0, major ≠ 0, > 15 Mb, chr1 not whole-LOH |
| 1 | 20 M–60 M | 1 | 1 | 40 M | no | minor ≠ 0 (het retained) |
| 2 | 0–10 M | 2 | 0 | 10 M | no | LOH but ≤ 15 Mb |
| 3 | 0–16 M | 1 | 0 | 16 M | no | whole chr3 all-LOH → excluded |
| 4 | 0–30 M | 0 | 0 | 30 M | no | homozygous deletion (major 0) |
| 5 | 0–15 M | 1 | 0 | 15 M | no | length exactly 15 Mb → not > 15 Mb |
| 5 | 15 M–50 M | 1 | 1 | 35 M | no | minor ≠ 0 |

**LOH fraction:** chr1 (0–20 M minor 0, 20–60 M minor 1) → 20/60 = **0.333…**; all-het chr → **0.0**;
all-LOH chr → **1.0**.

## Corner cases and assumptions

- **Homozygous deletion ≠ LOH:** minor 0 **and** major 0 is not counted (both alleles lost).
- **Strict 15 Mb boundary:** exactly 15 Mb → not counted (`>`, not `≥`).
- **Whole-chromosome LOH excluded** (Abkevich: whole-chromosome LOH does not correlate with HRD).
- **Empty / null input** → score 0, fraction 0. **Ordering invariance:** shuffling input segments
  does not change the count.
- **ASSUMPTION — input shape:** `DetectLOH` consumes allele-specific CN **segments** (chr, start,
  end, major CN, minor CN — the scarHRD `seg`-table shape), not raw tumour/normal VCFs; the upstream
  segmentation / BAF model is out of scope (handled by [[allele-specific-copy-number-ascat]] /
  ONCO-CNA-001).
- **ASSUMPTION — LOH-fraction definition** is the length-weighted 0–1 aggregation above (API choice,
  not source-mandated); the segment criterion is fully source-backed.

## Relation to the oncology family

`DetectLOH` reads its allele-specific major/minor CN segments off the upstream
[[allele-specific-copy-number-ascat]] layer — the same `AlleleSpecificSegment` substrate. Its HRD-LOH
count **is** the LOH term of the composite [[homologous-recombination-deficiency-score]]
(`HRD = LOH + TAI + LST`): that unit needs no centromere table for LOH precisely because this LOH
criterion is centromere-free (unlike its TAI/LST siblings). It is the **genome-wide** allele-specific
LOH caller, distinct from the **HLA-locus** specialization
[[hla-nomenclature-and-allele-specific-loh]] (LOHHLA), which calls LOH per HLA allele from a
continuous per-allele copy number plus an allelic-imbalance t-test rather than from integer segments.

## Scope and limitations

A [[scientific-rigor|research-grade]] correctness reference for genome-wide LOH detection and the
HRD-LOH count. The LOH-segment criterion, the 15 Mb strict cut-off, and the whole-chromosome
exclusion are Abkevich 2012 + scarHRD `calc.hrd.R` + oncoscanR `score_loh`; only the per-chromosome
LOH-fraction aggregation is a definitional/API choice. **Not for clinical or diagnostic use.** No
source contradictions — the primary paper and the two reference implementations each cover a disjoint
part and agree on the criterion, the cut-off, and the exclusion.
</content>
