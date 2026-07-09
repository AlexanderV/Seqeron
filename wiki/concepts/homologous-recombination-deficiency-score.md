---
type: concept
title: "HRD composite genomic-scar score (LOH + TAI + LST)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-HRD-001-Evidence.md
source_commit: ea6bdcb6f4ff447762681f38ff91dfacc4853d66
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-hrd-001-evidence
      evidence: "Test Unit ID: ONCO-HRD-001 ... Algorithm: Homologous Recombination Deficiency (HRD) composite genomic-scar score"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:allele-specific-copy-number-ascat
      source: onco-hrd-001-evidence
      evidence: "LOH, TAI and LST are all derived from AlleleSpecificSegment in DetectHRD(segments) ... reproducible from the repo's AlleleSpecificSegment (Major/Minor CN)"
      confidence: high
      status: current
---

# HRD composite genomic-scar score (LOH + TAI + LST)

The Oncology family's **homologous-recombination-deficiency (HRD) genomic-scar** unit
(**ONCO-HRD-001**): a tumour that has lost homologous-recombination repair accumulates
characteristic large-scale copy-number scars, and the HRD score sums three independent
scar counts into one instability index. The literature-traced record is
[[onco-hrd-001-evidence]]; [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the evidence-artifact pattern.

All three components are **derived per segment** from the allele-specific integer copy
numbers (major/minor CN) produced by the upstream copy-number layer
[[allele-specific-copy-number-ascat]] — `DetectHRD(segments)` consumes the same
`AlleleSpecificSegment` (Major/Minor CN) substrate. TAI and LST additionally need a
per-chromosome **centromere** coordinate table (embedded for GRCh38/GRCh37, see §4).

## 1. The composite score and cutoff (Telli 2016 / Stewart 2022)

The HRD score is the **unweighted sum** of the three genomic-scar counts:

```
HRD = LOH + TAI + LST          # no weighting (Telli et al. 2016, verbatim)
```

**HRD-high cutoff:** `HRD ≥ 42` (boundary **inclusive** — 42 is HRD-high, 41 is
HRD-negative). The 42 threshold is corroborated independently by Telli 2016 and the
Stewart 2022 review (myChoice CDx `gLOH + TAI + LST`). Each component is a **non-negative
integer event count**, so the sum is well-defined for any non-negative triple; a negative
input is invalid. The sum is commutative — component order does not matter.

## 2. Component definitions

- **HRD-LOH** (Abkevich 2012; oncoscanR `score_loh`) — count of **LOH regions longer than
  15 Mb but shorter than a whole chromosome**. A chromosome carrying a *global* (whole-chromosome)
  LOH alteration is **excluded** entirely. Reproduced by scarHRD `calc.hrd.R`: LOH segment =
  minor CN 0 & major CN ≠ 0, size `> 15 Mb`, drop any chromosome that is entirely LOH
  (`chrDel`). **Needs no centromere table.** This is the standalone LOH caller
  [[loss-of-heterozygosity-detection]] (ONCO-LOH-001, `DetectLOH`) reused here as the LOH term.
- **TAI** (Telomeric Allelic Imbalance; Birkbak 2012) — count of regions of allelic imbalance
  that **extend to a sub-telomere but do not cross the centromere**. Allelic imbalance =
  major CN ≠ minor CN (at least one allele present).
- **LST** (Large-scale State Transitions; Popova 2012) — count of **chromosomal breakpoints
  between adjacent ≥ 10 Mb regions** on the same arm, after smoothing out `< 3 Mb` small-scale
  CN variation.

The unweighted three-way sum is confirmed verbatim by scarHRD `scar_score.R`
(`sum_HRD0 <- res_lst + res_hrd + res_ai[1]`).

## 3. TAI derivation (scarHRD `calc.ai_new.R`, even/diploid path)

1. **Pre-filter:** drop segments shorter than `min.size` (default **1 Mb**).
2. **AI state:** `AI = 0` when major CN == minor CN (balanced); `AI = 2` (imbalanced) otherwise.
3. **Telomeric downgrade → TAI:** the **first** imbalanced segment whose **end < centromere
   start** (p-telomeric) and the **last** imbalanced segment whose **start > centromere end**
   (q-telomeric) are set to `AI = 1`. A first segment whose end **≥ centromere start** (crosses
   the centromere) is **not** telomeric.
4. **Whole-chromosome:** a single-segment imbalanced chromosome is whole-chromosome AI
   (`AI = 3`), **not** telomeric.
5. **TAI count** = number of `AI == 1` segments (`res_ai[1]`).

## 4. LST derivation (scarHRD `calc.lst.R`)

Autosomes only (sex chromosomes ignored). For each chromosome: split into p/q **arms** at the
centromere (`p.arm` start ≤ centromere start; `q.arm` end ≥ centromere end), clamp each arm to
the centromere, skip chromosomes with `< 2` segments. **Iteratively** remove arm segments
`< 3 Mb` (re-merging after each removal) until none remain. Flag each surviving segment `1` if
length `≥ 10 Mb` else `0`; count **one LST** per adjacent pair where **both flags are 1 AND the
gap between them is `< 3 Mb`**. Sample LST = sum of counted breaks over all arms.

**Centromere coordinate table.** The arm split depends on the exact centromere `[start, end]`
per chromosome. These are the UCSC cytoBand **`acen`** boundaries (p11 acen start → q11 acen
end), embedded in the Evidence as a published **GRCh38 + GRCh37** per-chromosome table and
cross-verified against the NCBI GRC modeled-centromere table (e.g. GRCh38 CEN1 modeled
122.0–125.2 Mb agrees to cytoband resolution with the acen 121.7–125.1 Mb). Retrieving these
public coordinates is what resolved the prior "centromere table unretrievable" blocker, making
TAI and LST derivable rather than un-computable.

## Worked oracles

| LOH | TAI | LST | Score | Status |
|-----|-----|-----|-------|--------|
| 14 | 14 | 14 | 42 | HRD-high (boundary) |
| 14 | 13 | 14 | 41 | HRD-negative |
| 20 | 15 | 12 | 47 | HRD-high |
| 5 | 4 | 3 | 12 | HRD-negative |
| 0 | 0 | 0 | 0 | HRD-negative (near-diploid) |

Component oracles: TAI — first-segment imbalance ending before centromere start (p-telomeric) and
last-segment imbalance starting after centromere end (q-telomeric) counted; interstitial and
single-segment whole-chromosome imbalance not counted; sub-1 Mb dropped. LST — two adjacent
≥ 10 Mb segments separated by < 3 Mb → 1 break; a < 10 Mb neighbour or a ≥ 3 Mb gap → 0;
iterative 3 Mb smoothing can expose a break; `< 2` segments skipped.

## Corner cases and assumptions

- **Boundary at 42:** inclusive (`≥ 42` HRD-high); 42 → high, 41 → negative.
- **Near-diploid / low signal:** small counts → low sum → HRD-negative (e.g. `0+0+0=0`).
- **Non-negative counts:** each component is an event count; a negative component / negative
  score is rejected (`ArgumentOutOfRangeException`).
- **ASSUMPTION — even-ploidy AI path for TAI:** `AlleleSpecificSegment` carries major/minor CN
  but **not** scarHRD's ASCAT per-sample ploidy / aberrant-cell-fraction columns, so the
  implemented AI rule is scarHRD's default **even/diploid path** (AI present ⟺ major ≠ minor,
  the literal `seg[,7]==seg[,8]` test). The odd-ploidy ploidy-normalised branch is not
  reproduced (it needs the absent ploidy column). This is the dominant path and matches
  Birkbak's "regions of allelic imbalance".

## Relation to the copy-number family

HRD is a **downstream summary** over the allele-specific segments of
[[allele-specific-copy-number-ascat]]: LOH, TAI and LST are all read off the same major/minor
CN segments, so this unit is the genomic-scar aggregation layer sitting above ASCAT. It is
distinct from the whole-chromosome / arm-scale [[aneuploidy-detection]] (a total-CN depth-ratio
caller with no allelic contrast) — HRD counts *sub-chromosomal* allele-specific scars. The
centromere coordinates it embeds are the same biological landmarks that the sequence-level
[[centromere-analysis]] unit detects de novo, but here they enter as a published coordinate
lookup rather than a detection target.

## Scope and limitations

A [[scientific-rigor|research-grade]] correctness reference for the HRD genomic-scar score and
its LOH/TAI/LST derivation. The composite sum + 42 cutoff are Telli 2016 / Stewart 2022; the
component definitions are Abkevich 2012 (LOH), Birkbak 2012 (TAI), Popova 2012 (LST), with the
scarHRD reference R (`calc.hrd.R` / `calc.ai_new.R` / `calc.lst.R` / `scar_score.R`) supplying
the operational algorithm, and UCSC cytoBand + NCBI GRC supplying the centromere coordinates.
**Not for clinical or diagnostic use.** No source contradictions — the primary papers, the
reference implementation, and the coordinate databases each cover a disjoint part and agree on
the sum, the cutoff, and the component rules.
</content>
</invoke>
