# Homologous Recombination Deficiency (HRD) Score

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-HRD-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-23 |

## 1. Overview

The HRD score is a composite genomic-instability biomarker used to identify tumours with homologous-recombination repair deficiency (e.g. *BRCA1/2*, *PALB2*), which predicts response to platinum chemotherapy and PARP inhibitors. It is the unweighted sum of three independent genomic-scar metrics: loss of heterozygosity (LOH), telomeric allelic imbalance (TAI), and large-scale state transitions (LST) [1]. A tumour is classified HRD-high when the sum is at or above a fixed clinical cutoff of 42 [1]. This implementation **derives all three components end-to-end from allele-specific copy-number segments** — LOH (Abkevich 2012 / scarHRD `calc.hrd` [2][7]), TAI (Birkbak 2012 / scarHRD `calc.ai_new` [3][7]), and LST (Popova 2012 / scarHRD `calc.lst` [4][7]) — using an embedded per-build centromere coordinate table (UCSC cytoBand `acen` [8], cross-verified vs NCBI GRC [9]), then computes the composite sum and the threshold classification. A caller-supplied-TAI/LST overload remains for externally computed components. The status is "Simplified" because TAI uses scarHRD's default even-ploidy allelic-imbalance rule and does not reproduce the odd-ploidy ploidy-renormalization branch (see §5.3).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Homologous-recombination-deficient tumours accumulate characteristic genome-wide "scars". Three independently developed metrics quantify these scars: HRD-LOH (Abkevich 2012) [2], NtAI/TAI (Birkbak 2012) [3], and LST (Popova 2012) [4]. Telli et al. (2016) combined them into a single score and established the clinical HRD-high cutoff used by the Myriad myChoice CDx assay [1][5].

### 2.2 Core Model

The combined HRD score is the **unweighted sum** of the three component counts [1]:

```
HRD = LOH + TAI + LST
```

where each component is a non-negative integer count of genomic events:

- **LOH** — number of LOH regions longer than 15 Mb but shorter than a whole chromosome [2][6].
- **TAI** — number of allelic-imbalance regions that extend to a sub-telomere but do not cross the centromere; allelic imbalance is when "the copy number of the two alleles were not equal, and at least one allele was present" [3]. Operationally (scarHRD `calc.ai_new` [7]): per chromosome, after dropping < 1 Mb segments and merging same-allele-state runs, the first segment counts when it is imbalanced and its end is before the centromere start (p-telomeric), and the last segment counts when it is imbalanced and its start is after the centromere end (q-telomeric); a single imbalanced segment is whole-chromosome AI, not telomeric.
- **LST** — number of chromosomal breaks between two adjacent regions each ≥ 10 Mb, after iteratively smoothing away regions < 3 Mb, with the gap between the pair < 3 Mb, counted per chromosome arm (Popova 2012 [4]; scarHRD `calc.lst` [7]).

Classification [1]:

```
status = HRD-high   if HRD ≥ 42
         HRD-negative otherwise
```

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | HRD = LOH + TAI + LST | unweighted sum [1] |
| INV-02 | Sum is order-independent (commutative) | integer addition; "unweighted" [1] |
| INV-03 | status = HRD-high iff score ≥ 42 (inclusive) | Telli 2016 cutoff "≥42" [1] |
| INV-04 | score ≥ 0 and each component ≥ 0 | components are event counts [2][3][4] |
| INV-05 | `DetectHRD(segments, tai, lst).Components.Loh == DetectLOH(segments).Score` | LOH derived from segments, not supplied [2][7] |
| INV-06 | TAI and LST are ≥ 0; sex-chromosome segments contribute 0 to both | autosome-only centromere table; scarHRD excludes X/Y from LST [7] |
| INV-07 | `DetectHRD(segments).Components == (DetectLOH(segments).Score, CalculateHrdTaiScore(segments), CalculateHrdLstScore(segments))` | all-derived overload sums the three standalone derivations [7] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| loh | int | required | HRD-LOH component count | ≥ 0 |
| tai | int | required | TAI/NtAI component count | ≥ 0 |
| lst | int | required | LST component count | ≥ 0 |
| score | int | required | combined HRD score for classification | ≥ 0 |
| components | HrdComponents | required | bundle of (Loh, Tai, Lst) for `DetectHRD` | each ≥ 0 |
| segments | IEnumerable&lt;AlleleSpecificSegment&gt; | required | allele-specific CN segments the components are derived from | non-null; End > Start; CN ≥ 0 |
| genome | ReferenceGenome | GRCh38 | reference assembly supplying the centromere table for TAI/LST | GRCh38 or GRCh37 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| CalculateHRDScore | int | combined score = loh + tai + lst |
| ClassifyHRDStatus | HrdStatus | HrdHigh (≥ 42) or HrdNegative |
| DetectHRD | HrdResult | components, summed score, and status |
| DetectHRD(segments,tai,lst) | HrdResult | LOH derived from segments + caller-supplied TAI/LST, summed score, and status |
| DetectHRD(segments,genome) | HrdResult | all three components (LOH+TAI+LST) derived from segments, summed score, and status |
| CalculateHrdTaiScore(segments,genome) | int | derived HRD-TAI count (≥ 0) |
| CalculateHrdLstScore(segments,genome) | int | derived HRD-LST count (≥ 0) |

### 3.3 Preconditions and Validation

Any negative component count throws `ArgumentOutOfRangeException` (`CalculateHRDScore`, `DetectHRD`). A negative score throws `ArgumentOutOfRangeException` (`ClassifyHRDStatus`). The cutoff comparison is inclusive at 42. The segment-driven `DetectHRD(segments, tai, lst)` throws `ArgumentNullException` for null segments, `ArgumentException` for a segment with non-positive length or negative copy number (via `DetectLOH`), and `ArgumentOutOfRangeException` for a negative `tai`/`lst`.

## 4. Algorithm

### 4.1 High-Level Steps

For the supplied-counts path:

1. Validate that each component count (and the score, for classification) is non-negative.
2. Sum the three components: `score = loh + tai + lst`.
3. Classify: `HrdHigh` if `score ≥ 42`, else `HrdNegative`.

For the all-derived path `DetectHRD(segments, genome)`:

1. Derive HRD-LOH via `DetectLOH` (Abkevich/scarHRD `calc.hrd`).
2. Derive HRD-TAI via `CalculateHrdTaiScore`: per autosome, drop < 1 Mb segments, merge same-allele-state runs; count the first segment if imbalanced with end < centromere start (p-telomeric) and the last if imbalanced with start > centromere end (q-telomeric); a single imbalanced segment is whole-chromosome AI (not counted).
3. Derive HRD-LST via `CalculateHrdLstScore`: per autosome, split into p-arm (start ≤ centromere start) and q-arm (end ≥ centromere end), clamp arm edges to the centromere, merge same-state, iteratively remove < 3 Mb segments (re-merging), then count adjacent pairs both ≥ 10 Mb with a < 3 Mb gap.
4. Sum and classify as above.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

| Constant | Value | Source |
|----------|-------|--------|
| HRD-high cutoff | 42 (inclusive) | Telli 2016 [1]; Stewart 2022 [5] |
| TAI minimum segment length | 1 Mb | scarHRD `calc.ai_new` `min.size=1e6` [7] |
| LST smoothing / adjacency gap | 3 Mb | Popova 2012 [4]; scarHRD `calc.lst` `3e6` [7] |
| LST large-segment threshold | 10 Mb | Popova 2012 [4]; scarHRD `calc.lst` `10e6` [7] |

**Centromere coordinate table** (p-arm end = centromere start; q-arm start = centromere end), embedded for chromosomes 1–22 from the UCSC cytoBand `acen` bands of hg38 and hg19 [8], cross-verified against the NCBI GRC modeled-centromere table [9]. Example (GRCh38, bp): chr1 121,700,000–125,100,000; chr2 91,800,000–96,000,000; chr17 22,700,000–27,400,000; chr21 10,900,000–13,000,000. The table is autosome-only — sex chromosomes are excluded from TAI and LST (scarHRD excludes chr23/24/X/Y [7]).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateHRDScore / ClassifyHRDStatus / DetectHRD(components) | O(1) | O(1) | three integer additions + one comparison |
| CalculateHrdTaiScore | O(n log n) | O(n) | per-chromosome sort/merge over n segments |
| CalculateHrdLstScore | O(n²) worst case | O(n) | iterative < 3 Mb smoothing with re-merge per removal |
| DetectHRD(segments) | O(n²) worst case | O(n) | dominated by LST smoothing |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.CalculateHRDScore(int, int, int)`: returns the unweighted sum.
- `OncologyAnalyzer.ClassifyHRDStatus(int)`: returns `HrdStatus` against the 42 cutoff.
- `OncologyAnalyzer.DetectHRD(HrdComponents)`: end-to-end sum + classification from supplied counts.
- `OncologyAnalyzer.DetectHRD(IEnumerable<AlleleSpecificSegment>, int tai, int lst)`: derives the HRD-LOH count from segments via `DetectLOH` (Abkevich 2012 / scarHRD `calc.hrd`), then sums with caller-supplied TAI/LST and classifies.
- `OncologyAnalyzer.DetectHRD(IEnumerable<AlleleSpecificSegment>, ReferenceGenome)`: derives all three components (LOH, TAI, LST) from segments, sums (`sum_HRD0`), and classifies.
- `OncologyAnalyzer.CalculateHrdTaiScore(IEnumerable<AlleleSpecificSegment>, ReferenceGenome)`: derives the HRD-TAI count (scarHRD `calc.ai_new`, even-ploidy path).
- `OncologyAnalyzer.CalculateHrdLstScore(IEnumerable<AlleleSpecificSegment>, ReferenceGenome)`: derives the HRD-LST count (scarHRD `calc.lst`, `chr.arm='no'`).

### 5.2 Current Behavior

`CalculateHRDScore` / `ClassifyHRDStatus` / `DetectHRD(HrdComponents)` take the three counts as already-computed inputs. `DetectHRD(segments, genome)` derives all three components end-to-end from allele-specific segments using the embedded per-build centromere table; `DetectHRD(segments, tai, lst)` derives only LOH and accepts caller-supplied TAI/LST (for externally computed components). The cutoff is exposed as the public constant `HrdHighScoreThreshold = 42`. TAI/LST derivation reuses the per-chromosome grouping + a same-allele-state merge (equal major AND minor) and the `ReferenceGenome` enum / `TryGetAutosomeNumber` helper shared with the WGD/ploidy code. No substring/pattern search is involved, so the repository suffix tree is not applicable (N/A).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- HRD = LOH + TAI + LST as an unweighted sum [1] (scarHRD `scar_score.R`: `sum_HRD0 <- res_lst + res_hrd + res_ai[1]` [7]).
- HRD-high classification at score ≥ 42, boundary inclusive [1][5].
- **HRD-LOH derived from allele-specific segments** via `DetectLOH` — LOH regions with minor==0 & major≠0, length > 15 Mb, whole-chromosome-LOH excluded, after same-state merge (Abkevich 2012 [2]; scarHRD `calc.hrd.R` [7]; oncoscanR `score_loh` [6]).
- **HRD-TAI derived from segments** (`CalculateHrdTaiScore`): the telomeric-downgrade rule of scarHRD `calc.ai_new` [7] — < 1 Mb filter, same-state merge, first-segment-end < centromere start → p-telomeric, last-segment-start > centromere end → q-telomeric, single imbalanced segment → whole-chr (not counted) [3].
- **HRD-LST derived from segments** (`CalculateHrdLstScore`): the per-arm rule of scarHRD `calc.lst` [7] — p/q-arm split at the centromere with edge clamping, iterative < 3 Mb smoothing, adjacent ≥ 10 Mb pair with < 3 Mb gap counted as a break [4].
- **Centromere coordinate table** embedded from UCSC cytoBand `acen` (hg38 + hg19) [8], cross-verified vs NCBI GRC modeled centromeres [9], replacing scarHRD's binary `chrominfo` object with citable public coordinates.

**Intentionally simplified:**

- **TAI even-ploidy path only.** Allelic imbalance is assigned by the literal even/diploid rule major ≠ minor (scarHRD `calc.ai_new` `seg[,7]==seg[,8]` path [7]). **Consequence:** on odd-ploidy chromosomes scarHRD applies an additional ploidy-renormalized AI rule using its ASCAT per-sample ploidy and aberrant-cell-fraction columns, which `AlleleSpecificSegment` does not carry; that branch is not reproduced, so TAI may differ from scarHRD on odd-ploidy chromosomes.

**Not implemented:**

- **scarHRD's ASCAT ploidy/cellularity preprocessing** (per-sample ploidy and aberrant-cell-fraction columns); **users should rely on:** supplying already-segmented allele-specific copy number (major/minor CN) — the contamination filter (`cont`) and ploidy-by-chromosome renormalization are out of scope because the inputs they need are not part of `AlleleSpecificSegment`.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Score exactly 42 | HRD-high | inclusive cutoff [1] |
| Score 41 | HRD-negative | below cutoff [1] |
| All components 0 (near-diploid) | score 0, HRD-negative | low-signal tumour; sum well-defined |
| Negative component or score | ArgumentOutOfRangeException | counts are non-negative |

### 6.2 Limitations

All three components (LOH, TAI, LST) are derived from segments. TAI uses scarHRD's even-ploidy allelic-imbalance rule and does not reproduce the odd-ploidy ploidy-renormalization branch (which needs ASCAT ploidy/cellularity columns absent from `AlleleSpecificSegment`; see §5.3). TAI and LST are computed on autosomes only (the centromere table is autosome-only, matching scarHRD's exclusion of chr23/24/X/Y). Centromere coordinates are the UCSC cytoBand `acen` band boundaries (band-resolution, the same lineage scarHRD's `chrominfo` derives from); the all-derived path is therefore a faithful reproduction of scarHRD's diploid/even-ploidy behavior rather than a bit-for-bit match of a specific scarHRD release's internal coordinates. The 42 cutoff is the myChoice/Telli-2016 value validated primarily in breast and ovarian carcinoma [1][5]; tissue-specific cutoffs are out of scope.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// Supplied-components path:
var result = OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(Loh: 20, Tai: 15, Lst: 12));
// result.Score  == 47
// result.Status == OncologyAnalyzer.HrdStatus.HrdHigh   (47 ≥ 42)

// All-derived path: LOH, TAI and LST all derived from allele-specific segments:
var allDerived = OncologyAnalyzer.DetectHRD(segments); // default GRCh38
// allDerived.Components.Loh / .Tai / .Lst are each derived from the segments
// (scarHRD calc.hrd / calc.ai_new / calc.lst); use ReferenceGenome.GRCh37 for hg19.

// Caller-supplied TAI/LST path (externally computed components), LOH still derived:
var fromSegments = OncologyAnalyzer.DetectHRD(segments, tai: 25, lst: 16);
```

**Performance baseline:** TAI/LST are O(n log n)/O(n²)-worst-case over the (small) per-chromosome segment count; segmentation tables have at most a few thousand segments genome-wide, so the worst-case smoothing loop is negligible. The 42-test fixture (which exercises all derivation paths on multi-segment chromosomes) runs in ~11 ms total, dominated by NUnit setup; no separate BenchmarkDotNet baseline is warranted at this input scale.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_CalculateHRDScore_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_CalculateHRDScore_Tests.cs) — covers `INV-01`–`INV-07`
- Evidence: [ONCO-HRD-001-Evidence.md](../../../docs/Evidence/ONCO-HRD-001-Evidence.md)
- Related algorithm: [Loss of Heterozygosity (HRD-LOH)](Loss_Of_Heterozygosity.md) (ONCO-LOH-001) — the derived LOH component

## 8. References

1. Telli ML, Timms KM, Reid J, et al. 2016. Homologous Recombination Deficiency (HRD) Score Predicts Response to Platinum-Containing Neoadjuvant Chemotherapy in Patients with Triple-Negative Breast Cancer. Clin Cancer Res 22(15):3764–3773. https://pubmed.ncbi.nlm.nih.gov/26957554/
2. Abkevich V, Timms KM, Hennessy BT, et al. 2012. Patterns of genomic loss of heterozygosity predict homologous recombination repair defects in epithelial ovarian cancer. Br J Cancer 107(10):1776–1782. https://www.nature.com/articles/bjc2012451
3. Birkbak NJ, Wang ZC, Kim JY, et al. 2012. Telomeric allelic imbalance indicates defective DNA repair and sensitivity to DNA-damaging agents. Cancer Discov 2(4):366–375. https://pmc.ncbi.nlm.nih.gov/articles/PMC3806629/
4. Popova T, Manié E, Rieunier G, et al. 2012. Ploidy and large-scale genomic instability consistently identify basal-like breast carcinomas with BRCA1/2 inactivation. Cancer Res 72(21):5454–5462. https://aacrjournals.org/cancerres/article/72/21/5454/576090/
5. Stewart MD, Merino Vega D, Arend RC, et al. 2022. Homologous Recombination Deficiency: Concepts, Definitions, and Assays. Oncologist 27(3):167–174. https://pmc.ncbi.nlm.nih.gov/articles/PMC8914493/
6. Christinat Y. oncoscanR `score_loh` documentation. https://rdrr.io/github/yannchristinat/oncoscanR/man/score_loh.html
7. Sztupinszki Z, Diossy M, Krzystanek M, et al. 2018. scarHRD — reference implementation (`R/calc.hrd.R`, `R/calc.ai_new.R`, `R/calc.lst.R`, `R/scar_score.R`). GitHub. https://github.com/sztup/scarHRD
8. UCSC Genome Browser. cytoBand track (`acen` bands = centromere), hg38 and hg19. https://api.genome.ucsc.edu/getData/track?genome=hg38;track=cytoBand
9. NCBI Genome Reference Consortium. GRCh38 modeled centromeres. https://www.ncbi.nlm.nih.gov/grc/human
