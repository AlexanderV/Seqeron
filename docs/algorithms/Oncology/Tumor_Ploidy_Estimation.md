# Tumor Ploidy Estimation

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-PLOIDY-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Estimates the average ploidy of a tumour genome from allele-specific copy-number segments and classifies whether the genome has undergone whole-genome doubling (WGD). Average ploidy ψ is the segment-length-weighted mean of per-segment total copy number [1]; on the n-scale a pure-diploid genome has ψ = 2.0 and elevated values (e.g. >2.7n) mark aneuploidy [2]. WGD is a binary call: a genome is doubled when more than half of it (by length) carries a major-allele copy number ≥ 2 [3][4]. Both computations are exact, specification-driven aggregations over the input segments.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Tumours frequently deviate from the normal diploid (2n) state through aneuploidy and whole-genome doubling. After allele-specific copy-number segmentation (e.g. by ASCAT or FACETS) each genomic segment carries a major-allele copy number n_A and a minor-allele copy number n_B; the segment total copy number is n_A + n_B [2]. The average tumour ploidy summarises the overall copy-number burden, and WGD detection identifies the macro-evolutionary doubling event associated with poor prognosis across cancer types [3].

### 2.2 Core Model

**Average ploidy.** "The average ploidy, PloidyTum, is the average total copy number of all genomic segments weighted by segment length" [1]:

ψ = Σ_i (CN_i · L_i) / Σ_i L_i

where CN_i = n_{A,i} + n_{B,i} is the segment total copy number and L_i = End_i − Start_i is the segment length. The originating allele-specific method is ASCAT, which reports a final tumour ploidy on the n-scale (2n = diploid) [2].

**Whole-genome doubling.** WGD is called when the fraction of the genome (by length) with major-allele copy number ≥ 2 strictly exceeds 0.5 [3][4]:

frac_elevated_mcn = Σ_{i: mcn_i ≥ 2} L_i / Σ_i L_i;  WGD ⇔ frac_elevated_mcn > 0.5

where mcn_i = tcn_i − lcn_i is the major-allele copy number (total minus minor) [4]. The reference implementation uses `treshold = 0.5` and the strict comparison `frac_elevated_mcn > treshold` [4].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | ψ > 0 for any non-empty valid genome with at least one positive copy number | weighted mean of non-negative CN with positive total length [1] |
| INV-02 | a genome of pure 1:1 segments has ψ = 2.0 exactly | every CN_i = 2, so the weighted mean is 2 (n-scale 2n diploid) [1][2] |
| INV-03 | min_i CN_i ≤ ψ ≤ max_i CN_i | a length-weighted mean lies within the value range [1] |
| INV-04 | WGD = true ⇔ (Σ length where major CN ≥ 2) / (Σ length) > 0.5 | direct definition; strict threshold [3][4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| segments | `IEnumerable<AlleleSpecificSegment>` | required | Allele-specific copy-number segments; total CN = Major+Minor, length = End−Start | non-null, non-empty; each segment End > Start, Major ≥ 0, Minor ≥ 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `EstimatePloidy` → ψ | `double` | Length-weighted average ploidy on the n-scale (2n = diploid) |
| `DetectWholeGenomeDoubling` → wgd | `bool` | `true` when > 50% of the genome (by length) has major CN ≥ 2 |

### 3.3 Preconditions and Validation

Coordinates are half-open [Start, End) with length End − Start in base pairs (per the shared `AlleleSpecificSegment`). `segments` null → `ArgumentNullException`. Empty segment set, a segment with End ≤ Start (length ≤ 0), or a negative copy number → `ArgumentException` (ploidy and the WGD fraction are undefined for an empty/zero-length genome). Both methods share the same validation (`ValidateSegment`).

## 4. Algorithm

### 4.1 High-Level Steps

1. For each segment, validate (End > Start, non-negative CN) and accumulate the length and (for ploidy) the product CN_i · L_i; (for WGD) accumulate length where major CN ≥ 2.
2. If total length is 0 (empty input), reject.
3. Ploidy: return Σ(CN·L) / Σ(L). WGD: return (elevated length / total length) > 0.5.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Major-CN-elevation cutoff: mcn ≥ 2 (`WholeGenomeDoublingMajorCopyNumber = 2`) [4].
- WGD fraction threshold: strict > 0.5 (`WholeGenomeDoublingFractionThreshold = 0.5`) [3][4].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| EstimatePloidy / DetectWholeGenomeDoubling | O(n) | O(1) | single pass over n segments; no auxiliary storage |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.EstimatePloidy(IEnumerable<AlleleSpecificSegment>)`: length-weighted average ploidy ψ.
- `OncologyAnalyzer.DetectWholeGenomeDoubling(IEnumerable<AlleleSpecificSegment>)`: WGD flag via the facets-suite major-CN≥2 / >50% rule.

### 5.2 Current Behavior

Both methods stream the input in a single pass and reuse the existing `AlleleSpecificSegment` record and `ValidateSegment` helper introduced for ONCO-LOH-001 / ONCO-HRD-001, so segment validation and the total-copy-number semantics (Major+Minor) are shared across the oncology copy-number units. This is not a search/matching operation, so the repository suffix tree is not applicable.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Average ploidy ψ = Σ(CN_i · L_i) / Σ(L_i), CN_i = Major+Minor — Patchwork length-weighted mean of total copy number [1].
- WGD ⇔ fraction of genome with major CN ≥ 2 (mcn = tcn − lcn) strictly > 0.5 — facets-suite `is_genome_doubled` (PMID 30013179) [3][4].

**Intentionally simplified:**

- WGD fraction denominator: the cited reference divides the elevated length by the **autosomal** genome length from a chromosome-size table; with no chromosome-size table in scope this implementation uses the total length of the supplied segments (the interrogated genome). **Consequence:** identical results when the caller supplies autosomal segments covering the genome; the ≥2 / >0.5 rule is unchanged.

**Not implemented:**

- Inference of allele-specific copy numbers themselves (grid search over purity/ploidy); **users should rely on:** an upstream segmenter (ASCAT/FACETS) to produce `AlleleSpecificSegment` inputs.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | WGD fraction over supplied-segment length, not an external autosomal genome size | Assumption | Fraction denominator differs if caller passes a partial genome | accepted | Matches reference exactly for whole-autosome inputs; see 5.3 |
| 2 | Registry lists `DetectWholeGenomeDoubling(ploidy)` (scalar); canonical method takes segments | Deviation | API shape differs from registry stub | accepted | The cited WGD definition (major CN ≥ 2 over >50% genome) requires per-segment data, not a scalar ploidy |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty segment set | `ArgumentException` | weighted mean / fraction undefined (Σ L = 0) [1] |
| Segment End ≤ Start | `ArgumentException` | non-positive length is invalid input |
| Negative copy number | `ArgumentException` | invalid input |
| All 1:1 segments | ψ = 2.0; WGD = false | total CN 2 but major CN 1 < 2 [1][4] |
| Exactly 50% major CN ≥ 2 | WGD = false | strict `>` 0.5 [4] |
| 2:0 (LOH) over >50% | WGD = true | doubling uses major (not total) CN [4] |

### 6.2 Limitations

Average ploidy and WGD are summary statistics over whatever segments are supplied; they do not infer copy numbers, do not separate clonal/subclonal copy-number states, and (per assumption 1) assume the supplied segments represent the genome of interest. The n-scale interpretation (2n = diploid) presumes a diploid reference [2].

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical walk-through:** segments with total CN 2 (1:1, 100 Mb), 4 (2:2, 100 Mb), 3 (2:1, 50 Mb):
ψ = (2·100 + 4·100 + 3·50) Mb / (100+100+50) Mb = 750 / 250 = **3.0** [1].
WGD: segments B (major 2) and C (major 2) total 150 Mb of 250 Mb = 0.60 > 0.5 → **doubled** [4].

**API usage example:**

```csharp
var segments = new[]
{
    new OncologyAnalyzer.AlleleSpecificSegment("1", 0, 100_000_000, 1, 1),
    new OncologyAnalyzer.AlleleSpecificSegment("2", 0, 100_000_000, 2, 2),
    new OncologyAnalyzer.AlleleSpecificSegment("3", 0,  50_000_000, 2, 1),
};
double ploidy = OncologyAnalyzer.EstimatePloidy(segments);          // 3.0
bool wgd = OncologyAnalyzer.DetectWholeGenomeDoubling(segments);    // true
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_EstimatePloidy_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_EstimatePloidy_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [ONCO-PLOIDY-001-Evidence.md](../../../docs/Evidence/ONCO-PLOIDY-001-Evidence.md)
- Related algorithms: [Copy_Number_Alteration_Classification](Copy_Number_Alteration_Classification.md), [HRD_Score](HRD_Score.md)

## 8. References

1. Mayrhofer M, et al. Patchwork: allele-specific copy number analysis of whole-genome sequenced tumor tissue. *Genome Biology* (PMC4053982). https://pmc.ncbi.nlm.nih.gov/articles/PMC4053982/
2. Van Loo P, Nordgard SH, Lingjærde OC, et al. 2010. Allele-specific copy number analysis of tumors. *PNAS* 107(39):16910–16915. https://doi.org/10.1073/pnas.1009843107
3. Bielski CM, Zehir A, Penson AV, et al. 2018. Genome doubling shapes the evolution and prognosis of advanced cancers. *Nature Genetics* 50(8):1189–1195. https://doi.org/10.1038/s41588-018-0165-1
4. facets-suite (MSKCC). `R/copy-number-scores.R`, `is_genome_doubled` (treshold = 0.5, mcn = tcn − lcn, PMID 30013179). https://github.com/mskcc/facets-suite/blob/master/R/copy-number-scores.R
</content>
