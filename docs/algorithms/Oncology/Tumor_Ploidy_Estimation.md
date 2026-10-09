# Tumor Ploidy Estimation

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-PLOIDY-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Production |
| Last Reviewed | 2026-10-09 |

## 1. Overview

Estimates the average ploidy of a tumour genome from allele-specific copy-number segments and classifies whether the genome has undergone whole-genome doubling (WGD). Average ploidy ψ is the segment-length-weighted mean of per-segment total copy number [1]; on the n-scale a pure-diploid genome has ψ = 2.0 and elevated values (e.g. >2.7n) mark aneuploidy [2]. WGD is a binary call: a genome is doubled when more than half of it (by length) carries a major-allele copy number ≥ 2 [3][4]. Both computations are exact, specification-driven aggregations over the input segments.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Tumours frequently deviate from the normal diploid (2n) state through aneuploidy and whole-genome doubling. After allele-specific copy-number segmentation (e.g. by ASCAT or FACETS) each genomic segment carries a major-allele copy number n_A and a minor-allele copy number n_B; the segment total copy number is n_A + n_B [2]. The average tumour ploidy summarises the overall copy-number burden, and WGD detection identifies the macro-evolutionary doubling event associated with poor prognosis across cancer types [3].

### 2.2 Core Model

**Average ploidy.** "The average ploidy, PloidyTum, is the average total copy number of all genomic segments weighted by segment length" [1]:

ψ = Σ_i (CN_i · L_i) / Σ_i L_i

where CN_i = n_{A,i} + n_{B,i} is the segment total copy number and L_i = End_i − Start_i is the segment length. The originating allele-specific method is ASCAT, which reports a final tumour ploidy on the n-scale (2n = diploid) [2].

**Whole-genome doubling.** WGD is called when more than half of the **autosomal genome** (chromosomes 1–22) has major-allele copy number ≥ 2 [3][4]:

frac_elevated_mcn = Σ_{i: mcn_i ≥ 2, chrom_i ∈ 1..22} L_i / G_autosomal;  WGD ⇔ frac_elevated_mcn > 0.5

where mcn_i = tcn_i − lcn_i is the major-allele copy number (total minus the lesser/minor allele, i.e. max of the two allele copy numbers). The reference implementation (facets-suite `copy-number-scores.R`) computes `autosomal_genome = sum(chrom_info$size[chr %in% 1:22])` with `chrom_info = get_sample_genome(segs, genome)`, whose per-chromosome `size = max(end) − min(start)` is the **interrogated span of the sample's own segments** (the `genome` build only supplies centromeres), and uses `treshold = 0.5` with the strict comparison `frac_elevated_mcn > treshold` [4]. Two denominators are therefore provided:

- **facets-suite (default)** (`DetectWholeGenomeDoubling(segments)`, identical to `DetectWholeGenomeDoublingFromSuppliedLength`): G_autosomal = Σ over autosomes present of (max End − min Start); gaps between a chromosome's first and last segment count in the denominator only. This is the reference implementation's call and is the default since review 2026-09 F32.
- **reference-assembly option** (`DetectWholeGenomeDoubling(segments, ReferenceGenome)`, explicit): G_autosomal = Σ_{c=1..22} of the reference chromosome size — 2,875,001,522 bp (GRCh38) / 2,881,033,286 bp (GRCh37), UCSC `*.chrom.sizes` [5]. Robust to partial-genome inputs, but not the facets-suite call: the two agree only when the segments span every autosome end-to-end (R-confirmed, Evidence "F32").

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | ψ > 0 for any non-empty valid genome with at least one positive copy number | weighted mean of non-negative CN with positive total length [1] |
| INV-02 | a genome of pure 1:1 segments has ψ = 2.0 exactly | every CN_i = 2, so the weighted mean is 2 (n-scale 2n diploid) [1][2] |
| INV-03 | min_i CN_i ≤ ψ ≤ max_i CN_i | a length-weighted mean lies within the value range [1] |
| INV-04 | WGD = true ⇔ (Σ autosomal length where major CN ≥ 2) / G_autosomal > 0.5 | direct definition; strict threshold; G_autosomal = interrogated autosomal span (default, facets-suite) or the reference chromosome-size table (explicit option) [3][4][5] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| segments | `IEnumerable<AlleleSpecificSegment>` | required | Allele-specific copy-number segments; total CN = Major+Minor, length = End−Start | non-null; each segment End > Start, Major ≥ 0, Minor ≥ 0. `EstimatePloidy` requires non-empty (ψ undefined for Σ L = 0); default `DetectWholeGenomeDoubling` requires ≥ 1 autosomal segment (facets 0/0 → NA); the `ReferenceGenome` overload accepts empty (fixed reference denominator) |
| genome | `ReferenceGenome` | — (overload) | Selects the explicit reference-assembly WGD denominator (GRCh38 or GRCh37 autosomal chromosome-size table); omit for the facets-suite default | enum value; `DetectWholeGenomeDoubling(segments, genome)` only |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `EstimatePloidy` → ψ | `double` | Length-weighted average ploidy on the n-scale (2n = diploid) |
| `DetectWholeGenomeDoubling` → wgd | `bool` | `true` when > 50% of the interrogated autosomal genome (default; or of the reference autosomal genome with a `ReferenceGenome`) has major CN ≥ 2 |

### 3.3 Preconditions and Validation

Coordinates are half-open [Start, End) with length End − Start in base pairs (per the shared `AlleleSpecificSegment`). `segments` null → `ArgumentNullException`. A segment with End ≤ Start (length ≤ 0) or a negative copy number → `ArgumentException`; both methods share the same per-segment validation (`ValidateSegment`). `EstimatePloidy` additionally rejects an empty segment set (ψ is undefined for Σ L = 0). The default `DetectWholeGenomeDoubling(segments)` and `DetectWholeGenomeDoublingFromSuppliedLength` (facets-suite) reject input with no autosomal segment (empty included) with `ArgumentException` — the interrogated-span denominator is 0 (R: 0/0 → `NA`). `DetectWholeGenomeDoubling(segments, ReferenceGenome)` divides by the fixed reference autosomal genome length, so an empty set yields numerator 0 → fraction 0 → `false` (no exception); an undefined `ReferenceGenome` value → `ArgumentOutOfRangeException`.

## 4. Algorithm

### 4.1 High-Level Steps

1. For each segment, validate (End > Start, non-negative CN) and accumulate the length and (for ploidy) the product CN_i · L_i; (for WGD) accumulate length where major CN ≥ 2 **on autosomes (chr1–22 only)**.
2. Ploidy: if total length is 0 (empty input), reject. WGD: G_autosomal = Σ per-autosome interrogated span max(End) − min(Start) (default, facets-suite `get_sample_genome`; 0 → reject), or the embedded chromosome-size table of the selected `ReferenceGenome` (explicit option).
3. Ploidy: return Σ(CN·L) / Σ(L). WGD: return (elevated autosomal length / G_autosomal) > 0.5.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Major-CN-elevation cutoff: mcn ≥ 2 (`WholeGenomeDoublingMajorCopyNumber = 2`) [4].
- WGD fraction threshold: strict > 0.5 (`WholeGenomeDoublingFractionThreshold = 0.5`) [3][4].
- Major CN: mcn = tcn − lcn = max(Major, Minor) (`IsElevatedMajorCopyNumber`), independent of allele labelling [4].
- Length sums are accumulated in `double` (facets-suite `as.numeric`), so Σ L cannot overflow Int64.
- facets-suite interrogated denominator: per-autosome span max(End) − min(Start) (`get_sample_genome`); no autosomal segment → undefined (R `NA`) → `ArgumentException` [4].
- Reference autosomal genome length: embedded chromosome-size tables `GRCh38AutosomeLengths` / `GRCh37AutosomeLengths` (chr1–22 from UCSC `*.chrom.sizes`); summed by `GetAutosomalGenomeLength` (GRCh38 = 2,875,001,522 bp, GRCh37 = 2,881,033,286 bp) [5].
- Autosome restriction: numerator counts only segments whose chromosome parses to 1–22 (`chrom %in% 1:22`), accepting both "7" and "chr7" forms [4][5].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| EstimatePloidy / DetectWholeGenomeDoubling | O(n) | O(1) | single pass over n segments; no auxiliary storage |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.CopyNumberPloidy.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.CopyNumberPloidy.cs)

- `OncologyAnalyzer.EstimatePloidy(IEnumerable<AlleleSpecificSegment>)`: length-weighted average ploidy ψ.
- `OncologyAnalyzer.EstimatePloidy(IEnumerable<AlleleSpecificSegment>, IEnumerable<int> probeCounts)`: ASCAT probe-count-weighted ploidy Σ(CN_i·n_i)/Σn_i (runASCAT) [6]; also used by `FitPurityPloidy` for `PurityPloidyFit.Ploidy`.
- `OncologyAnalyzer.ComputeAscatGenomeMetrics(IEnumerable<AlleleSpecificSegment>)` → `AscatGenomeMetrics(ModeMinorAllele, ModeMajorAllele, WgdStatus, GenomicInstability, LossOfHeterozygosity)`: ASCAT `ascat.metrics` `mode_minA`/`mode_majA`/`WGD` (NA/0/1/1+)/`GI`/`LOH` [6].
- `OncologyAnalyzer.DetectWholeGenomeDoubling(IEnumerable<AlleleSpecificSegment>)` (default): WGD flag exactly as facets-suite `is_genome_doubled(segs, get_sample_genome(segs))` — denominator Σ over autosomes of the interrogated span (max End − min Start); non-autosomal segments ignored (F32).
- `OncologyAnalyzer.DetectWholeGenomeDoubling(IEnumerable<AlleleSpecificSegment>, ReferenceGenome)`: explicit option — same numerator/threshold against the reference autosomal chromosome-size table.
- `OncologyAnalyzer.DetectWholeGenomeDoublingFromSuppliedLength(IEnumerable<AlleleSpecificSegment>)`: the same facets-suite call as the default (kept for existing callers).
- `OncologyAnalyzer.GetAutosomeLengths(ReferenceGenome)` / `GetAutosomalGenomeLength(ReferenceGenome)`: the embedded reference chromosome-size table and its autosomal sum.

### 5.2 Current Behavior

Both methods stream the input in a single pass and reuse the existing `AlleleSpecificSegment` record and `ValidateSegment` helper introduced for ONCO-LOH-001 / ONCO-HRD-001, so segment validation and the total-copy-number semantics (Major+Minor) are shared across the oncology copy-number units. This is not a search/matching operation, so the repository suffix tree is not applicable.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Average ploidy ψ = Σ(CN_i · L_i) / Σ(L_i), CN_i = Major+Minor — Patchwork length-weighted mean of total copy number [1].
- WGD ⇔ fraction of genome with major CN ≥ 2 (mcn = tcn − lcn) strictly > 0.5 — facets-suite `is_genome_doubled` (PMID 30013179) [3][4].
- `DetectWholeGenomeDoubling(segments)` (default, F32) = `DetectWholeGenomeDoublingFromSuppliedLength`: facets-suite `autosomal_genome = sum(chrom_info$size[chr %in% 1:22])` with `chrom_info = get_sample_genome(segs)` (size = max(end) − min(start)), numerator restricted to autosomes — line-by-line equal to the R code (Python port cross-check, review 2026-09; the real R functions sourced from facets-suite master on 6 genomes incl. gapped / telomere-trimmed / chrX inputs, review 2026-09 F32 — Evidence) [4].
- `DetectWholeGenomeDoubling(segments, ReferenceGenome)` (explicit option): same numerator/threshold; denominator = reference autosomal length from the embedded UCSC `*.chrom.sizes` tables (GRCh38/GRCh37) [5] — see §5.4 #1.
- `EstimatePloidy(segments, probeCounts)` (review 2026-09, F31): ASCAT runASCAT `ploidy = sum((nA+nB)*s[,"length"])/sum(s[,"length"])` (l. 283, length = #probes of each `make_segments` segment) = the reported `mean(nA+nB, na.rm=TRUE)` over probes (l. 98) [6]. R cross-check: `mean(rep(c(2,4,3), c(1000,200,800)))` = 2.6000000000000001 (bp-weighted default on the same segments: 3.0); 3000 random genomes bit-identical.
- `ComputeAscatGenomeMetrics` (review 2026-09, F30): ASCAT `ascat.metrics` verbatim [6] — autosomes only (`setdiff(chrs, sexchromosomes)`, X/Y); `modeAllele`: `round`, cap 5, weight `(endpos−startpos)/1e6` (no +1) summed per value, stable decreasing `order` + `which.max` ⇒ ties go to the smaller copy number; WGD: mode_majA 0 → NA, 1 → "0", 2 → "1", 3–5 → "1+"; `computeGIscore`: `round(1 − Σsize[nMajor=b ∧ nMinor=b]/Σsize, 4)`, size = `endpos−startpos+1`, b = 1 (WGD 0) or 2 (WGD 1 / 1+); LOH = `round(Σsize[nMinor=0]/Σsize, 4)`. Bit-identical to the R function on 9 hand-built genomes (tie, NA, cap, sex-chromosome cases) and 3000 random genomes (108 exact mode ties).

**Intentionally simplified:**

- (none). Note: the default `EstimatePloidy` weights by base pairs (Patchwork [1]); ASCAT's probe-count weighting is implemented as the `probeCounts` overload (F31; identical for uniform probe density). ASCAT's own WGD metric (`ascat.metrics`: mode of nMajor) is implemented as `ComputeAscatGenomeMetrics` (F30), alongside the facets-suite/Bielski rule.

**Not implemented:**

- (none for this unit's scope.) Allele-specific copy numbers themselves are inferred by ONCO-ASCAT-001, not here: `SegmentAlleleSpecificAspcf` (ASCAT ASPCF segmentation) → `FitPurityPloidy` (ASCAT runASCAT purity/ploidy grid search; `PurityPloidyFit.Segments` are `AlleleSpecificSegment`s ready for `EstimatePloidy` / `DetectWholeGenomeDoubling`). An external segmenter (ASCAT/FACETS) can still supply the segments.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | WGD denominator: default = facets-suite interrogated span; reference-assembly length only via the explicit `ReferenceGenome` overload | Resolved | Default call now equals facets-suite `is_genome_doubled(segs, get_sample_genome(segs))`; the reference-assembly option differs when segments do not span each autosome end-to-end | resolved 2026-10-09 (review 2026-09 F32) | The 2026-06-22 reference-table default rested on a misreading of facets-suite (`chrom_info` is `get_sample_genome(segs)`, F9); F32 made the reference's denominator the default (R-confirmed) and kept the reference-assembly denominator as an explicit option |
| 2 | Registry lists `DetectWholeGenomeDoubling(ploidy)` (scalar); canonical method takes segments | Deviation | API shape differs from registry stub | accepted | The cited WGD definition (major CN ≥ 2 over >50% genome) requires per-segment data, not a scalar ploidy |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Only non-autosomal segments | default `DetectWholeGenomeDoubling` / `…FromSuppliedLength` → `ArgumentException` | facets-suite 0/0 → `NA` [4] |
| Gap between segments of one chromosome | default (facets): gap counted in denominator (60 Mb elevated over a 140 Mb span → 0.43 → false) | `get_sample_genome` size = max(end) − min(start) [4] |
| Allele labels swapped (Major < Minor) | mcn = max(Major, Minor) | mcn = tcn − lcn, lcn = lesser allele [4] |
| Empty segment set | `EstimatePloidy` → `ArgumentException`; default `DetectWholeGenomeDoubling` → `ArgumentException`; `ReferenceGenome` overload → `false` | ψ undefined (Σ L = 0); facets 0/0 → NA; numerator 0 over a fixed reference denominator [1][4] |
| Segment End ≤ Start | `ArgumentException` | non-positive length is invalid input |
| Negative copy number | `ArgumentException` | invalid input |
| All 1:1 autosomal segments | ψ = 2.0; WGD = false | total CN 2 but major CN 1 < 2 [1][4] |
| Exactly half the GRCh38 autosomal genome at major CN ≥ 2 (`ReferenceGenome` overload) | WGD = false | strict `>` 0.5 against G_autosomal [4][5] |
| Just over half the reference autosomal genome at major CN ≥ 2 (`ReferenceGenome` overload) | WGD = true | frac > 0.5 [4][5] |
| Small fully-amplified region (e.g. 100 Mb, genome not tiled) | default (facets): WGD = true (100/100 Mb); `ReferenceGenome` overload: false | facets divides by the interrogated span; 100 Mb / 2.875 Gb ≈ 0.035 < 0.5 against the reference [4][5] |
| chrX/chrY amplified, no autosomal elevation | WGD = false | numerator restricted to autosomes (chr1–22) [4][5] |
| 2:0 (LOH) over >50% of the autosomal genome | WGD = true | doubling uses major (not total) CN [4] |

### 6.2 Limitations

Average ploidy and WGD are summary statistics over the supplied segments; they do not infer copy numbers and do not separate clonal/subclonal copy-number states. `EstimatePloidy` is still a length-weighted mean over whatever segments are supplied (it assumes they represent the genome of interest). The default WGD call divides by the interrogated autosomal span exactly as facets-suite does, so it judges only the genome the segments cover (a partial-genome input is judged on its own span); callers who want partial inputs judged against the whole autosomal genome pass the `ReferenceGenome` matching their coordinates to the explicit overload. The n-scale interpretation (2n = diploid) presumes a diploid reference [2].

## 7. Examples and Related Material

### 7.1 Worked Example

**Ploidy walk-through:** segments with total CN 2 (1:1, 100 Mb), 4 (2:2, 100 Mb), 3 (2:1, 50 Mb):
ψ = (2·100 + 4·100 + 3·50) Mb / (100+100+50) Mb = 750 / 250 = **3.0** [1].

**WGD walk-through (default = facets-suite span):** the three segments above span chr1 100 + chr2 100 + chr3 50 = 250 Mb; 150 Mb has major CN ≥ 2 → 150/250 = 0.60 > 0.5 → **doubled** [4].

**WGD walk-through (explicit reference-assembly option):** the GRCh38 autosomal genome is G = 2,875,001,522 bp (half = 1,437,500,761 bp) [5]. A single autosomal segment of 1,437,500,762 bp at major CN 2 gives 1,437,500,762 / G > 0.5 → **doubled**; at exactly 1,437,500,761 bp the fraction is 0.5, not > 0.5 → **not doubled** (strict) [4]. The same 150 Mb of major-CN≥2 from the ploidy example is only 150 Mb / 2.875 Gb ≈ 0.052 → **not doubled**, because the segments do not tile the genome.

**API usage example:**

```csharp
var segments = new[]
{
    new OncologyAnalyzer.AlleleSpecificSegment("1", 0, 100_000_000, 1, 1),
    new OncologyAnalyzer.AlleleSpecificSegment("2", 0, 100_000_000, 2, 2),
    new OncologyAnalyzer.AlleleSpecificSegment("3", 0,  50_000_000, 2, 1),
};
double ploidy = OncologyAnalyzer.EstimatePloidy(segments);          // 3.0
// Default = facets-suite (interrogated spans chr1 100 + chr2 100 + chr3 50 Mb): 150/250 = 0.60 > 0.5 → doubled.
bool wgd = OncologyAnalyzer.DetectWholeGenomeDoubling(segments);    // true
// Explicit reference-assembly option: 150 Mb / 2.875 Gb (GRCh38 autosomes) ≈ 5% → not doubled.
bool wgdRef = OncologyAnalyzer.DetectWholeGenomeDoubling(segments, OncologyAnalyzer.ReferenceGenome.GRCh38); // false
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_EstimatePloidy_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_EstimatePloidy_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [ONCO-PLOIDY-001-Evidence.md](../../../docs/Evidence/ONCO-PLOIDY-001-Evidence.md)
- Related algorithms: [Copy_Number_Alteration_Classification](Copy_Number_Alteration_Classification.md), [HRD_Score](HRD_Score.md)

## 8. References

1. Mayrhofer M, et al. Patchwork: allele-specific copy number analysis of whole-genome sequenced tumor tissue. *Genome Biology* (PMC4053982). https://pmc.ncbi.nlm.nih.gov/articles/PMC4053982/
2. Van Loo P, Nordgard SH, Lingjærde OC, et al. 2010. Allele-specific copy number analysis of tumors. *PNAS* 107(39):16910–16915. https://doi.org/10.1073/pnas.1009843107
3. Bielski CM, Zehir A, Penson AV, et al. 2018. Genome doubling shapes the evolution and prognosis of advanced cancers. *Nature Genetics* 50(8):1189–1195. https://doi.org/10.1038/s41588-018-0165-1
4. facets-suite (MSKCC). `R/copy-number-scores.R`, `is_genome_doubled` (treshold = 0.5, mcn = tcn − lcn, `autosomal_genome = sum(chrom_info$size[chr %in% 1:22])`, PMID 30013179), `get_sample_genome` (size = max(end) − min(start) per chromosome), `parse_segs`, `calculate_fraction_cna`. https://raw.githubusercontent.com/mskcc/facets-suite/master/R/copy-number-scores.R (re-read 2026-09-28)
5. UCSC Genome Browser. `hg38.chrom.sizes` (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.chrom.sizes) and `hg19.chrom.sizes` (https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes); GRCh38 chromosome lengths cross-verified against Ensembl REST GRCh38.p14 (https://rest.ensembl.org/info/assembly/homo_sapiens). Accessed 2026-06-22.
</content>
6. ASCAT (VanLoo-lab). `R/ascat.metrics.R` (`ascat.metrics`, `modeAllele`, `computeGIscore`) and `R/ascat.runAscat.R` (`ploidy = sum((nA+nB)*s[,"length"])/sum(s[,"length"])`, l. 283; `ploidy = mean(nA+nB, na.rm=TRUE)`, l. 98). https://raw.githubusercontent.com/VanLoo-lab/ascat/master/ASCAT/R/ascat.metrics.R, …/ascat.runAscat.R (read and executed in R 4.3.3, 2026-10-09)
