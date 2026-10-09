# Test Specification: ONCO-PLOIDY-001

**Test Unit ID:** ONCO-PLOIDY-001
**Area:** Oncology
**Algorithm:** Tumor Ploidy Estimation (length-weighted mean segment copy number) + Whole-Genome-Doubling detection
**Status:** ☐ In Progress (limitation fix — pending re-validation)
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-09-28

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Patchwork (Genome Biology) — verbatim ploidy definition | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC4053982/ | 2026-06-14 |
| 2 | Van Loo et al. ASCAT (PNAS 2010) | 1 | https://doi.org/10.1073/pnas.1009843107 | 2026-06-14 |
| 3 | Bielski et al. (Nature Genetics 2018) — WGD | 1 | https://doi.org/10.1038/s41588-018-0165-1 | 2026-06-14 |
| 4 | facets-suite `copy-number-scores.R` `is_genome_doubled` (reference impl, PMID 30013179) | 3 | https://github.com/mskcc/facets-suite/blob/master/R/copy-number-scores.R | 2026-06-22 |
| 5 | UCSC `hg38.chrom.sizes` / `hg19.chrom.sizes` (reference chromosome-size tables); Ensembl GRCh38.p14 cross-verification | 5 | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.chrom.sizes ; https://rest.ensembl.org/info/assembly/homo_sapiens | 2026-06-22 |

### 1.2 Key Evidence Points

1. Average tumour ploidy is "the average total copy number of all genomic segments weighted by segment length" → ψ = Σ(CN_i · L_i) / Σ(L_i) — Patchwork, PMC4053982.
2. Ploidy is reported on the n-scale (2n = diploid); ">2.7n" marks aneuploidy / near-triploid genomes — Van Loo et al. 2010, PNAS abstract.
3. WGD is called when the autosome-restricted fraction of genome with **major copy number ≥ 2** is strictly greater than 0.5: `frac_elevated_mcn > treshold` (treshold = 0.5) — facets-suite `is_genome_doubled` (PMID 30013179).
4. Major copy number `mcn = tcn - lcn` (total − minor); WGD uses the major allele CN ≥ 2, not total CN ≥ 2 — facets-suite `parse_segs`.
5. facets-suite WGD denominator: `autosomal_genome = sum(chrom_info$size[chr %in% 1:22])` with `chrom_info = get_sample_genome(segs)` — per-autosome interrogated span max(end) − min(start); numerator restricted to autosomes (`chrom %in% 1:22`) — implemented exactly by `DetectWholeGenomeDoublingFromSuppliedLength`. The canonical `DetectWholeGenomeDoubling` instead divides by the reference autosomal length GRCh38 Σ(chr1–22) = 2,875,001,522 bp / GRCh37 = 2,881,033,286 bp (UCSC `*.chrom.sizes`) — a documented deviation (review 2026-09, F9).

### 1.3 Documented Corner Cases

- Empty segment set → Σ(L) = 0, ploidy undefined (Patchwork weighted mean); WGD against the fixed reference denominator returns false (numerator 0).
- Segment with Length ≤ 0 or negative copy number → invalid input.
- WGD threshold is strict (`>` 0.5): exactly half the reference autosomal genome at major CN ≥ 2 is NOT doubled (facets-suite).
- WGD uses **major** CN: an all-1:1 genome (total CN 2) is NOT doubled (facets-suite `mcn >= 2`).
- WGD numerator is autosome-restricted: chrX/chrY/contig segments do not contribute (facets-suite `chrom %in% 1:22`); a fully-amplified region that does not tile the genome is NOT doubled (reference denominator removes supplied-segment bias).

### 1.4 Known Failure Modes / Pitfalls

1. Using a plain (unweighted) mean of per-segment copy numbers instead of length-weighting it — Patchwork ("weighted by segment length").
2. Calling WGD on total CN ≥ 2 instead of major CN ≥ 2 (would mis-call balanced diploids) — facets-suite.
3. Using `≥ 0.5` instead of strict `> 0.5` for the fraction — facets-suite `> treshold`.
4. Using the supplied segments' total length as the WGD denominator instead of the reference autosomal genome length — would over-call WGD for partial-genome inputs (the limitation this fix resolves) — facets-suite `autosomal_genome`.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `EstimatePloidy(IEnumerable<AlleleSpecificSegment>)` | OncologyAnalyzer | Canonical | ψ = Σ(CN·L)/Σ(L), CN = Major+Minor |
| `DetectWholeGenomeDoubling(IEnumerable<AlleleSpecificSegment>, ReferenceGenome=GRCh38)` | OncologyAnalyzer | Canonical | facets-suite rule: frac(autosomal major CN ≥ 2 length) / reference autosomal genome > 0.5 |
| `DetectWholeGenomeDoublingFromSuppliedLength(IEnumerable<AlleleSpecificSegment>)` | OncologyAnalyzer | Variant | facets-suite exact: denominator = Σ autosomal interrogated span (get_sample_genome) |
| `GetAutosomeLengths(ReferenceGenome)` / `GetAutosomalGenomeLength(ReferenceGenome)` | OncologyAnalyzer | Canonical | embedded reference chromosome-size table + autosomal sum |
| `EstimatePloidy(IEnumerable<AlleleSpecificSegment>, IEnumerable<int> probeCounts)` | OncologyAnalyzer | Variant | ASCAT runASCAT probe-count-weighted ψ = Σ(CN·n)/Σn (F31) |
| `ComputeAscatGenomeMetrics(IEnumerable<AlleleSpecificSegment>)` | OncologyAnalyzer | Variant | ASCAT `ascat.metrics` mode_minA / mode_majA / WGD (NA/0/1/1+) / GI / LOH (F30) |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | ploidy > 0 for any non-empty valid segment set with at least one positive copy number | Yes | Patchwork weighted mean; registry invariant |
| INV-2 | a genome of pure 1:1 (total CN 2) segments has ploidy exactly 2.0 | Yes | n-scale 2n diploid (ASCAT/Patchwork) |
| INV-3 | ploidy is length-weighted: min(CN_i) ≤ ψ ≤ max(CN_i) | Yes | weighted mean lies within the value range (Patchwork) |
| INV-4 | WGD = true ⇔ (Σ autosomal length where major CN ≥ 2) / G_autosomal > 0.5, G_autosomal from the reference chromosome-size table | Yes | facets-suite `is_genome_doubled` |
| INV-5 | embedded GRCh38/GRCh37 autosome length tables equal the authoritative UCSC `*.chrom.sizes` values exactly; sums = 2,875,001,522 / 2,881,033,286 bp | Yes | UCSC chrom.sizes; Ensembl GRCh38.p14 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Ploidy worked example | CN 2 (1:1)/4 (2:2)/3 (2:1), lengths 100/100/50 Mb | ψ = 750M/250M = 3.0 | Patchwork weighted mean |
| M2 | Pure diploid | all 1:1 segments | ψ = 2.0 exactly | n-scale 2n |
| M3 | Length weighting dominates | long 1:1 (300 Mb) + short 2:2 (10 Mb) | ψ = (2·300+4·10)/310 = 640/310 ≈ 2.0645 | "weighted by segment length" |
| M4 | Single segment | one 2:1 segment (total 3) | ψ = 3.0 | weighted mean of one value |
| M5 | Empty segments → reject | no segments | ArgumentException | Σ(L)=0 undefined |
| M6 | Invalid segment length → reject | End ≤ Start (Length ≤ 0) | ArgumentException | invalid input |
| M7 | Negative copy number → reject | Major or Minor < 0 | ArgumentException | invalid input |
| M8 | WGD just over half of GRCh38 genome → true | autosomal major-CN≥2 length = (G/2)+1 = 1,437,500,762 bp | true | facets-suite > 0.5 vs G_autosomal |
| M9 | WGD exactly half of GRCh38 genome → false | autosomal major-CN≥2 length = G/2 = 1,437,500,761 bp | false | strict `>` 0.5 |
| M10 | WGD just under half of GRCh38 genome → false | length = (G/2)−1 = 1,437,500,760 bp | false | frac < 0.5 |
| M11 | WGD all 1:1 (total 2) → false | every autosomal segment major CN = 1 | false | mcn >= 2 (not total) |
| M12 | WGD small fully-amplified region → false | 100 Mb all major ≥ 2, genome not tiled | false | reference denominator removes supplied-segment bias |
| M13 | WGD invalid/null → reject | Length ≤ 0; negative CN; null | ArgumentException / ArgumentNullException | shared validation |
| M14 | GRCh38 autosome table matches UCSC | `GetAutosomeLengths(GRCh38)` | equals 22 UCSC hg38.chrom.sizes values exactly | UCSC hg38.chrom.sizes |
| M15 | GRCh37 autosome table matches UCSC | `GetAutosomeLengths(GRCh37)` | equals 22 UCSC hg19.chrom.sizes values exactly | UCSC hg19.chrom.sizes |
| M16 | autosomal genome sums | `GetAutosomalGenomeLength` | GRCh38 = 2,875,001,522; GRCh37 = 2,881,033,286 bp | Σ(chr1–22) |
| M17 | GRCh37 selector uses hg19 denominator | (G_hg19/2)+1 bp at major ≥ 2, both builds | true under GRCh37 and GRCh38 | build-dependent denominator |
| M18 | WGD empty set → false | no segments, reference denominator | false (numerator 0) | fixed reference denominator |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | WGD LOH counts as elevated | 2:0 segment (major 2, minor 0) over half the reference genome | true | major CN, not heterozygosity |
| S2 | Ploidy with a CN-0 (homozygous deletion) segment | 0 (0:0)/4 (2:2) equal lengths | ψ = 2.0 | weighted mean includes zeros |
| S3 | WGD excludes sex chromosomes | chrX/chrY amplified, no autosomal elevation | false | facets-suite `chrom %in% 1:22` |
| S4 | WGD recognises "chr"-prefixed autosomes | chr7 over half the reference genome at major ≥ 2 | true | autosome parser accepts chr-prefix |
| L1 | facets-exact WGD 60% → true | 60% of interrogated span at major CN ≥ 2 | true | `DetectWholeGenomeDoublingFromSuppliedLength` |
| L2 | facets-exact WGD exactly 50% → false | half interrogated span elevated | false | strict `>` 0.5 |
| L3 | facets-exact WGD empty → reject | no segments | ArgumentException | interrogated denominator 0 (R `NA`) |
| L5 | facets-exact gap inside chromosome | (1,0,60M,2:2),(1,100M,140M,1:1) | false (0.4286) | get_sample_genome span (F9) |
| L6 | facets-exact excludes chrX from both terms | (1,0,40M,2:2),(1,40M,100M,1:1),(X,0,100M,2:2) | false (0.4) | chrom %in% 1:22 (F9) |
| L7 | facets-exact span starts at min start | (chr1,10M,70M,2:2),(chr1,70M,110M,1:1) | true (0.6) | size = max(end) − min(start) (F9) |
| L8 | facets-exact only non-autosomal → reject | X / chrY only | ArgumentException | 0/0 → NA (F9) |
| S5 | swapped allele labels | Major 1 / Minor 2 | elevated (mcn = 3 − 1 = 2) in both WGD methods | mcn = tcn − lcn (F10) |
| S6 | Σ L beyond Int64 | two 5e18-bp segments | ψ = 2; WGD true (both) | as.numeric sums (F11) |
| L4 | Legacy WGD null → reject | null | ArgumentNullException | guard contract |
| A1–A9 | ASCAT ascat.metrics G1–G9 | Evidence table "ASCAT ascat.metrics" (diploid, WGD, 1+, tie → smaller value, mode 0 → NA, X excluded, cap 5, tie 3/4, +1 size / chr prefix) | R values (mode_minA, mode_majA, WGD, GI, LOH) exactly | F30, R-executed |
| A10 | ASCAT metrics swapped allele labels | Major 0 / Minor 2 | same as 2:0 | nMajor = max |
| A11 | ASCAT metrics only X/Y or null | — | ArgumentException / ArgumentNullException | R zero-length mode errors |
| P1 | ASCAT probe-weighted ploidy | CN 2/4/3, probes 1000/200/800 | 2.6000000000000001 (bp default 3.0 unchanged) | F31, R `mean(rep())` |
| P2 | ASCAT probe-weighted, non-terminating | CN 2/3/5, probes 1/1/1 | 3.3333333333333335 | F31 |
| P3 | Probe-weighted invalid input | count mismatch, 0 probes, empty, null | ArgumentException / ArgumentNullException | guard contract |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Near-triploid genome ploidy | mostly CN-3 segments | ψ ≈ 3 (>2.7n aneuploid direction) | Van Loo aneuploidy direction |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Limitation fix on an existing unit. `EstimatePloidy` is unchanged. `DetectWholeGenomeDoubling` now divides by the reference autosomal genome length (embedded UCSC `*.chrom.sizes`) and adds a `ReferenceGenome` parameter; the prior supplied-segment-length behaviour is preserved as `DetectWholeGenomeDoublingFromSuppliedLength`. New accessors `GetAutosomeLengths` / `GetAutosomalGenomeLength`. Existing canonical test file `OncologyAnalyzer_EstimatePloidy_Tests.cs` is updated; property (`OncologyProperties.cs`) and combinatorial (`OncologyCombinatorialTests.cs`) WGD assertions were re-pointed at the legacy overload (they encode the supplied-length semantics) to keep their oracles valid.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M7, S2, C1 (ploidy) | ✅ Covered | unchanged from prior version |
| M8–M13 (WGD) | 🔁 Rewritten | re-derived against the reference autosomal genome denominator |
| M14–M18, S3–S4 | ❌ Missing → implemented | new reference-table / autosome-restriction / build-selector cases |
| L1–L4 (legacy overload) | ❌ Missing → implemented | smoke verification of supplied-length variant |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_EstimatePloidy_Tests.cs` — all ploidy + WGD + reference-table tests.
- **Remove:** none. Re-point WGD calls in `OncologyProperties.cs` / `OncologyCombinatorialTests.cs` to the legacy overload.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `OncologyAnalyzer_EstimatePloidy_Tests.cs` | canonical | 30 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M8 | 🔁 Rewritten | DetectWholeGenomeDoubling_JustOverHalfOfReferenceGenome_ReturnsTrue | ✅ Done |
| 2 | M9 | 🔁 Rewritten | DetectWholeGenomeDoubling_ExactlyHalfOfReferenceGenome_ReturnsFalse | ✅ Done |
| 3 | M10 | 🔁 Rewritten | DetectWholeGenomeDoubling_JustUnderHalfOfReferenceGenome_ReturnsFalse | ✅ Done |
| 4 | M11 | 🔁 Rewritten | DetectWholeGenomeDoubling_AllBalancedDiploid_ReturnsFalse | ✅ Done |
| 5 | M12 | 🔁 Rewritten | DetectWholeGenomeDoubling_SmallFullyAmplifiedRegion_IsNotDoubledAgainstReferenceGenome | ✅ Done |
| 6 | M13 | 🔁 Rewritten | NonPositiveLength/NegativeCopyNumber/Null_Throws | ✅ Done |
| 7 | M14 | ❌ Missing | GetAutosomeLengths_GRCh38_MatchesUcscChromSizes | ✅ Done |
| 8 | M15 | ❌ Missing | GetAutosomeLengths_GRCh37_MatchesUcscChromSizes | ✅ Done |
| 9 | M16 | ❌ Missing | GetAutosomalGenomeLength_BothBuilds_MatchSummedChromSizes | ✅ Done |
| 10 | M17 | ❌ Missing | DetectWholeGenomeDoubling_GRCh37Selector_UsesHg19Denominator | ✅ Done |
| 11 | M18 | ❌ Missing | DetectWholeGenomeDoubling_EmptySegments_ReturnsFalse | ✅ Done |
| 12 | S1 | 🔁 Rewritten | DetectWholeGenomeDoubling_LohSegmentsOverHalfReference_CountAsElevated | ✅ Done |
| 13 | S3 | ❌ Missing | DetectWholeGenomeDoubling_SexChromosomeSegments_ExcludedFromNumerator | ✅ Done |
| 14 | S4 | ❌ Missing | DetectWholeGenomeDoubling_ChrPrefixedAutosomes_AreRecognised | ✅ Done |
| 15 | L1 | ❌ Missing | DetectWholeGenomeDoublingFromSuppliedLength_SixtyPercentElevated_ReturnsTrue | ✅ Done |
| 16 | L2 | ❌ Missing | DetectWholeGenomeDoublingFromSuppliedLength_ExactlyHalf_ReturnsFalse | ✅ Done |
| 17 | L3 | ❌ Missing | DetectWholeGenomeDoublingFromSuppliedLength_Empty_Throws | ✅ Done |
| 18 | L4 | ❌ Missing | DetectWholeGenomeDoublingFromSuppliedLength_Null_Throws | ✅ Done |

**Total items:** 18 (WGD/reference); ploidy cases M1–M7, S2, C1 unchanged.
**✅ Done:** 18 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | EstimatePloidy_WorkedExample_ReturnsThree |
| M2 | ✅ Covered | EstimatePloidy_PureDiploid_ReturnsTwo |
| M3 | ✅ Covered | EstimatePloidy_LongDiploidShortAmplified_IsLengthWeighted |
| M4 | ✅ Covered | EstimatePloidy_SingleSegment_ReturnsItsTotalCopyNumber |
| M5 | ✅ Covered | EstimatePloidy_EmptySegments_Throws |
| M6 | ✅ Covered | EstimatePloidy_NonPositiveLength_Throws |
| M7 | ✅ Covered | EstimatePloidy_NegativeCopyNumber_Throws |
| M8 | ✅ Covered | DetectWholeGenomeDoubling_JustOverHalfOfReferenceGenome_ReturnsTrue |
| M9 | ✅ Covered | DetectWholeGenomeDoubling_ExactlyHalfOfReferenceGenome_ReturnsFalse |
| M10 | ✅ Covered | DetectWholeGenomeDoubling_JustUnderHalfOfReferenceGenome_ReturnsFalse |
| M11 | ✅ Covered | DetectWholeGenomeDoubling_AllBalancedDiploid_ReturnsFalse |
| M12 | ✅ Covered | DetectWholeGenomeDoubling_SmallFullyAmplifiedRegion_IsNotDoubledAgainstReferenceGenome |
| M13 | ✅ Covered | DetectWholeGenomeDoubling_NonPositiveLength/NegativeCopyNumber/Null_Throws |
| M14 | ✅ Covered | GetAutosomeLengths_GRCh38_MatchesUcscChromSizes |
| M15 | ✅ Covered | GetAutosomeLengths_GRCh37_MatchesUcscChromSizes |
| M16 | ✅ Covered | GetAutosomalGenomeLength_BothBuilds_MatchSummedChromSizes |
| M17 | ✅ Covered | DetectWholeGenomeDoubling_GRCh37Selector_UsesHg19Denominator |
| M18 | ✅ Covered | DetectWholeGenomeDoubling_EmptySegments_ReturnsFalse |
| S1 | ✅ Covered | DetectWholeGenomeDoubling_LohSegmentsOverHalfReference_CountAsElevated |
| S2 | ✅ Covered | EstimatePloidy_WithHomozygousDeletionSegment_IncludesZeros |
| S3 | ✅ Covered | DetectWholeGenomeDoubling_SexChromosomeSegments_ExcludedFromNumerator |
| S4 | ✅ Covered | DetectWholeGenomeDoubling_ChrPrefixedAutosomes_AreRecognised |
| L1 | ✅ Covered | DetectWholeGenomeDoublingFromSuppliedLength_SixtyPercentElevated_ReturnsTrue |
| L2 | ✅ Covered | DetectWholeGenomeDoublingFromSuppliedLength_ExactlyHalf_ReturnsFalse |
| L3 | ✅ Covered | DetectWholeGenomeDoublingFromSuppliedLength_Empty_Throws |
| L4 | ✅ Covered | DetectWholeGenomeDoublingFromSuppliedLength_Null_Throws |
| C1 | ✅ Covered | EstimatePloidy_NearTriploidGenome_ExceedsAneuploidyDirection |
| A1–A9 | ✅ Covered | ComputeAscatGenomeMetrics_MatchesAscatR (AscatMetrics_G1…G9) |
| A10 | ✅ Covered | ComputeAscatGenomeMetrics_SwappedAlleleLabels_SameResult |
| A11 | ✅ Covered | ComputeAscatGenomeMetrics_OnlySexChromosomes_Throws / _Null_Throws |
| P1 | ✅ Covered | EstimatePloidy_ProbeCounts_MatchesAscatR |
| P2 | ✅ Covered | EstimatePloidy_ProbeCounts_NonTerminating_MatchesAscatR |
| P3 | ✅ Covered | EstimatePloidy_ProbeCounts_InvalidInput_Throws |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Per-segment total CN supplied as `AlleleSpecificSegment` (total = Major+Minor; length = End−Start) | EstimatePloidy, DetectWholeGenomeDoubling |

*Resolved 2026-06-22:* the previous Assumption 2 (WGD denominator = supplied-segment length) is removed — the WGD fraction now divides by the reference autosomal genome length from the embedded UCSC `*.chrom.sizes` tables (facets-suite `autosomal_genome`), selected by `ReferenceGenome` (GRCh38 default).

---

## 7. Open Questions / Decisions

1. Registry lists `DetectWholeGenomeDoubling(ploidy)` (scalar). The authoritative facets-suite/Bielski WGD definition is the major-CN≥2 / >50%-of-genome rule, which requires segments, not a scalar ploidy. The canonical method therefore takes segments. Decision recorded; registry method-signature note updated in the algorithm doc.
2. The WGD denominator is now the reference autosomal genome length (chromosome-size table) per facets-suite `autosomal_genome`. The caller selects `ReferenceGenome` (GRCh38 default) to match the coordinate system of its segments. The legacy supplied-segment-length behaviour is retained as `DetectWholeGenomeDoublingFromSuppliedLength` for callers whose segments already tile the genome.
</content>
