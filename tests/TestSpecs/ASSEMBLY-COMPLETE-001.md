# Test Specification: ASSEMBLY-COMPLETE-001

**Test Unit ID:** ASSEMBLY-COMPLETE-001
**Area:** Assembly
**Algorithm:** Assembly completeness — BUSCO-like marker-gene classification and k-mer spectrum evaluation (Merqury completeness / QV, GenomeScope error rate / genome size)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect (B10 sweep U2)
**Last Updated:** 2026-11-25

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Simão et al. 2015 BUSCO; Manni et al. 2021 BUSCO v5 — `busco_tools/hmmer.py` (`_sort_matches`, `_get_busco_percentages`, `record_results`) | 3 (reference implementation) | 10.1093/bioinformatics/btv351; https://gitlab.com/ezlab/busco/-/raw/master/src/busco/busco_tools/hmmer.py | 2026-11-25 |
| 2 | Huang & Li 2023 compleasm — `compleasm.py` v0.2.6 (`pass_tids`, `record_1st_gene_label`, `get_region_clusters`) and main | 3 | 10.1093/bioinformatics/btad595; https://raw.githubusercontent.com/huangnengCSU/compleasm/v0.2.6/compleasm.py | 2026-11-25 |
| 3 | Rhie et al. 2020 Merqury — `eval/qv.sh`, `eval/spectra-cn.sh`, `build/filt.sh`, `HistToPloidyDepth.java` | 3 | 10.1186/s13059-020-02134-9; https://raw.githubusercontent.com/marbl/merqury/master/ | 2026-11-25 |
| 4 | Vurture et al. 2017 GenomeScope — `genomescope.R` (`error_rate` uniroot, `total_len`) | 3 | 10.1093/bioinformatics/btx153; https://raw.githubusercontent.com/schatzlab/genomescope/master/genomescope.R | 2026-11-25 |
| 5 | BLAST+ defaults (megablast word_size 28; blastn reward 2 / penalty −3 / gap 5,2 = `SequenceAligner.BlastDna`) | 2 | BLAST+ User Manual | 2026-11-25 |
| 6 | Biopython 1.88 `PairwiseAligner` (local), scipy 1.18.1 `brentq` | 3 | pip-installed | 2026-11-25 |

### 1.2 Key Evidence Points

1. BUSCO classes Complete (= Single-copy + Duplicated), Fragmented, Missing; every percentage is over n; summary `C:x%[S:x%,D:x%],F:x%,M:x%,n:N` with `round(·, 1)` (source 1).
2. compleasm (aligned-fraction mode): valid hit = mapped rate ≥ min_length_percent (0.6) OR identity ≥ min_identity; complete = mapped rate ≥ min_complete (0.9); several complete hits → Duplicated on different contigs or > 1 overlap cluster (half-open, abutting merges) (source 2).
3. Merqury: solid = count > filt, filt = first histogram valley of `kmerHistToPloidyDepth`; completeness = |asm ∩ solid| / |solid|; QV = −10·log10(1 − (1 − ASM_ONLY/TOTAL)^(1/k)) with ASM_ONLY counted against all read k-mers (source 3).
4. GenomeScope: read error rate = root of Σ C(k,i)x^i(1−x)^(k−i)·T = E, i.e. 1 − (1 − E/T)^(1/k); genome size = (T − E)/(2·kcov) (source 4).

### 1.3 Documented Corner Cases

- Marker on the reverse strand → '-' hit, forward marker coordinates.
- Abutting complete copies on one contig → Single (compleasm cluster merge); spaced copies → Duplicated.
- Marker shorter than 28 bp → seeded by its whole length (exact occurrence required).
- Histogram with a zero-slope plateau before the rise → no valley, filt 0 (Merqury); histogram rising from its first line → filt = second depth − 1.
- No assembly-only k-mer → QV +∞ (MCP null); assembly shorter than k → QV/E NaN (Merqury "asm has no kmers").
- No assembly → completeness NaN (MCP null).

### 1.4 Known Failure Modes / Pitfalls

1. Forward-strand-only marker search (old) misses reverse-strand genes.
2. D% relative to C instead of n (old).
3. Per-distinct-k-mer singleton fraction reported as a per-base error rate (old).
4. Including error k-mers in the genome-size numerator (old).
5. Percent evaluated as count·100/n instead of BUSCO's count/n·100 changes `round(·,1)` at ties (199/2000).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `AssessCompleteness(assembly, markerGenes, identityThreshold = 0.9, coverageThreshold = 0.9, minLengthFraction = 0.6)` | GenomeAssemblyAnalyzer | Canonical | compleasm classification, BUSCO percentages |
| `FindMarkerHits(marker, assembly)` | GenomeAssemblyAnalyzer | Canonical | Seeded SW-Gotoh (`SequenceAligner.LocalAlignAffine`, `BlastDna`) |
| `CompletenessResult.BuscoNotation` / `SingleCopyPercent` / `FragmentedPercent` / `MissingPercent` | GenomeAssemblyAnalyzer | Derived | BUSCO `record_results` |
| `EstimateCompletenessFromKmers(spectrum, assembly?, expectedCoverage = 0, canonical = true)` | GenomeAssemblyAnalyzer | Canonical | Merqury + GenomeScope; assembly k-mers via `KmerAnalyzer.CountKmers(AcgtOnly, Canonical)` |
| `EstimateCompletenessFromKmers(spectrum, expectedCoverage = 0)` | GenomeAssemblyAnalyzer | Legacy tuple | Delegates; completeness NaN |
| MCP `assess_completeness`, `estimate_completeness_from_kmers` | ChromosomeTools | Wrapper | Delegate |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | S + D + F + M = n; C = S + D; all percentages = class / n · 100 | Yes | Source 1 |
| INV-2 | Hit coverage = (MarkerEnd − MarkerStart + 1)/m ∈ (0, 1]; identity ∈ [0, 1]; 0 ≤ Start ≤ End < contig length | Yes | Source 2 |
| INV-3 | 0 ≤ completeness ≤ 1; AssemblySolidKmers ≤ SolidKmers; ASM_ONLY ≤ TOTAL | Yes | Source 3 |
| INV-4 | ErrorRate = 1 − (1 − E/T)^(1/k) equals GenomeScope's uniroot (|Δ| ≤ 2.3e-16) | Yes | Source 4, oracle |
| INV-5 | Solid k-mers have count > filt, error k-mers count ≤ filt | Yes | Source 3 |

---

## 4. Test Cases

### 4.1 MUST Tests

| ID | Test Case | Expected Outcome | Test |
|----|-----------|------------------|------|
| M1 | 13 oracle scenarios (strands, duplicates, tandem, fragments, indel, substitutions, lower case/N) — hits and class | Evidence table | `GenomeAssemblyAnalyzer_Completeness_Tests.FindMarkerHits_AndClassification_MatchOracle` |
| M2 | Reverse-strand coordinates | c1 − 100–202, marker 97–199, 102/103, 103/300 | `FindMarkerHits_ReverseStrand_ReportsForwardCoordinates` |
| M3 | Seed limit (substitution every 20 bp) | Missing | `AssessCompleteness_SeedRequirement_DivergedMarkerWithoutExactWordIsMissing` |
| M4 | Short marker (< 28 bp) | exact hit, Single | `AssessCompleteness_ShortMarker_IsSeededByItsWholeLength` |
| M5 | compleasm pass filter / inclusive thresholds | Fragmented / Missing / Fragmented / Single | `AssessCompleteness_Thresholds_FollowCompleasmPassFilter` |
| M6 | BUSCO percentages over n (multi-marker) | S1 D1 F1 M1, D% 25, notation | `AssessCompleteness_MultiMarker_BuscoPercentagesOverN`; MCP `AssessCompleteness_DuplicatedAndMissing_BuscoPercentagesOverN` |
| M7 | BUSCO notation rounding | 4 oracle strings | `BuscoNotation_MatchesBuscoOneLineSummary` (4) |
| M8 | Validation (both methods) | ArgumentNull / ArgumentOutOfRange / Argument | `AssessCompleteness_InvalidArguments_Throw`, `EstimateCompletenessFromKmers_InvalidArguments_Throw`; MCP `*_Schema_ValidatesCorrectly` |
| M9 | Merqury filt / GenomeScope error & size / Merqury completeness & QV | valley table row 1 | `EstimateCompletenessFromKmers_MerqurySolidThreshold_GenomeScopeErrorAndSize_MerquryCompletenessAndQv`; MCP `EstimateCompletenessFromKmers_WithAssembly_MerquryCompletenessAndQv` |
| M10 | canonical = false | 1/12, QV 4.3175 | `EstimateCompletenessFromKmers_NonCanonical_CountsAssemblyForwardOnly` |
| M11 | expectedCoverage override, no assembly | size 15, completeness NaN | `EstimateCompletenessFromKmers_ExpectedCoverage_OverridesPeak`; MCP `EstimateCompletenessFromKmers_SpectrumOnly_ErrorRateAndSize_NoCompleteness` |
| M12 | Merqury quirks (plateau, rising start) | filt 0 / 2 | `EstimateCompletenessFromKmers_MerquryQuirks_PlateauAndRisingStart` |
| M13 | QV +∞ / NaN | | `EstimateCompletenessFromKmers_NoAssemblyOnlyKmers_QvIsInfinite` |
| M14 | Random cross-check | oracle 400/400 (both); naive reference 400 + planted copies 150 | Fuzz `AssemblyCompletenessFuzzTests` (2) |

### 4.2 SHOULD Tests

| ID | Test Case | Expected Outcome | Test |
|----|-----------|------------------|------|
| S1 | Legacy tuple overload values | {4,4,2,1}: error 0, size 3; 1:100/30:1000 6-mers: error 5.5448e-4, size 1000 | `GenomeAssemblyAnalyzer_MutationKillers2_Tests.EstimateCompletenessFromKmers_PeakErrorAndSize`, `GenomeAssemblyAnalyzerTests.EstimateCompletenessFromKmers_ValidSpectrum_ReturnsEstimates` |
| S2 | One-line histogram, expected coverage half-even | size 2; round(5.5) = 6 | `GenomeAssemblyAnalyzer_MutationKillers4_Tests.EstimateCompletenessFromKmers_AllSingletons_AndExpectedCoverage` |
| S3 | Periodic self-RC marker stays single-copy | Single | MCP `AssessCompleteness_MarkerPresent_ReportsComplete` |

### 4.3 COULD Tests

| ID | Test Case | Expected Outcome | Test |
|----|-----------|------------------|------|
| C1 | Empty marker list | n 0, notation all 0.0 | `AssessCompleteness_InvalidArguments_Throw` |

---

## 5. Audit of Existing Tests

| Test | Status | Action |
|------|--------|--------|
| `GenomeAssemblyAnalyzer_MutationKillers2_Tests.EstimateCompletenessFromKmers_PeakErrorAndSize` | Locked singleton-ratio error 0.25, size 2, completeness 1 (unsourced) | Corrected to oracle: error 0, size 3, completeness NaN |
| `GenomeAssemblyAnalyzer_MutationKillers4_Tests.EstimateCompletenessFromKmers_AllSingletons_AndExpectedCoverage` | Locked (0, 1.0, 0) and size 5 (integer division) | Corrected: (NaN, 0, 2), size round(5.5) = 6 |
| `GenomeAssemblyAnalyzerTests.EstimateCompletenessFromKmers_ValidSpectrum_ReturnsEstimates` | Non-DNA k-mers of mixed length; weak bounds | Real 6-mers; exact oracle values |
| `GenomeAssemblyAnalyzer_MutationKillers2_Tests.AssessCompleteness_CompleteAndMissingMarkers` | "Absent" G80 is found on the reverse strand (C40) | Absent marker (GT)40 |
| `GenomeAssemblyAnalyzer_MutationKillers3/4/6` completeness tests, `GenomeAssemblyAnalyzerTests.AssessCompleteness_*` | Consistent with sources (D% equal for n = C) | Kept |
| MCP `EstimateCompletenessFromKmersTests` | Mixed-length pseudo k-mers, completeness 1 without assembly | Rewritten with oracle values |

## R1-WP1 addendum (2026-11-25): ASCII fold, span unit lookup, `IsValidDna` (A15-7, A11-6)

`GenomeAssemblyAnalyzer`: five `ToUpperInvariant` folds -> `SequenceCase.ToUpperAscii`; `FindRepetitiveRegions` unit counts via `KmerEnumerator` + span lookup (new private `UnitCounts`); private `IsAcgt` deleted -> `IsValidDna` (identical on the upper-cased input). New tests in `Unit/Chromosome/GenomeAssemblyAnalyzer_AsciiFold_Tests.cs`: `FindRepetitiveRegions_NonAsciiLetter_ActsAsNonAcgtSymbol` (U+017F, U+0130 vs N), `FindRepetitiveRegions_LowerCase_EqualsUpperCase`, `EstimateCompletenessFromKmers_LongSAndS_AreDistinctKmers`, `EstimateCompletenessFromKmers_Canonical_NonAsciiKmer_RejectedAsNotAcgt`, `FindMarkerHits_NonAsciiMarkerBase_FullLengthHitFound`, `FindMarkerHits_LowerCaseMarkerAndAssembly_FullIdentity`, `FindSyntenicBlocks_NonAsciiQueryLetter_ActsAsNonAcgt`, `CalculateLocalQuality_NonAsciiLetter_ExcludedFromGc`.

