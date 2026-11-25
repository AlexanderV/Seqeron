# Test Specification: ASSEMBLY-REPEAT-001

**Test Unit ID:** ASSEMBLY-REPEAT-001
**Area:** Assembly
**Algorithm:** Assembly repeat analysis (k-mer-count repeat regions, tandem repeats, repeat content, window complexity)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect (B10 sweep U3)
**Last Updated:** 2026-11-25

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | URL | Accessed |
|---|--------|---------------|-----|----------|
| 1 | NCBI WindowMasker (Morgulis et al. 2006): `seq_masker.cpp` (`DoMask`), `seq_masker_window.cpp`, `seq_masker_score_mean.cpp`, `win_mask_gen_counts.cpp` | 3 (reference implementation) | https://raw.githubusercontent.com/ncbi/ncbi-cxx-toolkit-public/master/src/algo/winmask/ | 2026-11-25 |
| 2 | RepeatMasker `ProcessRepeats` (`generateTableOutput`, aggregate statistics) | 3 | https://raw.githubusercontent.com/rmhubley/RepeatMasker/master/ProcessRepeats | 2026-11-25 |
| 3 | Tandem Repeats Finder (Benson 1999, NAR 27:573; TRF 4.10.0 README) via REP-APPROX-001 | 1 / 3 | docs/Evidence/REP-APPROX-001-Evidence.md | (earlier, project) |
| 4 | Linguistic complexity (Troyanskaya 2002; Orlov & Potapov 2004; Rosalind LING) via the canonical `SequenceComplexity` | 1 | doi:10.1093/bioinformatics/18.5.679 | – |

### 1.2 Key Evidence Points

1. WindowMasker: counts are pooled over the genome, both strands (canonical unit; self-complementary unit counted twice); only all-ACGT windows are scored; window score = integer mean of its unit counts; `DoMask` branches on T_extend / T_threshold build inclusive intervals (source 1).
2. RepeatMasker "bases masked" counts each base once; overlapping annotations contribute only the uncovered bases to their class; class lengths add up to the total (source 2).
3. TRF reports an array by alignment score (match 2, mismatch −7, indel −7, Minscore 50), with consensus pattern, fractional copy number and percent matches between adjacent copies (source 3).
4. Window complexity = canonical linguistic complexity, m = 10, N a fifth symbol (source 4).

### 1.3 Documented Corner Cases

- A window that contains a non-ACGT symbol (N, IUPAC, `-`) is skipped; units containing one are not counted (Jellyfish / WindowMasker rule).
- Palindromic even-k unit: both-strand count is 2 × forward count (`win_mask_gen_counts.cpp`).
- TRF: N never matches; a short N-run inside an array is bridged when the array's score outweighs the −7 per N column.
- Repeat content: contained annotations add 0 bases; the same coordinates on different sequences are different bases.

### 1.4 Known Failure Modes / Pitfalls

1. Old window heuristic (step `windowSize/10`) never advanced for `windowSize` < 10 (infinite loop).
2. 80 %-identity-to-first-copy tandem rule: no indels, integer copies, short arrays accepted that TRF rejects.
3. Summing overlapping annotations gave repeat content > 100 %.
4. Private k-mer complexity ≠ published linguistic complexity.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindRepetitiveRegions(seqs, k=15, minCopies=3, windowSize=100, bothStrands=true, extendCopies=null)` | GenomeAssemblyAnalyzer | Canonical (WindowMasker model) | counting via `KmerAnalyzer.CountKmers` (`AcgtOnly`) |
| `FindTandemRepeats(seqs, minUnit=2, maxUnit=50, minCopies=3, minScore=50)` | GenomeAssemblyAnalyzer | Delegating | `RepeatFinder.FindApproximateTandemRepeats` (TRF) |
| `CalculateRepeatContent(repeats, genomeLength)` | GenomeAssemblyAnalyzer | Canonical (RepeatMasker accounting) | |
| `CalculateLocalQuality(...).Complexity` | GenomeAssemblyAnalyzer | Delegating | `SequenceComplexity.CalculateLinguisticComplexity` |
| MCP `find_repetitive_regions`, `find_tandem_repeats`, `repeat_content` | ChromosomeTools | Wrapper | delegate |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Regions are inclusive [Start, End] of triggering windows, never overlap, ordered per sequence | Yes | Source 1 |
| INV-2 | Region ⇒ at least one window with mean unit count ≥ T_threshold | Yes | Source 1 |
| INV-3 | Total = distinct annotated bases; Σ class lengths = Total | Yes | Source 2 |
| INV-4 | FindTandemRepeats rows = TRF rows filtered by period / copies | Yes | Source 3 |
| INV-5 | Complexity = canonical LC of the window | Yes | Source 4 |

---

## 4. Test Cases

### 4.1 MUST Tests (`GenomeAssemblyAnalyzer_RepeatAnalysis_Tests`)

| ID | Test Case | Expected Outcome | Test |
|----|-----------|------------------|------|
| M1 | WindowMasker oracle (14 fixtures of 3000/3000) | intervals + copies identical | `FindRepetitiveRegions_MatchesWindowMaskerOracle` |
| M2 | (CAG)×200 | [0,599], copies 199 | `FindRepetitiveRegions_TandemArray_OneRegionWithMaxUnitCount`, `GenomeAssemblyAnalyzer_MutationKillers5_Tests.FindRepetitiveRegions_TandemSpansWholeSequence` |
| M3 | Both strands, pooled over sequences | forward-only none; both strands 3 regions | `FindRepetitiveRegions_CountsArePooledOverSequencesAndBothStrands` |
| M4 | IUPAC R excludes a unit | no region | `FindRepetitiveRegions_KmersWithAnyNonAcgtSymbolAreNotCounted` |
| M5 | TRF perfect / substituted / flanked arrays | (AT)×30 score 120; (ATG)×7·ATC·(ATG)×7 15 copies, 40/42; flanked → [9,36] GCA 28/3 | `FindTandemRepeats_PerfectArray_TrfValues`, `…_SubstitutedCopy_IsAlignedByTrf`, `…_SingleSubstitution_ExactTrfPurity`, `…_FlankedArray_FractionalCopyNumberAndPhase`, `…_PerfectDinucleotideRepeat` |
| M6 | Equals RepeatFinder TRF rows | identical | `FindTandemRepeats_MatchesRepeatFinderTrfRows` |
| M7 | Repeat content once-per-base | overlap/contained → 150, LINE 100 / SINE 50 | `CalculateRepeatContent_OverlappingAnnotations_AreCountedOnce`, per-base oracle `…_RandomAnnotations_MatchPerBaseOracle` (2000) |
| M8 | Window complexity canonical | = `SequenceComplexity`; "ACGTACGT" 26/32; polyA 10/3301; "ACNGT" 1 | `CalculateLocalQuality_Complexity_*`, `CalculateLocalQuality_ExactWindowMetrics`, `FindSuspiciousRegions_PolyA_ExactLowComplexityScore` |
| M9 | Validation | eager exceptions | `FindRepetitiveRegions_InvalidArguments_ThrowEagerly`, `FindTandemRepeats_InvalidArguments_ThrowEagerly`, `CalculateRepeatContent_InvalidInput_Throws` |

### 4.2 SHOULD Tests

| ID | Test Case | Expected Outcome | Test |
|----|-----------|------------------|------|
| S1 | minScore / minCopies / period filters | filtered | `FindTandemRepeats_BelowMinScore_NotReported`, `…_MinCopiesAndPeriodRange_Filter` |
| S2 | N never matches (bridged vs split) | 1 repeat vs [0,44],[85,129] | `FindTandemRepeats_NNeverMatches` |
| S3 | Case-insensitive; short/empty input; unique sequence | same / none | `FindRepetitiveRegions_IsCaseInsensitive`, `…_ShortOrEmptySequences_YieldNothing`, `…_UniqueSequence_NoRegion` |
| S4 | Separate sequences / adjacency / order | counted separately | `CalculateRepeatContent_SameCoordinatesOnDifferentSequences_AreSeparate`, `…_AdjacentAnnotations_AddUp_…` |

### 4.3 MCP (`Seqeron.Mcp.Chromosome.Tests`)

`FindRepetitiveRegions_Binding_InvokesSuccessfully` (ACGTACGT×100 → [0,799], 398), `FindRepetitiveRegions_StrandAndExtendOptions_AreForwarded`, `FindTandemRepeats_SubstitutedArray_ReportsTrfFractionalPurityAndMinScore`, `FindTandemRepeats_Schema_ValidatesCorrectly` (minScore), `RepeatContent_OverlappingAnnotations_CountEachBaseOnce`.

---

## 5. Audit of Existing Tests

| Test | Status | Action |
|------|--------|--------|
| MK2 `FindTandemRepeats_PerfectDinucleotideRepeat` | locked an 8-bp (AT)×4 repeat and int copies | now (AT)×30 (TRF) |
| MK3 `FindTandemRepeats_ImperfectCopyStillCountedAt80Percent` | locked the 80 % rule | replaced by TRF values (renamed `…_SubstitutedCopy_IsAlignedByTrf`) |
| MK4 `FindTandemRepeats_ImperfectCopy_ExactPurity` | locked 14/15 of the 80 % model | TRF 93/95 (renamed `…_SingleSubstitution_ExactTrfPurity`) |
| MK5 `FindRepetitiveRegions_TandemSpansWholeSequence` | locked End 593 ("truncated final windows") | [0,599], 199 |
| MK2 `CalculateLocalQuality_ExactWindowMetrics`, MK4 `FindSuspiciousRegions_PolyA_ExactLowComplexityScore` | locked the private complexity | canonical LC values |
| MCP `FindRepetitiveRegions_Binding_InvokesSuccessfully` | End 789 / 199 | 799 / 398 |

---

## 6. Deviations and Assumptions

| # | Item | Source behaviour | Seqeron behaviour |
|---|------|------------------|-------------------|
| 1 | T_threshold / T_extend | percentiles of the count distribution | counts (minCopies / extendCopies) |
| 2 | WindowMasker merge pass, discontiguous units, T_low/T_high | optional | not implemented |
| 3 | `minCopies` of FindTandemRepeats | TRF reports ≥ 1.9 copies | extra library filter (default 3) |
| 4 | Overlap between classes | RepeatMasker resolved annotation order | first in (start, input order) |
