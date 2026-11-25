# Test Specification: ASSEMBLY-QC-001

**Test Unit ID:** ASSEMBLY-QC-001
**Area:** Assembly
**Algorithm:** Assembly local-quality windows and sequence-only suspicious-region screen
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect (B10 sweep U5)
**Last Updated:** 2026-11-25

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | URL | Accessed |
|---|--------|---------------|-----|----------|
| 1 | Biopython `Bio.SeqUtils.gc_fraction` / `GC_skew` (windowed slicing `seq[i:i+w]`, partial last window) 1.88 | 3 | local install | 2026-11-25 |
| 2 | QUAST `basic_stats.GC_content` (via ASSEMBLY-STATS-001, F10): (G+C)/(A+C+G+T) | 3 | docs/Evidence/ASSEMBLY-STATS-001-Evidence.md | – |
| 3 | lh3/sdust `sdust.c` (`sdust_core`, `shift_window`, `find_perfect`, `save_masked_regions`) = Morgulis et al. 2006 SDUST | 3 | https://raw.githubusercontent.com/lh3/sdust/master/sdust.c | 2026-11-25 |
| 4 | Linguistic complexity (Orlov & Potapov 2004; Troyanskaya 2002; Rosalind LING) via canonical `SequenceComplexity` | 1 | docs/Evidence/ASSEMBLY-REPEAT-001-Evidence.md | – |
| 5 | QUAST / GAGE misassembly definitions (alignment based, no sequence-only threshold) | 3 | docs/Evidence/ASSEMBLY-QC-001-Evidence.md | 2026-11-25 |

### 1.2 Key Evidence Points

1. Windowed GC is a Core primitive (`GcSkewCalculator.CalculateWindowedGcContent`); the trailing partial window follows Biopython `GC_skew` (new additive `includePartialWindow` overload, B01 R3 request).
2. GC denominator = A+C+G+T(+U) only: N and IUPAC codes are excluded (F10 / QUAST).
3. SDUST perfect intervals are the sourced definition of "low complexity" (DUST score Σ c_t(c_t−1)/2 / (ℓ−1), triplets, W = 64, level 20); the former `linguistic complexity < 0.3` cut-off has no source.
4. No published threshold defines a suspicious GC or N window: the cut-offs are parameters, the composite score is gone (reported values are the measured criterion values).

### 1.3 Documented Corner Cases

- A window without any A/C/G/T/U base has GC 0 by the canonical helper but no GC: it is not a GC outlier (all-N window is N-rich only).
- A truncated trailing window already covered by the previous window is not tested (avoids tiny-window noise); a sequence shorter than the window is one window.
- Pooled GC over all sequences, not per sequence; no valid base anywhere → no GC test.
- Regions are the union of flagged windows/intervals: End is the last flagged base (the old code ended at the next clean window start − 1, inside the last flagged overlapping window).

### 1.4 Known Failure Modes / Pitfalls

1. `windowSize = 1` ⇒ step `1/2 = 0` ⇒ infinite loop; `windowSize ≤ 0` ⇒ infinite loop / exception inside the iterator.
2. Global GC denominator counted every non-N symbol (IUPAC included) while window GC used A/C/G/T ⇒ inconsistent deviation.
3. Unsourced score (`|ΔGC|/0.15`, `1 − LC/0.3`, N fraction) mixed units and was presented as one number.
4. `minComplexity` on linguistic complexity (windows of 500 → LC ≈ 0.9 for random sequence) is not comparable to any published threshold.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateLocalQuality(seqs, windowSize=1000, stepSize=null)` | GenomeAssemblyAnalyzer | Delegating | GC: `GcSkewCalculator.CalculateWindowedGcContent(includePartialWindow)`; complexity: `SequenceComplexity.CalculateLinguisticComplexity` |
| `FindSuspiciousRegions(seqs, gcDeviation=0.15, dustThreshold=2.0, maxNFraction=0.1, windowSize=500)` | GenomeAssemblyAnalyzer | Documented screen | low complexity: `SequenceComplexity.FindLowComplexityIntervals` (SDUST) |
| `GcSkewCalculator.CalculateWindowedGcContent(seq, w, step, fraction, includePartialWindow)` | GcSkewCalculator (Analysis) | Additive overload | |
| `CalculateLengthDistribution` | GenomeAssemblyAnalyzer | re-checked (F13) | |
| MCP `local_quality` (+`stepSize`), `find_suspicious_regions` | ChromosomeTools | Wrapper | delegate |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Window k starts at k·step, length min(w, n−start); every base covered | Yes | Source 1 |
| INV-2 | GC = Core windowed GC; = Biopython gc_fraction on ACGT-only windows | Yes | Source 1, 2 |
| INV-3 | Low-complexity region = SDUST interval (exact), score = DUST score of it | Yes | Source 3 |
| INV-4 | Regions sorted, non-overlapping, non-adjacent per sequence; Reason order N, GC, LC | Yes | – |
| INV-5 | Complexity = canonical LC of the (upper-cased) window | Yes | Source 4 |

---

## 4. Test Cases (`GenomeAssemblyAnalyzer_QualityScreen_Tests`)

### 4.1 MUST

| ID | Test Case | Expected Outcome | Test |
|----|-----------|------------------|------|
| M1 | 1500-bp random, w 1000 | windows (0,1000),(500,1000),(1000,500); GC / LC = oracle | `CalculateLocalQuality_MatchesBiopythonGcAndLinguisticComplexityOracle` |
| M2 | N run, w 300, step 100 | 8 windows incl. truncated; N counts 0/0/50/50/50/0/0/0; LC with a = 5 | `CalculateLocalQuality_NRunAndExplicitStep_MatchOracle` |
| M3 | GC denominator | "ACGTRYSWNN" → 0.5 | `CalculateLocalQuality_GcDenominatorExcludesNAndIupac_AsQuast` |
| M4 | Delegation to Core | equal to Core windowed GC | `CalculateLocalQuality_GcIsTheCanonicalWindowedGcProfile` |
| M5 | Core `includePartialWindow` | starts 0,3,6,9; ends 3,6,9,9; GC 1,.25,.5,1 | `WindowedGcContent_IncludePartialWindow_*` |
| M6 | Suspicious oracle fixtures (poly-A, poly-A small, GC island, N run, multi-sequence) | regions + measured values = oracle | `FindSuspiciousRegions_*MatchesOracle` (5) |
| M7 | poly(A)600 | [0,599] "Low complexity", DUST 299 | `FindSuspiciousRegions_PolyA600_SingleDustIntervalWithExactScore` |
| M8 | Pooled GC | gc/at sequences → ΔGC 0.75 / 0.25 | `FindSuspiciousRegions_GcDeviationUsesPooledGcOverAllSequences` |
| M9 | Validation | exceptions | `CalculateLocalQuality_InvalidArguments_Throw`, `FindSuspiciousRegions_InvalidArguments_Throw`, `CalculateLocalQuality_WindowOfOne_Terminates…` |

### 4.2 SHOULD

| ID | Test Case | Test |
|----|-----------|------|
| S1 | all-N window not a GC outlier | `FindSuspiciousRegions_WindowWithoutAcgt_IsNotAGcOutlier` |
| S2 | covered truncated window skipped; short sequence one window | `…_TruncatedTrailingWindow…`, `…_SequenceShorterThanWindow…` |
| S3 | dustThreshold forwarded; reason order | `…_DustThreshold_IsPassedToSdust`, `…_ReasonsAreJoinedInFixedOrder` |
| S4 | CalculateLengthDistribution half-open bins | `CalculateLengthDistribution_HalfOpenBinsAndOverflow` |

### 4.3 MCP (`Seqeron.Mcp.Chromosome.Tests`)

`LocalQuality_StepSizeAndTruncatedLastWindow`, `LocalQuality_WindowOfOne_Terminates`, `FindSuspiciousRegions_PolyA_ReportsMeasuredDustScoreAndNoScore`, `FindSuspiciousRegions_NRun_IsFlaggedAndParametersAreForwarded`.

---

## 5. Audit of Existing Tests

| Test | Status | Action |
|------|--------|--------|
| `GenomeAssemblyAnalyzerTests.FindSuspiciousRegions_LowComplexity_Detected` | used `minComplexity` | now `dustThreshold` |
| MK4 `FindSuspiciousRegions_PolyA_ExactLowComplexityScore` | locked the unsourced `1 − LC/0.3` score | replaced by exact DUST score 299 |
| MCP `FindSuspiciousRegions_LowComplexityRun_IsFlagged` | param 0.3 | dust threshold 2.0 |

---

## 6. Deviations and Assumptions

| # | Item | Source behaviour | Seqeron behaviour |
|---|------|------------------|-------------------|
| 1 | Window / step (500, w/2) and N/GC cut-offs | none published | caller parameters; documented as a screen |
| 2 | Pooled GC as reference | none | assembly-wide QUAST GC |
| 3 | Misassembly detection | alignment based (QUAST/GAGE) | not a misassembly caller; out of scope without a reference |
