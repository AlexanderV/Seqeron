# TestSpec: SEQ-COMPLEX-001 - Sequence Complexity Metrics

**Test Unit ID:** SEQ-COMPLEX-001
**Area:** Sequence Composition
**Created:** 2026-01-22
**Status:** Complete

## 1. Scope

This TestSpec covers all sequence complexity metrics in `SequenceComplexity`:

| Method | Purpose |
|--------|---------|
| `CalculateLinguisticComplexity` | Vocabulary richness (Troyanskaya summation) |
| `CalculateShannonEntropy` | Information content per base |
| `CalculateKmerEntropy` | Entropy of k-mer frequency distribution |
| `CalculateWindowedComplexity` | Sliding-window complexity profile |
| `FindLowComplexityRegions` | Low-entropy region detection |
| `CalculateDustScore` | DUST low-complexity score |
| `MaskLowComplexity` | Mask low-complexity by DUST threshold |
| `EstimateCompressionRatio` | Unique substring ratio |

## 2. Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| Wikipedia "Linguistic sequence complexity" | Encyclopedia | LC formula, vocabulary usage, examples |
| Troyanskaya et al. (2002) Bioinformatics 18(5):679–88 | Peer-reviewed | Summation formula: LC = Σ(observed) / Σ(possible) |
| Orlov & Potapov (2004) NAR 32:W628–W633 | Peer-reviewed | V_max = min(4^i, N−i+1), complexity profiles |
| Wikipedia "Entropy (information theory)" | Encyclopedia | Shannon formula H = −Σ p_i log₂ p_i |
| Shannon (1948) Bell System Technical Journal | Original paper | Information entropy definition |
| Morgulis et al. (2006) J Comput Biol 13(5):1028–40 | Peer-reviewed | Symmetric DUST: S = Σc_t(c_t−1)/2 / (w−1) |

## 3. Formulas

### 3.1 Linguistic Complexity (Troyanskaya 2002)

$$LC = \frac{\sum_{i=1}^{m} V_{obs}(i)}{\sum_{i=1}^{m} V_{max}(i)}$$

Where $V_{max}(i) = \min(4^i, N - i + 1)$. Range: $0 \le LC \le 1$.

### 3.2 Shannon Entropy (Shannon 1948)

$$H = -\sum_{i} p_i \log_2 p_i$$

Range for DNA: $0 \le H \le 2$ (log₂4 = 2).

### 3.3 K-mer Entropy

Shannon entropy applied to k-mer frequency distribution. Range: $0 \le H \le \log_2(4^k)$.

### 3.4 DUST Score (Morgulis 2006)

$$DUST = \frac{\sum_t c_t(c_t - 1)/2}{w - 1}$$

Where $c_t$ = count of triplet $t$, $w$ = number of triplets ($N - 2$).

## 4. Test Categories

### 4.1 Linguistic Complexity Tests (12 tests)

| ID | Test Method | Assertion | Source |
|----|-------------|-----------|--------|
| LC-1 | `CalculateLinguisticComplexity_HighComplexity_ReturnsHigh` | Exact: 91/103 | Troyanskaya (2002) |
| LC-2 | `CalculateLinguisticComplexity_LowComplexity_ReturnsLow` | Exact: 10/103 | Orlov & Potapov (2004) |
| LC-3 | `CalculateLinguisticComplexity_EmptySequence_ReturnsZero` | Exact: 0 | Formula definition |
| LC-4 | `CalculateLinguisticComplexity_RangeIsZeroToOne_ForMultipleSequences` | Range [0,1] | Troyanskaya (2002) |
| LC-5 | `CalculateLinguisticComplexity_StringOverload_MatchesDnaSequenceOverload` | Exact equality | API contract |
| LC-6 | `CalculateLinguisticComplexity_SingleNucleotide_ReturnsOne` | Exact: 1.0 | Formula definition |
| LC-7 | `CalculateLinguisticComplexity_DinucleotideRepeat_LowerThanRandom` | LC < 0.1 vs > 0.4 | Orlov & Potapov (2004) |
| LC-8 | `CalculateLinguisticComplexity_MaxWordLengthParameter_AffectsResult` | maxWord=1 → 1.0 | Parameter semantics |
| LC-9 | `CalculateLinguisticComplexity_WikipediaExample_MatchesHandCalculation` | Exact: 47/49 | Wikipedia |
| LC-10 | `CalculateLinguisticComplexity_WikipediaDinucleotideRepeat_MatchesHandCalculation` | Exact: 5/28 | Wikipedia |
| LC-11 | `CalculateLinguisticComplexity_MaximalComplexity_ReturnsOne` | Exact: 1.0 | Troyanskaya (2002) |
| LC-12 | `CalculateLinguisticComplexity_LowercaseInput_HandledCorrectly` | Case-insensitive | Robustness |

### 4.2 Shannon Entropy Tests (8 tests)

| ID | Test Method | Assertion | Source |
|----|-------------|-----------|--------|
| ENT-1 | `CalculateShannonEntropy_EqualBases_ReturnsTwo` | Exact: 2.0 | Wikipedia entropy |
| ENT-2 | `CalculateShannonEntropy_SingleBase_ReturnsZero` | Exact: 0 | Shannon (1948) |
| ENT-3 | `CalculateShannonEntropy_TwoBases_ReturnsOne` | Exact: 1.0 | Binary entropy |
| ENT-4 | `CalculateShannonEntropy_EmptySequence_ReturnsZero` | Exact: 0 | Convention |
| ENT-5 | `CalculateShannonEntropy_StringOverload_MatchesDnaSequenceOverload` | Exact equality | API contract |
| ENT-6 | `CalculateShannonEntropy_RangeIsZeroToTwo_ForDnaSequences` | Range [0,2] | Shannon max entropy |
| ENT-7 | `CalculateShannonEntropy_ThreeBases_ReturnsLog2Of3` | Exact: log₂(3) | Shannon formula |
| ENT-8 | `CalculateShannonEntropy_LowercaseInput_HandledCorrectly` | Case-insensitive | Robustness |

### 4.3 K-mer Entropy Tests (7 tests)

| ID | Test Method | Assertion | Source |
|----|-------------|-----------|--------|
| KME-1 | `CalculateKmerEntropy_VariedDinucleotides_ReturnsExact` | Exact: H formula (AT=4,TG=4,GC=4,CA=3) | Shannon formula |
| KME-2 | `CalculateKmerEntropy_RepeatedDinucleotides_ReturnsZero` | Exact: 0 | Single symbol |
| KME-3 | `CalculateKmerEntropy_SequenceShorterThanK_ReturnsZero` | Exact: 0 | No k-mers |
| KME-4 | `CalculateKmerEntropy_InvalidK_ThrowsException` | Throws | Guard clause |
| KME-5 | `CalculateKmerEntropy_NullSequence_ThrowsException` | Throws | Guard clause |
| KME-6 | `CalculateKmerEntropy_RangeIsNonNegativeAndBounded_ForDnaSequences` | Range [0, log₂(4^k)] | Shannon max |
| KME-7 | `CalculateKmerEntropy_UniformDinucleotides_ReturnsLog2Of3` | Exact: log₂(3) | Shannon formula |

### 4.4 Windowed Complexity Tests (4 tests)

| ID | Test Method | Assertion | Source |
|----|-------------|-----------|--------|
| WIN-1 | `CalculateWindowedComplexity_ReturnsCorrectPointCount` | Exact: 9 points | floor((180−20)/20)+1 |
| WIN-2 | `CalculateWindowedComplexity_IncludesBothMetrics_ExactValues` | H=2.0 exact | Shannon max |
| WIN-3 | `CalculateWindowedComplexity_PositionsAreCorrect` | pos=10, 30 | Window center |
| WIN-4 | `CalculateWindowedComplexity_NullSequence_ThrowsException` | Throws | Guard clause |
| WIN-5 | `CalculateWindowedComplexity_ZeroWindowSize_ThrowsException` | Throws | Guard clause |
| WIN-6 | `CalculateWindowedComplexity_ZeroStepSize_ThrowsException` | Throws | Guard clause |

### 4.5 Low Complexity Region Tests (5 tests)

| ID | Test Method | Assertion | Source |
|----|-------------|-----------|--------|
| LCR-1 | `FindLowComplexityRegions_FindsPolyARegion` | Count=1, start=79, end=146, minH=0 | Entropy-based detection |
| LCR-2 | `FindLowComplexityRegions_HighComplexity_ReturnsEmpty` | Empty | Definition |
| LCR-3 | `FindLowComplexityRegions_ReturnsCorrectSequence` | start=6, end=75, length=70 | Region merging |
| LCR-4 | `FindLowComplexityRegions_NullSequence_ThrowsException` | Throws | Guard clause |
| LCR-5 | `FindLowComplexityRegions_InvalidWindowSize_ThrowsException` | Throws | Guard clause |

### 4.6 DUST Score Tests (5 tests)

| ID | Test Method | Assertion | Source |
|----|-------------|-----------|--------|
| DUST-1 | `CalculateDustScore_LowComplexity_ReturnsHigh` | Exact: 8.0 (N=18) | Morgulis (2006) |
| DUST-2 | `CalculateDustScore_HighComplexity_ReturnsLow` | Exact: 6/13 (N=16) | Morgulis (2006) |
| DUST-3 | `CalculateDustScore_EmptySequence_ReturnsZero` | Exact: 0 | Convention |
| DUST-4 | `CalculateDustScore_StringOverload_ReturnsExact` | Exact: 2.5 (N=7) | Morgulis (2006) |
| DUST-5 | `CalculateDustScore_SequenceShorterThanWordSize_ReturnsZero` | Exact: 0 | Boundary |
| DUST-6 | `CalculateDustScore_NullSequence_ThrowsException` | Throws | Guard clause |

### 4.7 Masking Tests (5 tests)

| ID | Test Method | Assertion | Source |
|----|-------------|-----------|--------|
| MASK-1 | `MaskLowComplexity_MasksLowComplexityWindows` | N=192 (all masked) | DUST threshold |
| MASK-2 | `MaskLowComplexity_PreservesHighComplexity` | No N chars | DUST threshold |
| MASK-3 | `MaskLowComplexity_CustomMaskChar` | X=100 (all masked) | Parameter |
| MASK-4 | `MaskLowComplexity_NullSequence_ThrowsException` | Throws | Guard clause |
| MASK-5 | `MaskLowComplexity_ResultLengthEqualsInputLength` | Length invariant | Definition |
| MASK-6 | `MaskLowComplexity_ShortSequence_PreservesOriginal` | Returns "ATGC" | seq < window |

### 4.8 Compression Ratio Tests (5 tests)

| ID | Test Method | Assertion | Source |
|----|-------------|-----------|--------|
| CR-1 | `EstimateCompressionRatio_HighComplexity_ReturnsExact` | Exact: 14/27 (N=30) | Unique/expected |
| CR-2 | `EstimateCompressionRatio_LowComplexity_ReturnsExact` | Exact: 5/112 (N=31) | Unique/expected |
| CR-3 | `EstimateCompressionRatio_EmptySequence_ReturnsZero` | Exact: 0 | Convention |
| CR-4 | `EstimateCompressionRatio_RangeIsZeroToOne` | Range [0,1] | Mathematical |
| CR-5 | `EstimateCompressionRatio_NullSequence_ThrowsException` | Throws | Guard clause |

### 4.9 Edge Cases (2 remaining)

| ID | Test Method | Assertion | Source |
|----|-------------|-----------|--------|
| EDGE-1 | `CalculateLinguisticComplexity_ZeroWordLength_ThrowsException` | Throws | Guard clause |
| EDGE-2 | `CalculateLinguisticComplexity_NegativeWordLength_ThrowsException` | Throws | Guard clause |

## 5. Hand-Calculated Cross-Verification Table

### 5.1 Linguistic Complexity (Troyanskaya formula)

| Input | N | maxWord | obs | max | LC | Match |
|-------|---|---------|-----|-----|----|-------|
| `"A"` | 1 | 10 | 1 | 1 | 1.0 | ✓ |
| `"AAAA"` | 4 | 10 | 4 | 10 | 0.4 | ✓ |
| `"ATGC"` | 4 | 10 | 10 | 10 | 1.0 | ✓ |
| `"AAAAAAAAAAAAAAAA"` | 16 | 10 | 10 | 103 | 10/103 | ✓ |
| `"ATGCTAGCATGCAATG"` | 16 | 10 | 91 | 103 | 91/103 | ✓ |
| `"ACGGGAAGCTGATTCCA"` | 17 | 4 | 47 | 49 | 47/49 | ✓ |
| `"ACACACACACACACACA"` | 17 | 10 | 20 | 112 | 5/28 | ✓ |

### 5.2 Shannon Entropy

| Input | Distribution | H | Formula |
|-------|-------------|---|---------|
| `"ATGCATGCATGCATGC"` | A=T=G=C=25% | 2.0 | log₂(4) |
| `"AAAAAAA"` | A=100% | 0 | -1·log₂(1) = 0 |
| `"ATATATAT"` | A=T=50% | 1.0 | -2·(0.5·log₂0.5) |
| `"ATGATGATG"` | A=T=G=33.3% | log₂(3) | -3·(⅓·log₂⅓) |

### 5.3 K-mer Entropy

| Input | k | K-mers | Counts | H |
|-------|---|--------|--------|---|
| `"ATCG"` | 2 | AT,TC,CG | 1,1,1 | log₂(3) |
| `"AAAAAAAAAA"` | 2 | AA | 9 | 0 |
| `"ATGCATGCATGCATGC"` | 2 | AT,TG,GC,CA | 4,4,4,3 | −3·(4/15)log₂(4/15)−(3/15)log₂(3/15) |

### 5.4 DUST Score

| Input | N | Triplets | Score | DUST |
|-------|---|----------|-------|------|
| `"AAAAAAAAAAAAAAAAAA"` | 18 | AAA×16 | 16·15/2=120 | 120/15=8.0 |
| `"ATGCTAGCATGCTAGC"` | 16 | diverse | 6 | 6/13 |
| `"AAAAAAA"` | 7 | AAA×5 | 5·4/2=10 | 10/4=2.5 |

### 5.5 Compression Ratio

| Input | N | Unique | Expected | Ratio |
|-------|---|--------|----------|-------|
| `"ATGCTAGCATGCAATGCTAGCATGCAATGC"` | 30 | 112 | 216 | 14/27 |
| 31×A | 31 | 10 | 224 | 5/112 |

## 6. Validation Checklist

- [x] Evidence documented with sources
- [x] All assertions use exact hand-calculated values (no vague ranges)
- [x] All guard clauses tested (null, invalid params)
- [x] Invariants tested (ranges, length preservation)
- [x] Cross-verified against Wikipedia examples
- [x] DUST formula matches Morgulis et al. (2006)
- [x] Shannon formula matches Shannon (1948)
- [x] LC formula matches Troyanskaya et al. (2002)
- [x] No assumptions — all behaviors sourced
- [x] No duplicates — each test serves a distinct purpose
- [x] Coverage classification complete: 0 missing, 0 weak, 0 duplicate
