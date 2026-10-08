# TestSpec: SEQ-COMPLEX-001 - Sequence Complexity Metrics

**Test Unit ID:** SEQ-COMPLEX-001
**Area:** Sequence Composition
**Created:** 2026-01-22
**Status:** Complete

## 1. Scope

This TestSpec covers all sequence complexity metrics in `SequenceComplexity`:

| Method | Purpose |
|--------|---------|
| `CalculateLinguisticComplexity` | Vocabulary richness (Orlov & Potapov 2004 summation, word length ≤ m; = Troyanskaya 2002 when m ≥ N) |
| `CalculateShannonEntropy` | Information content per base |
| `CalculateKmerEntropy` | Entropy of k-mer frequency distribution |
| `CalculateWindowedComplexity` | Sliding-window complexity profile |
| `FindLowComplexityRegions` | Low-entropy region detection |
| `CalculateDustScore` | DUST low-complexity score (Σc(c−1)/2 / (ℓ−1), triplets only; see SEQ-COMPLEX-DUST-001) |
| `MaskLowComplexity` | SDUST perfect-interval masking (port of lh3/sdust; see SEQ-COMPLEX-DUST-001) |
| `EstimateCompressionRatio` | Normalized Lempel–Ziv (1976) complexity (= `CalculateNormalizedLempelZivComplexity`; see SEQ-COMPLEX-COMPRESS-001) |

## 2. Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| Wikipedia "Linguistic sequence complexity" | Encyclopedia | LC formula, vocabulary usage, examples |
| Troyanskaya et al. (2002) Bioinformatics 18(5):679–88 | Peer-reviewed | Summation formula: LC = Σ(observed) / Σ(possible) |
| Orlov & Potapov (2004) NAR 32:W628–W633 | Peer-reviewed | V_max = min(4^i, N−i+1); CL = ΣV_i/ΣV_max,i with word length limited by m ≤ N |
| Rosalind LING "Linguistic Complexity of a Genome" | Problem set | lc(s) = sub(s)/m(a,n) over all lengths; sample `ATTTGGATT` → 0.875 |
| universalmotif `sequence_complexity` (R, bjmt/universalmotif, `R/sequence_complexity.R`) | Reference implementation | Trifonov (1990) **product** form C = ΠU_i; doc examples reproduce our V_i / V_max,i to 4 dp |
| Wikipedia "Entropy (information theory)" | Encyclopedia | Shannon formula H = −Σ p_i log₂ p_i |
| Shannon (1948) Bell System Technical Journal | Original paper | Information entropy definition |
| Morgulis et al. (2006) J Comput Biol 13(5):1028–40 | Peer-reviewed | Symmetric DUST originator: S = Σc_t(c_t−1)/2 / (ℓ−1), ℓ = number of triplets (N−2). Implemented exactly (lh3/sdust `find_perfect`, NCBI dustmasker `thresholds_[i] = i·level`; review-2026-09 B04 F2 replaced the earlier ℓ divisor taken from the longdust README). |

## 3. Formulas

### 3.1 Linguistic Complexity (Orlov & Potapov 2004; Troyanskaya 2002)

$$LC = \frac{\sum_{i=1}^{m} V_{obs}(i)}{\sum_{i=1}^{m} V_{max}(i)}$$

Where $V_{max}(i) = \min(a^i, N - i + 1)$, $a$ = alphabet size (4 for DNA/RNA; extended to the symbols present for other input, B04 F20; caller-supplied via the `alphabetSize` overloads, B04 F37) and $m$ = `maxWordLength` (clamped to N). Range: $0 \le LC \le 1$.
With $m \ge N$ this is exactly Troyanskaya et al. (2002) $LC = A(s)/M(s)$ over all lengths (Rosalind LING);
the default $m = 10$ is the Orlov & Potapov (2004) word-length-limited variant. This is **not** Trifonov's
(1990) product $C = \prod_i U_i$ (universalmotif "Trifonov"). For $m > 12$ the $V_i$ are counted from the
suffix tree (Troyanskaya 2002: $V_i$ = number of edges spanning depth $i$), linear in N; identical values.

### 3.2 Shannon Entropy (Shannon 1948)

$$H = -\sum_{i} p_i \log_2 p_i$$

Range for DNA: $0 \le H \le 2$ (log₂4 = 2).

### 3.3 K-mer Entropy

Shannon entropy applied to k-mer frequency distribution. Range: $0 \le H \le \log_2(4^k)$.

### 3.4 DUST Score (Morgulis 2006)

$$DUST = \frac{\sum_t c_t(c_t - 1)/2}{\ell(x) - 1}$$

Where $c_t$ = count of triplet $t$ and $\ell(x) = N - 2$ = the number of overlapping triplets
(fewer than 2 triplets → 0). Divisor $\ell - 1$ per Morgulis et al. (2006), lh3/sdust
`find_perfect` (`new_l = kdq_size(w) − i − 1`, mask iff `r·10 > T·l`) and NCBI dustmasker,
confirmed on the compiled sdust binary (review-2026-09 B04 F2). DUST is defined for triplets
only: `wordSize ≠ 3` is rejected (B04 F34); the sourced k-mer generalisation is longdust
(`CalculateLongdustScore`). *Superseded:* the 2026-06 revision of this spec used the divisor
$\ell(x)$ from the longdust README restatement — that was a regression (see SEQ-COMPLEX-DUST-001).

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
| LC-13 | `CalculateLinguisticComplexity_RosalindLingSample_AllWordLengths_Returns0875` | Exact: 0.875 (m=9, m=int.MaxValue) | Rosalind LING; Troyanskaya (2002) |
| LC-14 | `CalculateLinguisticComplexity_TroyanskayaFullLength_MatchesReference` | Exact: 28/31, 4/31, 33/140, 69/70 (m=N) | Troyanskaya (2002); Python reference |
| LC-15 | `CalculateLinguisticComplexity_OrlovSumForm_UniversalmotifSequences_MatchesReference` | Exact: 51/55, 52/55, 28/55, 9/11, 3/5, 4/11 (m=7) | Orlov & Potapov (2004); V_i cross-checked vs universalmotif |
| LC-16 | `CalculateLinguisticComplexity_WordLengthsBeyond4Pow31_NoOverflow` | Exact: 749/761 (N=m=40) | Formula; no 4^i overflow |
| LC-17 | `CalculateLinguisticComplexity_SuffixTreePath_RepeatRichSequence_MatchesReference` | Exact: 782/1209, 1405/1937, 1971/2251, 6427/6987 (N=120) | Troyanskaya suffix-tree counting; Python reference |
| LC-18 | `CalculateLinguisticComplexity_FullLengthLongHomopolymer_ExactAndLinearTime` | Exact closed form, N=200,000 | Troyanskaya (2002) linear-time claim |
| LC-19 | `CalculateLinguisticComplexity_NonAcgtSymbols_AlphabetExtended_MatchesBruteForce` | Exact: ACGTN 1.0 (was 15/14), ATGCATGCNN 22/25, N×8 8/33, IUPAC hash/suffix-tree 69/142, 345/442, RNA 6/7, ACGTU 0.8 | Troyanskaya (2002) / Rosalind LING alphabet size a; Python brute force (3000 random, 0 mismatches) |
| LC-20 | `CalculateLinguisticComplexity_RnaEqualsDnaCounterpart` | RNA = DNA spelling (a = 4) | IUPAC U ≡ T |
| LC-21 | `CalculateLinguisticComplexity_LargeAlphabet_NeverExceedsOne_NoOverflow` | 300 distinct symbols → 1.0 (m=10, m=300) | LC ≤ 1; saturating a^i |
| LC-22 | `CalculateLinguisticComplexity_FixedAlphabet_MatchesBruteForce` (8 cases, string/lower-case/`DnaSequence`) | ATTTGGATT a=4 → 7/8 (Rosalind LING); AACCAACC a=2 → 8/9, a=4 → 3/4; protein a=20 → 9/10; suffix-tree path a=4 → 138/173, a=6 → 46/59; GATTACA×3 → 41/70; A×8 a=1 → 1 | Troyanskaya 2002 fixed alphabet a; exact-Fraction brute force (4 001 random, 0 mismatches) |
| LC-23 | `CalculateLinguisticComplexity_FixedAlphabet4_EqualsInferredForDna` | a=4 ≡ inferred for pure DNA (m 1..100) | contract |
| LC-24 | `CalculateLinguisticComplexity_FixedAlphabet_InvalidArguments_Throw` / `_HugeAlphabet_NoOverflow` | a < distinct → ArgumentException; a<1, m<1 → ArgumentOutOfRange; null DnaSequence → ArgumentNull; empty/null string → 0; a = int.MaxValue → 1.0 | contract; saturation |

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
| LCR-1 | `FindLowComplexityRegions_FindsPolyARegion` | Count=1, start=79, end=145, minH=0 (F4) | Entropy-based detection |
| LCR-2 | `FindLowComplexityRegions_HighComplexity_ReturnsEmpty` | Empty | Definition |
| LCR-3 | `FindLowComplexityRegions_ReturnsCorrectSequence` | start=6, end=74, length=69 (F4) | Region merging |
| LCR-4 | `FindLowComplexityRegions_NullSequence_ThrowsException` | Throws | Guard clause |
| LCR-5 | `FindLowComplexityRegions_InvalidWindowSize_ThrowsException` | Throws | Guard clause |

### 4.6 DUST Score Tests (5 tests)

| ID | Test Method | Assertion | Source |
|----|-------------|-----------|--------|
| DUST-1 | `CalculateDustScore_LowComplexity_ReturnsHigh` | Exact: 8.0 (N=18, =120/15) | Morgulis (2006) ℓ−1; sdust `find_perfect` (B04 F2) |
| DUST-2 | `CalculateDustScore_HighComplexity_ReturnsLow` | Exact: 6/13 (N=16) | Morgulis (2006) ℓ−1; sdust (B04 F2) |
| DUST-3 | `CalculateDustScore_EmptySequence_ReturnsZero` | Exact: 0 | Convention |
| DUST-4 | `CalculateDustScore_StringOverload_ReturnsExact` | Exact: 2.5 (N=7, =10/4) | Morgulis (2006) ℓ−1; sdust (B04 F2) |
| DUST-5 | `CalculateDustScore_SequenceShorterThanWordSize_ReturnsZero` | Exact: 0 | Boundary |
| DUST-6 | `CalculateDustScore_NullSequence_ThrowsException` | Throws | Guard clause |

### 4.7 Masking Tests (5 tests)

| ID | Test Method | Assertion | Source |
|----|-------------|-----------|--------|
| MASK-1 | `MaskLowComplexity_MasksLowComplexityWindows` | N=192 (all masked) | lh3/sdust `-w 64 -t 20` → [0,192) |
| MASK-2 | `MaskLowComplexity_PreservesHighComplexity` | No N chars | lh3/sdust `-t 100` → no interval |
| MASK-3 | `MaskLowComplexity_CustomMaskChar` | X=100 (all masked) | lh3/sdust `-w 64 -t 10` → [0,100) |
| MASK-4 | `MaskLowComplexity_NullSequence_ThrowsException` | Throws | Guard clause |
| MASK-5 | `MaskLowComplexity_ResultLengthEqualsInputLength` | Length invariant | Definition |
| MASK-6 | `MaskLowComplexity_ShortSequence_PreservesOriginal` | Returns "ATGC" (threshold 0) | raw score 0 is not > 0 (sdust `-t 0`: no output) |

### 4.8 Compression Ratio Tests (SUPERSEDED — see SEQ-COMPLEX-COMPRESS-001)

`EstimateCompressionRatio` was re-implemented as the **normalized Lempel–Ziv (1976)
complexity** $c / (n / \log_b n)$ (reference: entropy/antropy `lziv_complexity`,
Zhang et al. 2009). The earlier "unique-substring ratio" heuristic (asserting 14/27 and
5/112) was removed; those assertions are obsolete. The Lempel–Ziv metric is fully
specified and tested under `tests/TestSpecs/SEQ-COMPLEX-COMPRESS-001.md` /
`SequenceComplexity_EstimateCompressionRatio_Tests.cs`. The only assertions surviving here
are the empty→0 and null→throw guards.

| ID | Test Method | Assertion | Source |
|----|-------------|-----------|--------|
| CR-3 | `EstimateCompressionRatio_EmptySequence_ReturnsZero` | Exact: 0 | Convention |
| CR-5 | `EstimateCompressionRatio_NullSequence_ThrowsException` | Throws | Guard clause |

### 4.9 Edge Cases (2 remaining)

| ID | Test Method | Assertion | Source |
|----|-------------|-----------|--------|
| EDGE-1 | `CalculateLinguisticComplexity_ZeroWordLength_ThrowsException` | Throws | Guard clause |
| EDGE-2 | `CalculateLinguisticComplexity_NegativeWordLength_ThrowsException` | Throws | Guard clause |

## 5. Hand-Calculated Cross-Verification Table

### 5.1 Linguistic Complexity (summation formula; m = maxWord)

| Input | N | maxWord | obs | max | LC | Match |
|-------|---|---------|-----|-----|----|-------|
| `"A"` | 1 | 10 | 1 | 1 | 1.0 | ✓ |
| `"AAAA"` | 4 | 10 | 4 | 10 | 0.4 | ✓ |
| `"ATGC"` | 4 | 10 | 10 | 10 | 1.0 | ✓ |
| `"AAAAAAAAAAAAAAAA"` | 16 | 10 | 10 | 103 | 10/103 | ✓ |
| `"ATGCTAGCATGCAATG"` | 16 | 10 | 91 | 103 | 91/103 | ✓ |
| `"ACGGGAAGCTGATTCCA"` | 17 | 4 | 47 | 49 | 47/49 | ✓ |
| `"ACACACACACACACACA"` | 17 | 10 | 20 | 112 | 5/28 | ✓ |
| `"ATTTGGATT"` (Rosalind) | 9 | ≥9 | 35 | 40 | 0.875 | ✓ |
| `"ATGCTAGCATGCAATG"` | 16 | 16 | 112 | 124 | 28/31 | ✓ |
| `"AAAAAAAAAAAAAAAA"` | 16 | 16 | 16 | 124 | 4/31 | ✓ |

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

| Input | N | Triplets ℓ | Score | DUST = Score/(ℓ−1) |
|-------|---|-----------|-------|----------------|
| `"AAAAAAAAAAAAAAAAAA"` | 18 | AAA×16 | 16·15/2=120 | 120/15=8.0 |
| `"ATGCTAGCATGCTAGC"` | 16 | 14 (6 dups ×2) | 6 | 6/13 |
| `"AAAAAAA"` | 7 | AAA×5 | 5·4/2=10 | 10/4=2.5 |

(Divisor = ℓ − 1 with ℓ = N − 2 triplets — Morgulis et al. 2006, lh3/sdust `find_perfect`
`new_l = kdq_size(w) − i − 1`, NCBI dustmasker; corrected by review-2026-09 B04 F2 from the
earlier ℓ divisor taken from the longdust README. Values = current tests.)

### 5.5 Compression Ratio (SUPERSEDED — normalized Lempel–Ziv, see SEQ-COMPLEX-COMPRESS-001)

`EstimateCompressionRatio` now returns the normalized Lempel–Ziv (1976) complexity
$c/(n/\log_b n)$. The old unique-substring values (14/27, 5/112) no longer apply.
Worked LZ76 values (antropy 0.2.2 `lziv_complexity` doctests / exhaustive-history definition;
corrected by review-2026-09 B04 F1 — the earlier 8 / 9 / 5 values were the Naereen LZ78
incremental parse):

| Input | n | b | c | Normalized LZ |
|-------|---|---|---|---------------|
| `"1001111011000010"` | 16 | 2 | 6 (1/0/01/1110/1100/0010) | 6/(16/log₂16) = 1.5 |
| `"ACGTACGTACGTACGT"` | 16 | 4 | 5 | 5/(16/log₄16) = 0.625 |
| `"0"×16` | 16 | (clamp 2) | 2 | 2/(16/log₂16) = 0.5 |

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
