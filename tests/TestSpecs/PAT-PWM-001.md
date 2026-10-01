# Test Specification: PAT-PWM-001

**Test Unit ID:** PAT-PWM-001
**Area:** Pattern Matching
**Algorithm:** Position Weight Matrix (PWM)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-03-01

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| Source | URL | Accessed |
|--------|-----|----------|
| Wikipedia: Position weight matrix | https://en.wikipedia.org/wiki/Position_weight_matrix | 2026-01-22 (verified 2026-03-01) |
| Kel et al. (2003) MATCH: A tool for searching TF binding sites | https://pmc.ncbi.nlm.nih.gov/articles/PMC169193/ | 2026-01-22 |
| Rosalind: Consensus and Profile | https://rosalind.info/problems/cons/ | 2026-01-22 (verified 2026-03-01) |
| Nishida et al. (2008) Pseudocounts for transcription factor binding sites | Nucleic Acids Research 37(3):939-944 | Reference |
| Stormo (2000) DNA binding sites: representation and discovery | Bioinformatics Review | Reference |

### 1.2 Algorithm Description

#### Position Weight Matrix (Wikipedia)

> A position weight matrix (PWM), also known as a position-specific weight matrix (PSWM)
> or position-specific scoring matrix (PSSM), is a commonly used representation of
> motifs (patterns) in biological sequences.

**Construction Process (Wikipedia):**

1. **Position Frequency Matrix (PFM):** Count occurrences of each nucleotide at each position
2. **Position Probability Matrix (PPM):** Normalize by dividing by number of sequences
3. **Position Weight Matrix (PWM):** Convert to log-odds using background model

**Log-odds Formula (Wikipedia):**
$$M_{k,j} = \log_2\left(\frac{M_{k,j}}{b_k}\right)$$

Where:
- $M_{k,j}$ is the probability of nucleotide k at position j
- $b_k$ is the background frequency (typically 0.25 for DNA)

**Pseudocounts (Wikipedia/Nishida et al.):**
> Pseudocounts (or Laplace estimators) are often applied when calculating PPMs if
> based on a small dataset, in order to avoid matrix entries having a value of 0.

**Scoring (Wikipedia):**
> When the PWM elements are calculated using log likelihoods, the score of a sequence
> can be calculated by adding (rather than multiplying) the relevant values at each
> position in the PWM. The sequence score gives an indication of how different the
> sequence is from a random sequence.

**Score Interpretation:**
- Score = 0: Equal probability of being functional vs random site
- Score > 0: More likely to be a functional site
- Score < 0: More likely to be a random site

#### Profile Matrix (Rosalind CONS Problem)

> A profile matrix is a 4×n matrix P in which P_{1,j} represents the number of
> times that 'A' occurs in the jth position of one of the strings.

**Consensus String (Rosalind):**
> A consensus string c is a string of length n formed from our collection by taking
> the most common symbol at each position.

### 1.3 Reference Examples from Evidence

#### Wikipedia PPM Example

**Input Sequences (10 sequences of length 9):**
```
GAGGTAAAC
TCCGTAAGT
CAGGTTGGA
ACAGTCAGT
TAGGTCATT
TAGGTACTG
ATGGTAACT
CAGGTATAC
TGTGTGAGT
AAGGTAAGT
```

**Expected PPM (Wikipedia):**
| Base | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|------|---|---|---|---|---|---|---|---|---|
| A | 0.3 | 0.6 | 0.1 | 0.0 | 0.0 | 0.6 | 0.7 | 0.2 | 0.1 |
| C | 0.2 | 0.2 | 0.1 | 0.0 | 0.0 | 0.2 | 0.1 | 0.1 | 0.2 |
| G | 0.1 | 0.1 | 0.7 | 1.0 | 0.0 | 0.1 | 0.1 | 0.5 | 0.1 |
| T | 0.4 | 0.1 | 0.1 | 0.0 | 1.0 | 0.1 | 0.1 | 0.2 | 0.6 |

#### Rosalind CONS Example

**Input Sequences (7 sequences of length 8):**
```
ATCCAGCT
GGGCAACT
ATGGATCT
AAGCAACC
TTGGAACT
ATGCCATT
ATGGCACT
```

**Expected Profile Matrix:**
| Base | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|------|---|---|---|---|---|---|---|---|
| A | 5 | 1 | 0 | 0 | 5 | 5 | 0 | 0 |
| C | 0 | 0 | 1 | 4 | 2 | 0 | 6 | 1 |
| G | 1 | 1 | 6 | 3 | 0 | 1 | 0 | 0 |
| T | 1 | 5 | 0 | 0 | 0 | 1 | 1 | 6 |

**Expected Consensus:** ATGCAACT

### 1.4 Edge Cases from Evidence

| Edge Case | Expected Behavior | Source |
|-----------|-------------------|--------|
| Empty alignment | Exception | Implementation contract |
| Single sequence | Valid PWM | Mathematical definition |
| Unequal length sequences | Exception | Wikipedia (same length required) |
| All same base at position | Log-odds = log2(1/0.25) ≈ 2.0 | Wikipedia formula |
| Pseudocount = 0 | Risk of -∞ for unseen bases | Wikipedia |
| Default pseudocount 0.25 | Avoids zero probabilities | Nishida et al. |
| Threshold at boundary | Include/exclude based on >= | Implementation |
| High threshold | Few or no matches | Definition |
| Sequence shorter than PWM | No matches | Definition |
| Non-ACGT characters in training | ArgumentException | IUPAC-IUB: only A,C,G,T are valid bases |
| Lowercase input | Case-insensitive matching | Guaranteed: CreatePwm calls ToUpperInvariant() |
| Invalid IUPAC code in training | ArgumentException | Only standard bases accepted |
| Null input | ArgumentNullException | Implementation contract |

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CreatePwm(IEnumerable<string>, double)` | MotifFinder | **Canonical** | PWM construction |
| `ScanWithPwm(DnaSequence, PWM, double)` | MotifFinder | **Canonical** | Sequence scanning |

### 2.1 Supporting Types

| Type | Description |
|------|-------------|
| `PositionWeightMatrix` | PWM data structure with Matrix, Length, Consensus, MaxScore, MinScore |
| `MotifMatch` | Match result with Position, MatchedSequence, Pattern, Score |

---

## 3. Test Classification

### 3.1 Must Tests (Evidence-Based)

| ID | Test Case | Expected | Source | Test Method |
|----|-----------|----------|--------|-------------|
| M1 | Single sequence creates valid PWM | Length = sequence length, Consensus = sequence | Mathematical definition | `CreatePwm_SingleSequence_CreatesValidMatrix` |
| M2 | Multiple identical sequences create same PWM | Consensus = common sequence | Definition | `CreatePwm_MultipleIdenticalSequences_CreatesSameConsensus` |
| M3 | Consensus derives from highest scoring base at each position | Follows max rule | Wikipedia/Rosalind | `CreatePwm_MixedSequences_ConsensusFollowsMaxRule` |
| M4 | Empty sequences throws ArgumentException | Exception | Contract | `CreatePwm_EmptySequences_ThrowsArgumentException` |
| M5 | Unequal length sequences throws ArgumentException | Exception | Wikipedia requirement | `CreatePwm_UnequalLengths_ThrowsArgumentException` |
| M6 | PWM MaxScore > MinScore for non-uniform matrix | MaxScore > MinScore (strict) | Definition | `Pwm_MaxScore_GreaterThanMinScore_ForNonUniform` |
| M7 | Log-odds scores: perfect match gets max score at position | Verified numerically | Wikipedia formula | `Pwm_LogOdds_PerfectMatchGetsMaxPositionalScore` |
| M8 | ScanWithPwm finds exact trained sequence | Match found | Definition | `ScanWithPwm_FindsTrainedSequence` |
| M9 | Threshold filters results correctly | Only scores >= threshold, actual removal verified | Definition | `ScanWithPwm_ThresholdFiltersResults` |
| M10 | Null sequence in CreatePwm throws | ArgumentNullException | Contract | `CreatePwm_NullSequences_ThrowsArgumentNullException` |
| M11 | Null sequence in ScanWithPwm throws | ArgumentNullException | Contract | `ScanWithPwm_NullSequence_ThrowsArgumentNullException` |
| M12 | Null PWM in ScanWithPwm throws | ArgumentNullException | Contract | `ScanWithPwm_NullPwm_ThrowsArgumentNullException` |
| M13 | PWM length matches input sequence length | pwm.Length = sequences[0].Length | Definition | `Pwm_Length_MatchesInputSequenceLength` |
| M14 | Rosalind CONS consensus test case | Consensus = "ATGCAACT" | Rosalind | `CreatePwm_RosalindCONS_TestCase` |
| M15 | Non-ACGT characters in training sequences throw | ArgumentException | IUPAC-IUB 1970: only A,C,G,T | `CreatePwm_NonAcgtCharacters_ThrowsArgumentException` |
| M16 | Sequence shorter than PWM returns no matches | Empty result | Definition | `ScanWithPwm_SequenceShorterThanPwm_ReturnsEmpty` |
| M17 | Log-odds formula numerical verification | Exact matrix values match hand-calculation | Wikipedia formula | `Pwm_LogOddsFormula_NumericalVerification` |
| M18 | Wikipedia 10-sequence example produces correct consensus and scores | Consensus="TAGGTAAGT", key log-odds verified | Wikipedia | `CreatePwm_WikipediaExample_ConsensusAndScoresMatchSource` |
| M19 | Scanning score equals sum of positional log-odds | Score = Σ PWM[base,pos] | Wikipedia | `ScanWithPwm_ScoreEqualsSumOfLogOdds` |
| M20 | Pseudocount = 0 produces -∞ for unseen bases | -∞ for unseen, 2.0 for observed | Wikipedia edge case | `CreatePwm_ZeroPseudocount_ProducesNegativeInfinity` |
| M21 | Threshold boundary: score == threshold is included | >= semantics, not > | Definition | `ScanWithPwm_ThresholdBoundary_ExactScoreIncluded` |

### 3.2 Should Tests (Recommended)

| ID | Test Case | Expected | Source | Test Method |
|----|-----------|----------|--------|-------------|
| S1 | Pseudocount prevents zero probabilities | All scores finite | Wikipedia/Nishida | `Pwm_Pseudocount_PreventsInfiniteScores` |
| S2 | Case-insensitive input handling | Uppercase normalization | Guaranteed | `CreatePwm_LowercaseInput_NormalizesToUppercase` |
| S3 | Multiple matches returned in position order | Sorted by position | Definition | `ScanWithPwm_ReturnsCorrectPositions` |
| S4 | MatchedSequence property populated correctly | Substring at match position | Definition | `ScanWithPwm_ReturnsMatchedSequence` |
| S5 | Score within valid range [MinScore, MaxScore] | All match scores in bounds | Wikipedia | `ScanWithPwm_ScoreWithinValidRange` |
| S6 | High threshold (MaxScore) returns only perfect matches | Strict count reduction verified | Definition | `ScanWithPwm_HighThreshold_ReturnsFewerMatches` |

### 3.3 Edge Case Tests

| ID | Test Case | Expected | Test Method |
|----|-----------|----------|-------------|
| E1 | Uniform sequences: exact log-odds values | N=3, freq(A)=log2(3.25), others=-2.0 | `CreatePwm_UniformSequence_ExactLogOddsValues` |
| E2 | Maximum entropy: equal frequency → all log-odds = 0 | MaxScore = MinScore = 0 | `CreatePwm_MaximumEntropy_AllLogOddsZero` |
| E3 | Multiple overlapping matches all returned | AA in AAAAAA → 5 matches | `ScanWithPwm_MultipleMatchesSameScore_AllReturned` |

### 3.4 Invariant Tests

| ID | Test Case | Test Method |
|----|-----------|-------------|
| I1 | All PWM invariants hold (Length, Consensus, MaxScore≥MinScore, Matrix dims, finite scores, valid bases) | `Pwm_AllInvariants_HoldForValidInput` |

### 3.5 Could Tests (Optional)

| ID | Test Case | Expected | Notes |
|----|-----------|----------|-------|
| C1 | Large PWM performance | Reasonable time | Performance |
| C2 | Genome-scale scanning | Memory efficient | Performance |

---

## 4. Invariants

| Invariant | Description | Test Method |
|-----------|-------------|-------------|
| **PWM Length** | PWM.Length equals input sequence length | Assert.That(pwm.Length, Is.EqualTo(seqLength)) |
| **Consensus Length** | Consensus length equals PWM length | Assert.That(pwm.Consensus.Length, Is.EqualTo(pwm.Length)) |
| **MaxScore >= MinScore** | Maximum always >= minimum | Assert.That(pwm.MaxScore, Is.GreaterThanOrEqualTo(pwm.MinScore)) |
| **Score Range** | All match scores between MinScore and MaxScore | Assert within bounds |
| **Matrix Dimensions** | Matrix is 4 x Length | Assert.That(pwm.Matrix.GetLength(0), Is.EqualTo(4)) |

---

## 5. Contracts

| ID | Contract | Guarantee Source |
|----|----------|------------------|
| C1 | Background frequency = 0.25 for each base | Wikipedia: "0.25 for nucleotides" (uniform) |
| C2 | Pseudocount is a configurable parameter (default 0.25) | API: `CreatePwm(sequences, pseudocount: 0.25)` |
| C3 | Non-ACGT characters in training sequences → ArgumentException | Strict validation (IUPAC-IUB: only A,C,G,T defined) |
| C4 | Case-insensitive input | Guaranteed: `ToUpperInvariant()` in CreatePwm |
| C5 | PWM formula: log2((count + p) / (N + 4p) / b) | Wikipedia log-odds; equals Biopython `counts.normalize(pseudocounts=p).log_odds(background)` |
| C7 | pseudocount finite and ≥ 0, else ArgumentOutOfRangeException; null element → ArgumentException | Nishida 2008 (pseudocounts are non-negative counts) |
| C8 | Background overload: 4 positive finite values, normalised to sum 1 | Biopython `log_odds(background=...)` |
| C9 | Forward-strand scan; reverse strand via `PositionWeightMatrix.ReverseComplement()` | Biopython `reverse_complement`, `search(both=True)` |
| C6 | Score = sum of positional log-odds | Wikipedia: "adding (rather than multiplying) the relevant values" |

---

## 6. Source Verification Audit (2026-03-01)

### 6.1 Sources Fetched

1. **Wikipedia: Position weight matrix** — Full article fetched. Confirms: PFM → PPM → PWM construction process, log-odds formula `log2(PPM / b_k)`, background b_k=0.25 for DNA, pseudocounts recommendation, additive scoring.
2. **Rosalind: Consensus and Profile (CONS)** — Problem page fetched. Confirms: profile matrix = 4×n count matrix, consensus = most common symbol at each position. Sample dataset (7 sequences of length 8) and expected output verified.

### 6.2 Cross-Reference Results

| Element | Wikipedia | Rosalind | Implementation | Tests |
|---------|-----------|----------|----------------|-------|
| PFM counting | ✅ indicator function | ✅ count matrix | ✅ count loop | ✅ Rosalind test verifies consensus |
| PPM normalization | ✅ divide by N | N/A (counts only) | ✅ `(count+p)/(N+4p)` | ✅ M17 numerical test |
| Log-odds formula | ✅ `log2(PPM/b_k)` | N/A | ✅ `Math.Log2(freq/0.25)` | ✅ M17, M18 exact values |
| Background = 0.25 | ✅ "0.25 for nucleotides" | N/A | ✅ hardcoded 0.25 | ✅ M17 verifies |
| Pseudocounts | ✅ "Laplace estimators" | N/A | ✅ configurable, default 0.25 | ✅ S1 finiteness |
| Additive scoring | ✅ "adding...the relevant values" | N/A | ✅ `score += pwm.Matrix[baseIndex, j]` | ✅ M19 sum verification |
| Consensus = max at each position | ✅ implied by PPM | ✅ "most common symbol" | ✅ GenerateConsensus | ✅ M14 Rosalind, M18 Wikipedia |
| Wikipedia PPM example | ✅ 10 sequences, 9 positions | N/A | N/A (reference only) | ✅ M18: 10 sequences, consensus + key scores |
| Rosalind CONS example | N/A | ✅ 7 seqs, consensus "ATGCAACT" | N/A | ✅ M14 exact match |

### 6.3 Discrepancies Found and Fixed

| Discrepancy | Source | Fix |
|-------------|--------|-----|
| Spec listed "7 sequences of length 9" for Wikipedia example | Wikipedia uses 10 sequences | Fixed spec to list correct 10 sequences |
| Spec sequences 5–7 were incorrectly transcribed | Verified against Wikipedia PFM | Fixed to match Wikipedia article |
| CreatePwm silently ignored non-ACGT characters (A3) | No source justifies silent skip | Changed to throw ArgumentException |
| Test M15 tested non-ACGT via DnaSequence (which rejects non-ACGT) | Dead test | Replaced with CreatePwm input validation test |
| Algorithm doc formula showed `(PPM + p) / b_k` | Pseudocount applies to PFM, not PPM | Fixed formula in Position_Weight_Matrix.md |

---

## 7. Coverage Classification Audit

### 7.1 Summary

| Category | Count | Actions Taken |
|----------|-------|---------------|
| ❌ Missing | 2 | Added M20 (pseudocount=0 → -∞) and M21 (threshold boundary >= semantics) |
| ⚠ Weak | 6 | Strengthened with exact values and concrete assertions |
| 🔁 Duplicate | 2 | Removed `Pwm_Consensus_HasCorrectLength` (⊂ I1) and `Pwm_Matrix_HasCorrectDimensions` (⊂ I1) |
| ✅ Covered | 21 | No changes needed |

### 7.2 Missing Tests — Implemented

| ID | Test | What Was Missing |
|----|------|-----------------|
| M20 | `CreatePwm_ZeroPseudocount_ProducesNegativeInfinity` | Edge case documented in spec but untested: pseudocount=0 → -∞ for unseen bases |
| M21 | `ScanWithPwm_ThresholdBoundary_ExactScoreIncluded` | Threshold >= boundary never explicitly verified; score == threshold must be included |

### 7.3 Weak Tests — Strengthened

| Test | Problem | Fix |
|------|---------|-----|
| `Pwm_MaxScore_GreaterThanMinScore_ForNonUniform` | Used `>=` instead of `>` for non-uniform matrices | Changed to strict `Is.GreaterThan` |
| `ScanWithPwm_ThresholdFiltersResults` | Never verified filtering actually removed matches | Added assertion: `filteredMatches.Count < allMatches.Count` with predictable data |
| `ScanWithPwm_HighThreshold_ReturnsFewerMatches` | `LessThanOrEqualTo` allowed equal count | Changed to strict `GreaterThan`, exact count assertions, exact positions |
| `ScanWithPwm_ReturnsCorrectPositions` | `EquivalentTo` (set equality) ignored position order | Changed to `Is.EqualTo` — verifies sorted position order |
| `CreatePwm_UniformSequence_ExactLogOddsValues` | Only checked `MaxScore - MinScore > 0` | Full exact values: log2(3.25) for observed, -2.0 for unseen, at every position |
| `CreatePwm_MaximumEntropy_AllLogOddsZero` | Only checked Length/Consensus.Length/MaxScore>=MinScore | All 16 matrix cells asserted = 0.0; MaxScore = MinScore = 0.0 |

### 7.4 Duplicate Tests — Removed

| Test | Reason |
|------|--------|
| `Pwm_Consensus_HasCorrectLength` | Subset of `Pwm_AllInvariants_HoldForValidInput` (Invariant 2) |
| `Pwm_Matrix_HasCorrectDimensions` | Subset of `Pwm_AllInvariants_HoldForValidInput` (Invariant 4) |

### 7.5 Final Test Count

| File | Category | Count |
|------|----------|-------|
| `MotifFinder_PWM_Tests.cs` | Canonical (PAT-PWM-001) | 31 |
| `MotifFinderTests.cs` | Smoke | 2 |
| `PatternMatchingProperties.cs` | Property | 3 |
| `PatternMatchingSnapshotTests.cs` | Snapshot | 1 |
| **Total** | | **37** |

---

## 8. Review 2026-09 (B05) — Biopython-locked tests

Reference: Biopython 1.88 `Bio.motifs` on the Wikipedia 10-sequence example, target `CCTAGGTAAGTAACAGGTCAGTGG`.

| ID | Test | Locked values |
|----|------|---------------|
| B1 | `CreatePwm_WikipediaExample_MatrixEqualsBiopython` | all 36 cells of `normalize(0.25).log_odds()`; max 11.707530265108911, min −14.88138790352414 |
| B2 | `ScanWithPwm_WikipediaExample_EqualsBiopythonSearch` | `search(both=False)` → (2, 11.7075), (6, 4.7792), (13, 9.3161); all 16 `calculate` scores |
| B3 | `ReverseComplement_WikipediaExample_EqualsBiopython` | `reverse_complement()` matrix; minus-strand hit at 8 (Biopython −16), 2.38769 |
| B4 | `CreatePwm_NonUniformBackground_EqualsBiopython` | `log_odds(background={A:.3,C:.2,G:.2,T:.3})` matrix; max 11.095108114768236, min −16.07877255458597; hits 2/6/13 |
| B5 | `CreatePwm_ZeroPseudocount_MaxMinEqualBiopython` | max 12.438028776954503, min −∞ |
| B6 | `CreatePwm_InvalidPseudocount_Throws` | −0.1/−0.5/NaN/+∞ → ArgumentOutOfRangeException (previously NaN cells) |
| B7 | `CreatePwm_NullElement_ThrowsArgumentException`, `CreatePwm_InvalidBackground_Throws`, `PositionWeightMatrix_Constructor_ValidatesShape`, `CreatePwm_UniformBackgroundOverload_EqualsDefault` | contracts |

## 9. Review 2026-09 (B05 follow-up) — both strands, calculate, FromCounts, score distribution

Tests: `Unit/Analysis/MotifFinder_PwmStrandsAndThresholds_Tests.cs`, `Metamorphic/MotifPwmMetamorphicTests.cs` (PWM-BOTH, PWM-DIST). Reference: Biopython 1.88 (`Bio/motifs/matrix.py`, `Bio/motifs/thresholds.py`, `Bio/motifs/jaspar`).

| ID | Test | Locked values / invariant |
|----|------|---------------------------|
| C1 | `ScanWithPwmBothStrands_WikipediaExample_EqualsBiopython` | `search(T, 0, both=True)` → (2, 11.70753), (6, 4.77916), (−16, 2.38769), (13, 9.31606) |
| C2 | `ScanWithPwmBothStrands_LowerThreshold_EqualsBiopython` | threshold −3 adds (−21, −0.79252), (−17, −0.64793) |
| C3 | `ScanWithPwmBothStrands_PalindromicPwm_ReportsBothStrandsPlusFirst` | palindromic PWM: both strands per window, '+' first (Biopython tie order is NumPy-argsort-defined) |
| C4 | `ScanWithPwmBothStrands_RandomCases_EqualBiopython` (seed 20260930, 3 cases) | positions + scores |
| C5 | `CalculatePwmScores_InvalidWindowsNaN_MixedCase_EqualsBiopython`, `_DnaSequence_EqualsBiopythonCalculate_AndShortSequenceEmpty` | `calculate` incl. NaN windows |
| C6 | `MeanStd_WikipediaExample_EqualBiopython` | mean 5.522241422369563, std 3.2186083897634568, mean(bg) 5.975137019521005 |
| C7 | `ScoreDistribution_WikipediaExample_ThresholdsEqualBiopython` | fpr .01 → 4.028388324862519, .001 → 7.10122696197535, fnr .1 → 1.2303323735684302, balanced 0.1430202404361971 (rate 0.06385040283203125), patser 2.5924271925194056 |
| C8 | `ScoreDistribution_BackgroundAndPrecision_EqualBiopython` | bg .3/.2/.2/.3, precision 100 |
| C9 | `PythonFloorDiv_MatchesCPython`, `ScoreDistribution_InvalidInputs_Throw`, `FromCounts_ValidationAndScalarPseudocount` | CPython `//`; contracts; FromCounts ≡ CreatePwm on the same counts |
| C10 | `PwmBothStrands_ReverseComplementMirror` (5 seeds), `PwmScoreDistribution_Monotone` (3 seeds) | mirror relation; FPR/FNR monotone, on-grid, densities sum to 1 |


## 10. Review 2026-09 (B05 audit group D) — exact p-values, per-base / JASPAR pseudocounts

Tests: `Unit/Analysis/MotifFinder_PwmPValue_Tests.cs`; MCP `PwmScorePValueTests`, `CreatePwmTests.CreatePwm_PseudocountOptions_EqualBiopython`.
References: exhaustive enumeration of all 4^L words (C oracle, left-to-right double window sums, long-double accumulation) for the Wikipedia PWM (L = 9) and Bucher TATA box POL012.1 (L = 15); TFM-Pvalue C++ (CRAN TFMPvalue 1.0.0 `src/Matrix.cpp`, driver = `testScoreToPvalue` / `testPvalueToScore` loops); Biopython 1.88 (`counts.normalize(pseudocounts=dict).log_odds(bg)`, `jaspar.calculate_pseudocounts`).

| ID | Test | Locked values / invariant |
|----|------|---------------------------|
| D1 | `PwmScorePValue_ToyMatrix_EqualsBinomialTail` | ±2/−1 matrix: P(S ≥ 6) = 1/64, ≥ 3 = 10/64, ≥ 0 = 37/64, ≥ −3 = 1, > max = 0 |
| D2 | `PwmScorePValue_Wikipedia(_Background)_EqualsExhaustiveEnumeration` | e.g. site 6 score 4.77915994208994 → 0.006160736083984375 (bg .3/.2/.2/.3: 0.007509384); consensus → 4^-9 |
| D3 | `PwmScorePValue_ObservedWindowScores_AreCounted` | scores from `CalculatePwmScores` count their own window (TFM-Pvalue: 0 / 0.006023406982421875) |
| D4 | `PwmScorePValue_BucherTataBox_EqualsExhaustiveEnumeration` | 13.282082872405402 → 6.019137799739838e-06; 0 → 0.02533565554767847 |
| D5 | `PwmScorePValue_Wikipedia_EqualsTfmPvalue`, `PwmScoreThresholdForPValue_PValuesEqualTfmPvalue` | non-tied thresholds = TFM-Pvalue sc2pv; pv2sc p-values identical (TFM score rounded: 7.11 vs 7.150252343998389) |
| D6 | `PwmScoreThresholdForPValue_*_EqualsExhaustiveEnumeration` | smallest word score t with P(S ≥ t) ≤ p, next lower word score has P > p (Wikipedia p = .01 → 4.028050165603465; next 4.028050165603464 one ulp below) |
| D7 | `_BelowConsensusProbability_IsPositiveInfinity`, `_RoundTripsThroughScorePValue`, `PwmScorePValue_ContrastsWithGridDistribution` | +∞/0; round trip; Biopython grid threshold 4.028388324862519 vs exact 4.028050165603465 |
| D8 | `PwmPValue_InvalidArguments_Throw`, `PwmPValue_EmptyMatrix_ScoresZero` | contracts (non-finite matrix → ArgumentException) |
| D9 | `CreatePwm_PerBasePseudocounts_EqualsBiopython`, `CreatePwmWithJasparPseudocounts_(Background|Uniform)_EqualsBiopython`, `CreatePwm_PerBase_EqualsScalarAndFromCounts`, `CreatePwm_PerBase_InvalidArguments_Throw` | 20-cell matrices to 1e-12; √7·q pseudocounts |
| D10 | `PwmScoreThresholdForPValue_LargestAchievablePValue_WhereTfmIsNotMaximal` | 8-column random matrix, p = 1e-4: 8.685556713812995 / 9.930678854038951e-05 (TFM-Pvalue: 8.6866 / 9.923731321091039e-05; next lower word has P > 1e-4) |
| D11 | `PwmScorePValue_BudgetExhausted_ReturnsCertifiedBounds` | 20-column random matrix at 0: IsExact = false, bounds [0.058402542608746444, 0.05840680768687889], PValue = upper |

## 11. Review 2026-09 (B05 audit group D, part 2) — generic-alphabet (protein) PWM

Tests: `Unit/Analysis/MotifFinder_AlphabetPwm_Tests.cs`; MCP `AlphabetPwmAndSigma70Tests` (create_alphabet_pwm, scan_with_alphabet_pwm).
Reference: Biopython 1.88 `motifs.create(instances, alphabet).counts.normalize(pseudocounts).log_odds(background)`, `consensus`, `anticonsensus`, `max`, `min`, `mean`, `std`; window scores = Σ_j `pssm[letter][j]` with NaN outside the alphabet (`_pwm.c` rule; Biopython `calculate` itself is DNA-only); scratch cross-check of 400 random cases (≤ 1.1e-14 relative; ACGT cases = Biopython `calculate` float32 exactly).

| ID | Test | Locked values / invariant |
|----|------|---------------------------|
| E1 | `ProteinPwm_ScalarPseudocount_EqualsBiopython` | 6 instances, p = 0.5: W[M,0] 2.7813597135246595, W[K,1] 2.4918530963296748, consensus MKVLAT, anticonsensus AAAACA, max 16.398651663952972, min −4.068431430675826, mean 3.809139711610947, std 3.8318293687371425 |
| E2 | `CalculateScores_ProteinSequence_NaNForUnknownSymbols_EqualsBiopythonSums` | 19 windows of GGMKVLATxxMRvLGTPPMKXLAS incl. 11 NaN; lower case scored |
| E3 | `Scan_Protein_ForwardHitsAtThreshold_EqualsBiopythonSearchBothFalseRule` | threshold 5 → positions 2, 10 (12.939220045315675) |
| E4 | `ProteinPwm_PerSymbolPseudocountsAndBackground_EqualsBiopython` | max 19.061964092904915, min −18.17848406036431, mean(bg) 10.253659399954454, std(bg) 4.8477154676833605 |
| E5 | `ProteinPwm_ZeroPseudocount_UnseenSymbolsNegativeInfinity_EqualsBiopython` | min −∞, max 24.03143403943405, mean 21.429827293694583, std 2.1524130640960037 |
| E6 | `IgnoreUnknownSymbols_SkipsGapsAndX_EqualsBiopythonCounts` | W[M,0] 2.6520766965796927, W[V,3] 3.192645077942396, all-unknown column uniform, consensus MKAV |
| E7 | `DnaAlphabet_ScoresEqualDnaPwm`, `FromCounts_EqualsCreateFromInstances` | ACGT generic ≡ DNA PWM (shared kernels) |
| E8 | `ShortSequence_EmptyScoresAndNoHits`, `GetMatrix_ReturnsDefensiveCopy`, `Guards` | contracts |
