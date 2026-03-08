# PHYLO-DIST-001: Phylogenetic Distance Matrix - Test Specification

## Test Unit Identification
- **ID:** PHYLO-DIST-001
- **Title:** Distance Matrix Calculation
- **Area:** Phylogenetic Analysis
- **Canonical Methods:**
  - `PhylogeneticAnalyzer.CalculateDistanceMatrix(seqs, method)`
  - `PhylogeneticAnalyzer.CalculatePairwiseDistance(s1, s2, method)`
- **Complexity:** O(n² × m) where n = number of sequences, m = sequence length

---

## Test Categories

### MUST Tests (Core Functionality - Evidence-Based)

| Test ID | Test Name | Description | Evidence Source |
|---------|-----------|-------------|-----------------|
| M01 | IdenticalSequences_ReturnsZero | d(seq, seq) = 0 for all methods | Mathematical definition |
| M02 | SymmetricMatrix | Matrix[i,j] = Matrix[j,i] | Time-reversibility (Wikipedia) |
| M03 | ZeroDiagonal | Matrix[i,i] = 0 | Distance matrix property |
| M04 | NonNegativeDistances | All distances ≥ 0 | Metric property |
| M05 | PDistance_ReturnsProportionDifferent | p = differences / sites | Wikipedia formula |
| M06 | Hamming_ReturnsRawCount | Count of mismatches | Definition |
| M07 | JukesCantor_GreaterThanPDistance | JC69 ≥ p-distance | Correction increases distance |
| M08 | JukesCantor_KnownValue | Verify formula with calculated value | JC69 formula |
| M09 | Kimura2P_TransitionsVsTransversions | Formula verification + distinguishes mutation types | K80 formula (Wikipedia, Kimura 1980) |
| M10 | UnequalLength_ThrowsException | Pre-condition validation | API contract |
| M11 | GapsIgnored | '-' positions excluded from comparison | Wikipedia |
| M12 | CaseInsensitive | 'acgt' same as 'ACGT' | Standard practice |
| M13 | HighDivergence_JC_ReturnsInfinity | When p ≥ 0.75, JC → +∞ | JC69 saturation |
| M14 | AmbiguousBases_SkippedLikeGaps | N, R, Y etc. excluded from comparison | IUPAC standard / pairwise deletion |
| M15 | NullSequence_ThrowsArgumentNullException | Pre-condition validation | API contract |

### SHOULD Tests (Quality and Robustness)

| Test ID | Test Name | Description | Rationale |
|---------|-----------|-------------|-----------|
| S01 | MatrixDimensions_MatchSequenceCount | n×n matrix for n sequences | API contract |
| S02 | ThreeSequences_CorrectMatrix | Verify complete 3×3 matrix | Integration test |
| S03 | Kimura_GreaterThanPDistance | K2P also increases distance | Correction property |
| S04 | AllGaps_ReturnsZero | No comparable sites → 0 distance | Edge case |
| S05 | SingleDifference_CalculatesCorrectly | 1 mismatch in 8 → p = 0.125 | Boundary case |
| S06 | Kimura2P_MixedChanges_MatchesFormula | K2P with S=1/8, V=1/8 | Formula verification (Kimura 1980) |
| S07 | Kimura2P_HighDivergence_ReturnsInfinity | K2P saturation: V≥0.5 → +∞ | K2P formula domain (Wikipedia) |

### COULD Tests (Extended Scenarios)

| Test ID | Test Name | Description | Rationale |
|---------|-----------|-------------|-----------|
| C01 | LargeMatrix_Performance | 100 sequences performance | Scalability |
| C02 | AllMethods_ConsistentOrder | JC ≥ p-dist for same input | Property consistency |
| C03 | MixedGaps_CorrectHandling | Multiple gaps in alignment | Complex gap patterns |

---

## Test Implementation Notes

### Current Test Coverage (Audit Complete)

**Canonical file:** `PhylogeneticAnalyzer_DistanceMatrix_Tests.cs` — 26 tests total

| Spec ID | Test Name | Status |
|---------|-----------|--------|
| M01 | CalculatePairwiseDistance_IdenticalSequences_ReturnsZeroForAllMethods | ✅ Covered |
| M02 | CalculateDistanceMatrix_IsSymmetric | ✅ Covered |
| M03 | CalculateDistanceMatrix_DiagonalIsZero | ✅ Covered |
| M04 | CalculateDistanceMatrix_AllValuesNonNegative | ✅ Covered (finite-distance sequences, NaN check) |
| M05 | CalculatePairwiseDistance_PDistance_ReturnsProportionDifferent | ✅ Covered |
| M06 | CalculatePairwiseDistance_Hamming_ReturnsRawCount | ✅ Covered |
| M07 | CalculatePairwiseDistance_JukesCantor_GreaterThanOrEqualToPDistance | ✅ Covered |
| M08 | CalculatePairwiseDistance_JukesCantor_MatchesFormula | ✅ Covered (formula + hardcoded reference 0.13674) |
| M09 | CalculatePairwiseDistance_Kimura2P_DistinguishesTransitionTypes | ✅ Covered (formula + hardcoded references 0.34657, 0.31713) |
| M10 | CalculatePairwiseDistance_UnequalLengths_ThrowsArgumentException | ✅ Covered |
| M11 | CalculatePairwiseDistance_GapsIgnored | ✅ Covered |
| M12 | CalculatePairwiseDistance_CaseInsensitive | ✅ Covered |
| M13 | CalculatePairwiseDistance_JukesCantor_HighDivergence_ReturnsInfinity | ✅ Covered |
| M14 | CalculatePairwiseDistance_AmbiguousBases_SkippedLikeGaps | ✅ Covered |
| M15 | CalculatePairwiseDistance_NullSequence_ThrowsArgumentNullException | ✅ Covered |
| S01 | CalculateDistanceMatrix_DimensionsMatchSequenceCount | ✅ Covered |
| S02 | CalculateDistanceMatrix_ThreeSequences_CorrectValues | ✅ Covered |
| S03 | CalculatePairwiseDistance_Kimura2P_GreaterThanOrEqualToPDistance | ✅ Covered |
| S04 | CalculatePairwiseDistance_AllGaps_ReturnsZero | ✅ Covered |
| S05 | CalculatePairwiseDistance_SingleDifference_CorrectPDistance | ✅ Covered |
| S06 | CalculatePairwiseDistance_Kimura2P_MixedChanges_MatchesFormula | ✅ Covered |
| S07 | CalculatePairwiseDistance_Kimura2P_HighDivergence_ReturnsInfinity | ✅ Covered |
| C01 | CalculateDistanceMatrix_100Sequences_CompletesSuccessfully | ✅ Covered |
| C02 | CalculatePairwiseDistance_AllMethods_ConsistentOrdering | ✅ Covered |
| C03 | CalculatePairwiseDistance_MixedGaps_CorrectHandling | ✅ Covered |

**Additional verification tests:** BothGapsAtSamePosition, JC_TwoDifferences (p=0.25), Pairwise Symmetric

### Remediation Summary
1. ~~**Weak:** `CompletelyDifferent_ReturnsHigh`~~ → Removed (replaced by M08 with exact formula + hardcoded value)
2. ~~**Weak:** `Kimura2Parameter_CalculatesCorrectly`~~ → Removed (replaced by M09 with formula + hardcoded values)
3. ~~**Duplicate:** 4 distance tests in `PhylogeneticProperties.cs`~~ → Removed (M01-M04 duplicates)
4. **Strengthened:** M04 uses finite-distance sequences (not saturated), M08/M09 include hardcoded reference values
5. **Added:** C01 (100-seq performance), C02 (JC≥p, K2P≥p ordering), C03 (complex gap patterns)
6. **Implementation fix:** Ambiguous IUPAC bases (N, R, Y, W, S, M, K, B, D, H, V) now skipped like gaps
7. **Implementation fix:** Added `ArgumentNullException` for null sequences
8. **Added:** M14 (ambiguous base handling), M15 (null validation)

---

## Invariants to Verify

1. **Symmetry:** For all i,j: `Matrix[i,j] == Matrix[j,i]` — time-reversibility (Wikipedia: Substitution model)
2. **Zero Diagonal:** For all i: `Matrix[i,i] == 0` — identity property
3. **Non-Negative:** For all i,j: `Matrix[i,j] >= 0` — metric property
4. **Correction Ordering:** `JC >= PDistance` and `K2P >= PDistance` for same sequences
5. **JC69 Saturation:** When p ≥ 3/4, JC returns +∞ (formula domain: 1-4p/3 ≤ 0)
6. **K2P Saturation:** When 1-2S-V ≤ 0 or 1-2V ≤ 0, K2P returns +∞

---

## Deviations and Assumptions

None. All formulas, behaviors, and edge cases are sourced from external references:

- **JC69 formula:** $d = -\frac{3}{4}\ln(1-\frac{4}{3}p)$ — Wikipedia (Models of DNA evolution, JC69 section); Jukes & Cantor (1969)
- **K2P formula:** $K = -\frac{1}{2}\ln((1-2p-q)\sqrt{1-2q})$ — Wikipedia (Models of DNA evolution, K80 section); Kimura (1980)
- **Transition classification:** A↔G (purine↔purine), C↔T (pyrimidine↔pyrimidine) — Wikipedia (Transition genetics)
- **Saturation:** JC69 at $p \geq 3/4$, K2P when $1-2p-q \leq 0$ or $1-2V \leq 0$ — formula domain constraints
- **Gap handling:** Gaps and ambiguous IUPAC bases excluded from comparison — Wikipedia (Distance matrices in phylogeny); only standard bases (A, C, G, T) are compared
- **Zero-site edge case:** Returns 0 when comparableSites = 0 — mathematical limit (0 differences / any denominator → 0)
- **Case insensitivity:** Standard bioinformatics practice

---

## References
- Evidence Document: [PHYLO-DIST-001-Evidence.md](../docs/Evidence/PHYLO-DIST-001-Evidence.md)
- Algorithm Doc: [docs/algorithms/Phylogenetics/Distance_Matrix.md](../docs/algorithms/Phylogenetics/Distance_Matrix.md)
