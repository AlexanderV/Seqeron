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
| M09 | Kimura2P_TransitionsVsTransversions | Distinguishes mutation types | K80 definition |
| M10 | UnequalLength_ThrowsException | Pre-condition validation | API contract |
| M11 | GapsIgnored | '-' positions excluded from comparison | Wikipedia |
| M12 | CaseInsensitive | 'acgt' same as 'ACGT' | Standard practice |
| M13 | HighDivergence_JC_ReturnsInfinity | When p ≥ 0.75, JC → +∞ | JC69 saturation |

### SHOULD Tests (Quality and Robustness)

| Test ID | Test Name | Description | Rationale |
|---------|-----------|-------------|-----------|
| S01 | MatrixDimensions_MatchSequenceCount | n×n matrix for n sequences | API contract |
| S02 | ThreeSequences_CorrectMatrix | Verify complete 3×3 matrix | Integration test |
| S03 | Kimura_GreaterThanPDistance | K2P also increases distance | Correction property |
| S04 | AllGaps_ReturnsZero | No comparable sites → 0 distance | Edge case |
| S05 | SingleDifference_CalculatesCorrectly | 1 mismatch in 8 → p = 0.125 | Boundary case |

### COULD Tests (Extended Scenarios)

| Test ID | Test Name | Description | Rationale |
|---------|-----------|-------------|-----------|
| C01 | LargeMatrix_Performance | 100 sequences performance | Scalability |
| C02 | AllMethods_ConsistentOrder | JC ≥ p-dist for same input | Property consistency |
| C03 | MixedGaps_CorrectHandling | Multiple gaps in alignment | Complex gap patterns |

---

## Test Implementation Notes

### Existing Tests Analysis (Audit)

| Existing Test | Status | Action |
|---------------|--------|--------|
| CalculatePairwiseDistance_IdenticalSequences_ReturnsZero | **Covered** | Keep as M01 |
| CalculatePairwiseDistance_CompletelyDifferent_ReturnsHigh | **Weak** | Replace with specific assertion |
| CalculatePairwiseDistance_PDistance_ReturnsProportionDifferent | **Covered** | Keep as M05 |
| CalculatePairwiseDistance_Hamming_ReturnsRawCount | **Covered** | Keep as M06 |
| CalculatePairwiseDistance_JukesCantor_GreaterThanPDistance | **Covered** | Keep as M07 |
| CalculatePairwiseDistance_Kimura2Parameter_CalculatesCorrectly | **Weak** | Replace with specific assertion |
| CalculatePairwiseDistance_WithGaps_IgnoresGapPositions | **Covered** | Keep as M11 |
| CalculateDistanceMatrix_ReturnsSymmetricMatrix | **Covered** | Keep as M02 |
| CalculateDistanceMatrix_DiagonalIsZero | **Covered** | Keep as M03 |
| CalculatePairwiseDistance_DifferentLength_Throws | **Covered** | Keep as M10 |

### Consolidation Plan
1. **Canonical file:** `PhylogeneticAnalyzer_DistanceMatrix_Tests.cs` (new file for PHYLO-DIST-001)
2. **Existing file:** `PhylogeneticAnalyzerTests.cs` contains distance tests mixed with tree tests
3. **Action:** Extract and refactor distance-related tests into dedicated test file
4. **Strengthen:** Add missing Must tests (M08, M09, M12, M13)
5. **Remove:** Weak tests with vague assertions

---

## Invariants to Verify

1. **Symmetry:** For all i,j: `Matrix[i,j] == Matrix[j,i]`
2. **Zero Diagonal:** For all i: `Matrix[i,i] == 0`
3. **Non-Negative:** For all i,j: `Matrix[i,j] >= 0`
4. **Correction Ordering:** `JC >= PDistance` for same sequences

---

## Open Questions / Decisions

| Question | Decision | Rationale |
|----------|----------|-----------|
| Empty sequences handling? | Return 0 | No differences possible |
| All-gap alignment? | Return 0 | No comparable sites |
| Very high p-distance (≥0.75)? | Return PositiveInfinity | Per JC69 formula saturation |

---

## Test Data Sets

### Standard Test Sequences
```csharp
// Identical
"ACGTACGT" vs "ACGTACGT" → d = 0

// Single difference (1/8 = 0.125 p-distance)
"ACGTACGT" vs "TCGTACGT" → Hamming = 1, p = 0.125, JC ≈ 0.137

// Two differences (2/8 = 0.25 p-distance)  
"ACGTACGT" vs "TCGTACGA" → Hamming = 2, p = 0.25, JC ≈ 0.304

// With gap (7 comparable sites, 0 differences)
"ACG-ACGT" vs "ACGTACGT" → Hamming = 0

// Transition (A↔G = purine↔purine)
"ACGT" vs "GCGT" → 1 transition at pos 0

// Transversion (A↔C = purine↔pyrimidine)  
"ACGT" vs "CCGT" → 1 transversion at pos 0
```

---

## References
- Evidence Document: [PHYLO-DIST-001-Evidence.md](../docs/Evidence/PHYLO-DIST-001-Evidence.md)
- Algorithm Doc: [docs/algorithms/Phylogenetics/Distance_Matrix.md](../docs/algorithms/Phylogenetics/Distance_Matrix.md)
