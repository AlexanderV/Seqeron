# Test Specification: POP-FST-001 (F-Statistics)

**Test Unit ID:** POP-FST-001
**Algorithm:** F-Statistics (Fst, Fis, Fit)
**Class:** PopulationGeneticsAnalyzer
**Date:** 2026-02-01
**Status:** Ready for Implementation

---

## 1. Scope

### 1.1 Methods Under Test

| Method | Type | Test Depth |
|--------|------|------------|
| `CalculateFst(pop1, pop2)` | Canonical | Deep |
| `CalculateFStatistics(pop1Name, pop2Name, data)` | Canonical | Deep |
| `CalculatePairwiseFst(populations)` | Derived | Moderate |

### 1.2 Out of Scope

- Selection scanning (ScanForSelection uses Fst but is separate)
- Individual-level inbreeding coefficients (ROH-based)

---

## 2. Test Categories

### 2.1 Must Tests (Evidence-Based)

| ID | Test Case | Expected Result | Evidence |
|----|-----------|-----------------|----------|
| M1 | Identical populations | Fst = 0 | Wikipedia: "zero implies complete panmixia" |
| M2 | Fixed differences (p1=1.0, p2=0.0) | Fst = 1.0 | Wikipedia: "value of one implies... complete differentiation"; math proof: pBar=0.5, var=0.25, het=0.25 |
| M3 | Value range invariant | 0 ≤ Fst ≤ 1 | Wikipedia: "values range from 0 to 1" |
| M4 | Different populations | Fst > 0 | Mathematical definition |
| M5 | Empty populations | Fst = 0 (graceful) | Implementation contract |
| M6 | Pairwise matrix diagonal | All diagonal = 0 | Mathematical property |
| M7 | Pairwise matrix symmetry | Fst[i,j] = Fst[j,i] | Mathematical property |
| M8 | F-statistics components returned | Fis, Fit, Fst present | API contract |

### 2.2 Should Tests

| ID | Test Case | Expected Result | Rationale |
|----|-----------|-----------------|-----------|
| S1 | Moderate differentiation scenario | 0.05 < Fst < 0.25 | Realistic use case |
| S2 | Weighted by sample size | Larger samples have more weight | Wright (1965) - variance weighted by subpopulation sizes |
| S3 | Multiple loci aggregation | Combines across loci | Algorithm design |
| S4 | F-statistics partition relation | (1-Fit) = (1-Fis)(1-Fst) exactly | Wright's algebraic identity: (Hi/Hs)(Hs/Ht) = Hi/Ht |
| S5 | Three populations matrix | 3x3 symmetric matrix | Multi-population use |

### 2.3 Could Tests

| ID | Test Case | Expected Result | Rationale |
|----|-----------|-----------------|-----------|
| C1 | Reference population values | Match literature | Validation against Cavalli-Sforza |
| C2 | Large number of loci | Performance acceptable | Scalability |
| C3 | Single locus calculation | Valid result | Minimum input |

---

## 3. Test Data

### 3.1 Identical Populations (M1)
```csharp
var pop1 = [(0.5, 100), (0.3, 100)];
var pop2 = [(0.5, 100), (0.3, 100)];
// Expected: Fst = 0
```

### 3.2 Fixed Differences (M2)
```csharp
var pop1 = [(1.0, 100)];
var pop2 = [(0.0, 100)];
// Expected: Fst = 1.0 (Wright's variance Fst: pBar=0.5, var=0.25, het=0.25)
```

### 3.3 Moderate Differentiation (M4)
```csharp
var pop1 = [(0.9, 100), (0.8, 100)];
var pop2 = [(0.1, 100), (0.2, 100)];
// Expected: Fst > 0
```

### 3.4 F-Statistics Data (M8)
```csharp
var data = [
    (HetObs1: 20, N1: 50, HetObs2: 25, N2: 50, AlleleFreq1: 0.4, AlleleFreq2: 0.5),
    (HetObs1: 30, N1: 50, HetObs2: 15, N2: 50, AlleleFreq1: 0.5, AlleleFreq2: 0.3)
];
// Expected: FStatistics with Fis, Fit, Fst values
```

### 3.5 Three Populations (S5)
```csharp
var populations = [
    ("Pop1", [(0.5, 100)]),
    ("Pop2", [(0.6, 100)]),
    ("Pop3", [(0.9, 100)])
];
// Expected: 3x3 symmetric matrix with diagonal = 0
```

---

## 4. Invariants to Assert

### 4.1 Value Range
```csharp
Assert.That(fst, Is.InRange(0.0, 1.0));
```

### 4.2 Matrix Properties
```csharp
Assert.Multiple(() => {
    Assert.That(matrix[0, 0], Is.EqualTo(0));  // Diagonal
    Assert.That(matrix[0, 1], Is.EqualTo(matrix[1, 0]));  // Symmetric
});
```

### 4.3 Non-negative for Different Populations
```csharp
Assert.That(fst, Is.GreaterThanOrEqualTo(0));
```

---

## 5. Implementation

### 5.1 Algorithm

`CalculateFst` implements Wright's variance-based Fst (Wright 1965):

$$F_{ST} = \frac{\sigma_S^2}{\bar{p}(1-\bar{p})}$$

Where $\sigma_S^2 = \sum c_i (p_i - \bar{p})^2$ with $c_i = n_i / N$ (population-size weights).

Multi-locus Fst uses the ratio-of-sums estimator: $F_{ST} = \sum_l \sigma_{S,l}^2 / \sum_l \bar{p}_l(1-\bar{p}_l)$.

`CalculateFStatistics` uses heterozygosity-based definitions (Wikipedia F-statistics §Definitions):
- $F_{IS} = 1 - H_I/H_S$
- $F_{IT} = 1 - H_I/H_T$
- $F_{ST} = 1 - H_S/H_T$

The partition $(1-F_{IT}) = (1-F_{IS})(1-F_{ST})$ is an algebraic identity for this computation.

---

## 6. Test Structure

```
PopulationGeneticsAnalyzer_FStatistics_Tests.cs
├── CalculateFst Tests
│   ├── CalculateFst_IdenticalPopulations_ReturnsZero
│   ├── CalculateFst_FixedDifferences_ReturnsOne
│   ├── CalculateFst_ValueRange_BetweenZeroAndOne
│   ├── CalculateFst_EmptyPopulations_ReturnsZero
│   ├── CalculateFst_SingleLocus_ExactValue
│   ├── CalculateFst_MultiLocus_ExactValue
│   └── CalculateFst_UnequalSampleSizes_WeightedCalculation
├── CalculatePairwiseFst Tests
│   ├── CalculatePairwiseFst_ThreePopulations_CorrectDimensions
│   ├── CalculatePairwiseFst_DiagonalIsZero
│   ├── CalculatePairwiseFst_SymmetricMatrix
│   ├── CalculatePairwiseFst_AllValuesInRange
│   └── CalculatePairwiseFst_ExactCellValues
├── CalculateFStatistics Tests
│   ├── CalculateFStatistics_ReturnsAllComponents
│   ├── CalculateFStatistics_PopulationNamesPreserved
│   ├── CalculateFStatistics_EmptyData_ReturnsZeroValues
│   ├── CalculateFStatistics_ComponentsInValidRange
│   ├── CalculateFStatistics_PartitionRelationship_ExactIdentity
│   └── CalculateFStatistics_HandCalculated_ExactValues
├── Reference Data Validation Tests
│   ├── CalculateFst_MultiLocusModerate_ExactValue
│   ├── CalculateFst_WrightInterpretationScale_ExactValues
│   ├── CalculateFst_FixedAlleles_MultiLoci_ReturnsOne
│   └── CalculateFst_IslandModelConsistency_ExactValuesAndMonotonic
└── Missing Coverage Tests
    ├── CalculateFst_MonomorphicSites_ReturnsZero
    ├── CalculateFst_BothFixedSameAllele_ReturnsZero
    ├── CalculatePairwiseFst_ExactCellValues
    └── CalculateFStatistics_ExcessHeterozygosity_NegativeFis
```

---

## 7. Decisions

| Decision | Rationale |
|----------|-----------|
| Wright's variance Fst (not Weir-Cockerham θ) | Direct implementation of the theoretical definition; no ANOVA components needed for known allele frequencies |
| Partition formula tested as exact identity | $(H_I/H_S)(H_S/H_T) = H_I/H_T$ — algebraic, not approximate |
| Hand-calculated verification tests | Ensures implementation matches formula, not just "green tests" |
| Assert Fst = 1.0 for fixed differences | Mathematical certainty: pBar=0.5, var=0.25, het=0.25, ratio=1.0 || Exact binary fractions for FP-sensitive tests | Avoids IEEE 754 representation errors in hand-calculated expected values |
| Excess heterozygosity test (Fis < 0) | Verifies negative Fis case documented in Wikipedia |
| Removed duplicate DifferentPopulations test | Subsumed by MultiLocus_ExactValue with same data and exact assertion |
---

## 8. Deviations and Assumptions

None.
