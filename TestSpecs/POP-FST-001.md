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
| M2 | Fixed differences (p1=1.0, p2=0.0) | Fst > 0.5 (high) | Wikipedia: "value of one implies... complete differentiation" |
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
| S2 | Weighted by sample size | Larger samples have more weight | Weir-Cockerham |
| S3 | Multiple loci aggregation | Combines across loci | Algorithm design |
| S4 | F-statistics partition relation | (1-Fit) ≈ (1-Fis)(1-Fst) | Wright's formula |
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
// Expected: Fst > 0.5 (high differentiation)
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

## 5. Audit of Existing Tests

### 5.1 Current Tests in PopulationGeneticsAnalyzerTests.cs

| Test | Coverage | Status | Action |
|------|----------|--------|--------|
| `CalculateFst_IdenticalPopulations_ReturnsZero` | M1 | ✓ Adequate | Keep, move to canonical file |
| `CalculateFst_DifferentPopulations_ReturnsPositive` | M4 | ✓ Adequate | Keep, move to canonical file |
| `CalculateFst_FixedDifferences_ReturnsHighFst` | M2 | ✓ Adequate | Keep, move to canonical file |
| `CalculatePairwiseFst_ThreePopulations_ReturnsMatrix` | S5, M6, M7 | ✓ Adequate | Keep, move to canonical file |
| `CalculateFStatistics_ReturnsAllComponents` | M8 | Weak | Strengthen assertions |
| `CalculateFst_EmptyPopulations_ReturnsZero` | M5 | ✓ Adequate | Keep, move to canonical file |

### 5.2 Missing Tests

| Test | Priority | Reason |
|------|----------|--------|
| Value range invariant check | Must | M3 not explicitly tested |
| F-statistics partition verification | Should | S4 not tested |
| Fis/Fit range validation | Should | Component ranges not tested |
| Single locus edge case | Could | C3 not tested |

### 5.3 Consolidation Plan

1. **Create:** `PopulationGeneticsAnalyzer_FStatistics_Tests.cs` as canonical test file
2. **Move:** All Fst-related tests from `PopulationGeneticsAnalyzerTests.cs`
3. **Add:** Missing Must tests (M3)
4. **Strengthen:** F-statistics component tests
5. **Remove:** Duplicate coverage after consolidation

---

## 6. Test Structure

```
PopulationGeneticsAnalyzer_FStatistics_Tests.cs
├── CalculateFst_Tests (region)
│   ├── CalculateFst_IdenticalPopulations_ReturnsZero
│   ├── CalculateFst_DifferentPopulations_ReturnsPositive
│   ├── CalculateFst_FixedDifferences_ReturnsHighFst
│   ├── CalculateFst_ValueRange_BetweenZeroAndOne
│   ├── CalculateFst_EmptyPopulations_ReturnsZero
│   └── CalculateFst_SingleLocus_ValidResult
├── CalculatePairwiseFst_Tests (region)
│   ├── CalculatePairwiseFst_ThreePopulations_CorrectDimensions
│   ├── CalculatePairwiseFst_DiagonalIsZero
│   └── CalculatePairwiseFst_SymmetricMatrix
└── CalculateFStatistics_Tests (region)
    ├── CalculateFStatistics_ReturnsAllComponents
    ├── CalculateFStatistics_PopulationNamesPreserved
    ├── CalculateFStatistics_EmptyData_ReturnsZeroValues
    └── CalculateFStatistics_ComponentsInValidRange
```

---

## 7. Decisions

| Decision | Rationale |
|----------|-----------|
| Create separate canonical test file | Follows project convention (POP-FREQ-001, POP-DIV-001, POP-HW-001) |
| Move existing tests, don't duplicate | Test pool policy |
| Remove tests from PopulationGeneticsAnalyzerTests.cs | Eliminate duplicates after migration |
| Use Assert.Multiple for invariants | Group related assertions |

---

## 8. Open Questions

None - specification complete.
