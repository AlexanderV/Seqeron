# Test Specification: CODON-CAI-001

## Test Unit Information
- **ID:** CODON-CAI-001
- **Title:** Codon Adaptation Index (CAI) Calculation
- **Canonical Method:** `CodonOptimizer.CalculateCAI(string, CodonUsageTable)`
- **Area:** Codon Optimization
- **Complexity:** O(n)

## Method Under Test

```csharp
public static double CalculateCAI(string codingSequence, CodonUsageTable table)
```

## Algorithm Summary

CAI = geometric mean of relative adaptiveness values:
- `w_i = f_i / max(f_j)` for synonymous codons
- `CAI = exp((1/L) × Σ ln(w_i))`

## Test Categories

### MUST Tests (Required for Completion)

| # | Test Name | Description | Evidence |
|---|-----------|-------------|----------|
| M1 | CalculateCAI_EmptySequence_ReturnsZero | Empty input returns 0 | Edge case convention |
| M2 | CalculateCAI_SingleMetCodon_ReturnsOne | AUG (only Met codon) has CAI=1.0 | Mathematical: w=1.0/1.0=1.0 |
| M3 | CalculateCAI_SingleTrpCodon_ReturnsOne | UGG (only Trp codon) has CAI=1.0 | Mathematical: w=1.0/1.0=1.0 |
| M4 | CalculateCAI_AllOptimalCodons_ReturnsOne | All optimal codons → CAI=1.0 | Sharp & Li (1987) |
| M5 | CalculateCAI_RareCodons_ReturnsLow | Rare codons → CAI < 0.5 | Sharp & Li (1987) |
| M6 | CalculateCAI_RangeIsZeroToOne | CAI always in [0, 1] | CAI definition |
| M7 | CalculateCAI_DifferentOrganisms_DifferentResults | Same sequence, different CAI per organism | Organism-specific bias |
| M8 | CalculateCAI_DnaInput_HandledCorrectly | T→U conversion works | Implementation requirement |
| M9 | CalculateCAI_LowercaseInput_Handled | Case-insensitive processing | Robustness |
| M10 | CalculateCAI_ExcludesStopCodons | Stop codons not counted in calculation | Standard practice |
| M11 | CalculateCAI_GeometricMeanProperty | Single rare codon significantly lowers CAI | Mathematical property |
| M12 | CalculateCAI_HandCalculatedValue_Matches | Verify against hand-calculated example | Validation |

### SHOULD Tests (Recommended)

| # | Test Name | Description | Evidence |
|---|-----------|-------------|----------|
| S1 | CalculateCAI_LongSequence_Completes | Performance for long sequences | Usability |
| S2 | CalculateCAI_MixedOptimalRare_IntermediateValue | Mix yields intermediate CAI | Mathematical property |
| S3 | CalculateCAI_OnlyStopCodons_ReturnsZero | Sequence of only stops → 0 | Edge case |

### COULD Tests (Optional)

| # | Test Name | Description | Evidence |
|---|-----------|-------------|----------|
| C1 | CalculateCAI_AllThreeOrganismTables_Valid | E. coli, Yeast, Human all work | API coverage |

## Invariants to Verify

1. **Range Invariant:** `0 ≤ CAI ≤ 1`
2. **Single-Codon AA:** Met/Trp codons contribute w=1.0 always
3. **Monotonicity:** Replacing rare codons with optimal ones increases CAI
4. **Idempotence:** Calculating twice gives same result

## Edge Cases

| Case | Input | Expected Output |
|------|-------|-----------------|
| Empty | `""` | 0 |
| Single Met | `"AUG"` | 1.0 |
| Single Trp | `"UGG"` | 1.0 |
| DNA format | `"ATGCTG"` | Same as `"AUGCUG"` |
| Lowercase | `"augcug"` | Same as `"AUGCUG"` |
| Incomplete codon | `"AUGC"` | Based on `"AUG"` only (1.0) |
| Only stop codons | `"UAAUAGUGA"` | 0 (no codons to evaluate) |

## Test Data

### Hand-Calculated Reference Values (E. coli K12)

**Test 1: All Optimal Codons**
```
Sequence: CUGCCGACC (Leu-Pro-Thr)
Codons:   CUG(0.47), CCG(0.49), ACC(0.40)
w values: 0.47/0.47=1.0, 0.49/0.49=1.0, 0.40/0.40=1.0
CAI = (1.0 × 1.0 × 1.0)^(1/3) = 1.0
```

**Test 2: Mixed Codons**
```
Sequence: AUGCUGACC (Met-Leu-Thr)
AUG: 1.0/1.0 = 1.0
CUG: 0.47/0.47 = 1.0
ACC: 0.40/0.40 = 1.0
CAI = (1.0 × 1.0 × 1.0)^(1/3) = 1.0
```

**Test 3: Rare Codons**
```
Sequence: CUAACU (Leu-Thr)
CUA: 0.04/0.47 ≈ 0.085
ACU: 0.19/0.40 ≈ 0.475
CAI = (0.085 × 0.475)^(1/2) ≈ 0.20
```

### Codon Frequencies (E. coli K12)

| AA | Codon | Freq | Max | w |
|----|-------|------|-----|---|
| Leu | CUG | 0.47 | 0.47 | 1.00 |
| Leu | CUA | 0.04 | 0.47 | 0.085 |
| Pro | CCG | 0.49 | 0.49 | 1.00 |
| Pro | CCA | 0.20 | 0.49 | 0.408 |
| Thr | ACC | 0.40 | 0.40 | 1.00 |
| Thr | ACU | 0.19 | 0.40 | 0.475 |

## Audit Notes

### Existing Test Coverage

**CodonOptimizerTests.cs (lines 50-90):**
- `CalculateCAI_AllOptimalCodons_HighCAI` - Partial coverage
- `CalculateCAI_RareCodons_LowerCAI` - Partial coverage
- `CalculateCAI_EmptySequence_ReturnsZero` - Covered
- `CalculateCAI_SingleCodon_ReturnsValue` - Covered
- `CalculateCAI_DifferentOrganisms_DifferentResults` - Covered

**CodonUsageAnalyzerTests.cs (lines 113-148):**
- Tests for `CodonUsageAnalyzer.CalculateCai` (different class)
- Duplicate/related tests for a wrapper method

### Consolidation Plan

1. **Keep in canonical file:** CodonOptimizer_CAI_Tests.cs (new file)
2. **Remove duplicates:** CodonOptimizerTests.cs CAI tests move to new file
3. **CodonUsageAnalyzer tests:** Keep as smoke tests (different API)

### Missing Coverage

- M4: Exact CAI=1.0 for all-optimal codons not verified precisely
- M10: Stop codon exclusion not explicitly tested
- M11: Geometric mean sensitivity not tested
- M12: No hand-calculated verification test

## Open Questions

None - algorithm well-documented in Sharp & Li (1987).

## Decisions

1. Create new canonical test file: `CodonOptimizer_CAI_Tests.cs`
2. Consolidate CAI tests from `CodonOptimizerTests.cs`
3. Keep `CodonUsageAnalyzerTests.cs` CAI tests as smoke verification (different class)
