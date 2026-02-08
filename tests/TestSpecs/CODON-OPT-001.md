# Test Specification: CODON-OPT-001

## Test Unit Information
- **ID:** CODON-OPT-001
- **Title:** Sequence Optimization
- **Canonical Method:** `CodonOptimizer.OptimizeSequence(...)`
- **Area:** Codon Optimization
- **Complexity:** O(n)
- **Status:** ☑ Complete

## Method Under Test

```csharp
public static OptimizationResult OptimizeSequence(
    string codingSequence,
    CodonUsageTable targetOrganism,
    OptimizationStrategy strategy = OptimizationStrategy.BalancedOptimization,
    double gcTargetMin = 0.40,
    double gcTargetMax = 0.60)
```

## Test Categories

### MUST Tests (Required for Completion)

| # | Test Name | Description | Evidence |
|---|-----------|-------------|----------|
| M1 | OptimizeSequence_PreservesProtein_AllStrategies | Optimized sequence must encode identical protein | Synonymous substitution definition |
| M2 | OptimizeSequence_EmptySequence_ReturnsEmptyResult | Empty input returns empty output with CAI=0 | Edge case definition |
| M3 | OptimizeSequence_ConvertsThymine_ToUracil | T is converted to U in RNA representation | RNA notation standard |
| M4 | OptimizeSequence_TrimsToCompleteCodons | Incomplete codons (length % 3 != 0) are trimmed | Codon definition |
| M5 | OptimizeSequence_MaximizeCAI_IncreasesOrMaintainsCAI | CAI after optimization >= CAI before | Sharp & Li (1987) |
| M6 | OptimizeSequence_SingleAminoAcidCodons_Unchanged | AUG (Met) and UGG (Trp) cannot be changed | Standard genetic code |
| M7 | OptimizeSequence_StopCodons_Preserved | Stop codons remain as stop codons | Translation termination |
| M8 | OptimizeSequence_DifferentOrganisms_DifferentResults | E. coli vs Yeast optimization differs | Organism-specific bias |
| M9 | OptimizeSequence_LowercaseInput_Handled | Case-insensitive processing | Robustness requirement |
| M10 | OptimizeSequence_ReturnsValidOptimizationResult | All result fields populated correctly | API contract |

### SHOULD Tests (Recommended)

| # | Test Name | Description | Evidence |
|---|-----------|-------------|----------|
| S1 | OptimizeSequence_AvoidRareCodons_OnlyReplacesRare | Only codons below threshold are changed | Strategy definition |
| S2 | OptimizeSequence_BalancedOptimization_AimsForTargetGc | GC content moves toward 40-60% range | Implementation spec |
| S3 | OptimizeSequence_TracksChanges_Correctly | Changes list accurately reflects modifications | API contract |
| S4 | OptimizeSequence_LongSequence_Completes | Performance acceptable for long sequences | Usability |
| S5 | OptimizeSequence_KnownSequence_VerifiedOutput | GFP or similar known sequence optimizes correctly | Integration test |

### COULD Tests (Optional)

| # | Test Name | Description | Evidence |
|---|-----------|-------------|----------|
| C1 | OptimizeSequence_HarmonizeExpression_MaintainsDistribution | Codon distribution matches host pattern | Mignon et al. (2018) |

## Invariants to Verify

1. **Protein Preservation**: `original.ProteinSequence == optimized.ProteinSequence`
2. **CAI Improvement**: `MaximizeCAI strategy → OptimizedCAI >= OriginalCAI`
3. **Codon Length**: `OptimizedSequence.Length % 3 == 0`
4. **Same Length**: `OriginalSequence.Length == OptimizedSequence.Length`

## Edge Cases

| Case | Input | Expected Output |
|------|-------|-----------------|
| Empty | `""` | Empty result, CAI=0 |
| Single codon (Met) | `"AUG"` | Unchanged |
| Single codon (Trp) | `"UGG"` | Unchanged |
| DNA input | `"ATGGCT"` | Converts T→U, then optimizes |
| Incomplete codon | `"AUGGCUA"` (7 nt) | Trims to `"AUGGCU"` (6 nt) |
| All optimal | Already optimal E. coli | No changes, high CAI |
| All rare | All rare codons | Significant changes, CAI improves |

## Test Data

### Reference Sequences

```
# Short test sequence (M-A-Stop)
AUGGCUUAA → Protein: MA*

# E. coli rare codons for Arginine
AGA, AGG → Frequency ~0.07, 0.04

# E. coli preferred codons for Leucine
CUG → Frequency 0.47

# Yeast preferred codons for Leucine
UUA, UUG → Frequency 0.28, 0.29
```

### Known CAI Values (ASSUMPTION)

- AUG alone: CAI = 1.0 (only codon for Met)
- All optimal E. coli codons: CAI > 0.8
- All rare E. coli codons: CAI < 0.3

## Audit Notes

### Existing Test Coverage Analysis

File: `CodonOptimizerTests.cs` (481 lines)

| Region | Tests | Coverage | Assessment |
|--------|-------|----------|------------|
| Standard Genetic Code | 4 | Good | Covers basic optimization |
| CAI Calculation | 5 | Good | Covers CAI scenarios |
| Optimization Strategy | 5 | Partial | Missing some strategy-specific tests |
| GC Content | 2 | Weak | Limited GC testing |
| Codon Usage Analysis | 3+ | Good | Covers usage functions |
| Edge Cases | 7+ | Good | Good edge case coverage |
| Integration | 2 | Adequate | Full workflow covered |

### Consolidation Plan

1. **Keep existing tests**: Well-structured with good coverage
2. ~~**Add invariant tests**: Assert.Multiple for protein preservation~~ ✅ Done
3. ~~**Add strategy-specific MUST tests**: Ensure all strategies are tested~~ ✅ Done
4. ~~**Strengthen CAI mathematical tests**: Verify formula correctness~~ ✅ Done
5. **Remove duplicates**: None identified

## Open Questions

None - existing implementation and tests are well-structured.

## Decisions

1. Tests for CODON-CAI-001, CODON-RARE-001, CODON-USAGE-001 exist in same file
   - These will be separated in their respective Test Units
   - CODON-OPT-001 focuses only on OptimizeSequence method

## Date
2026-02-04
