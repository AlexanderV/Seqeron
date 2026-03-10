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
    double gcTargetMax = 0.60,
    double rareCodonThreshold = 0.15)
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

# E. coli rare codons for Arginine (Kazusa species=316407, W3110 K-12 substrain)
AGA, AGG → Frequency ~0.04, 0.02

# E. coli preferred codons for Leucine
CUG → Frequency 0.50

# Yeast preferred codons for Leucine (Kazusa species=4932)
UUA, UUG → Frequency 0.28, 0.29
```

### Known CAI Values (hand-verified against Sharp & Li formula)

- AUG alone: CAI = 1.0 (only codon for Met, wi = 1.0)
- CUGCCGACC (L-P-T, all optimal E. coli): CAI = 1.0
- CUAAGACGA (L-R-R, rare E. coli): CAI ≈ 0.106

## Audit Notes

### Existing Test Coverage Analysis

File: `CodonOptimizer_OptimizeSequence_Tests.cs`

| Region | Tests | Coverage | Assessment |
|--------|-------|----------|------------|
| Protein Preservation | 4+ | Good | All 5 strategies tested (incl. MinimizeSecondary) |
| CAI Behavior | 5 | Good | Covers increase/maintain, exact formula verification |
| Special Codons | 5 | Good | Met, Trp, Stop codons |
| Organism Specificity | 2 | Good | E. coli vs Yeast vs Human, exact optimized sequences per Kazusa |
| Input Handling | 2 | Good | Lowercase, exact result field values (hand-computed CAI, GC) |
| Strategy Specifics | 5 | Good | All strategies covered, AvoidRareCodons asserts replacement |
| Invariants | 6 | Good | Protein, CAI range, Sharp & Li formula, MaximizeCAI→1.0 |

### Coverage Classification Result (2026-03-10)

| Classification | Count | Details |
|---------------|-------|---------|
| ❌ Missing → Added | 1 | M1: MinimizeSecondary added to AllStrategies test |
| ⚠ Weak → Strengthened | 5 | M8: exact codon assertions per Kazusa; M10: exact hand-computed values; S2: GC enters target range; C1: CAI validity |
| 🔁 Duplicate → Removed | 19 | All OptimizeSequence duplicates in CodonOptimizerTests.cs |
| ✅ Covered | 24 | All remaining tests |

### Code Bug Fixed

- **BalancedOptimization Changes list**: `Changes` and `ChangedCodons` were not updated after GC content balancing phase. Fixed by rebuilding changes from original vs final codons.

### Data Source Traceability

| Organism | Kazusa Species ID | Kazusa Name | Dataset Size |
|----------|------------------|-------------|--------------|
| E. coli K12 | 316407 | E. coli W3110 (K-12 substrain) | 4332 CDS |
| S. cerevisiae | 4932 | Saccharomyces cerevisiae | 14411 CDS |
| H. sapiens | 9606 | Homo sapiens | 93487 CDS |

---

## Deviations and Assumptions

- **BalancedOptimization Changes rebuild (fixed 2026-03-10)**: Previously, `Changes` list only reflected the initial optimization pass, missing GC content balancing modifications. Fixed to rebuild changes by comparing original vs final codons.
- **Codon usage tables**: All three tables (E. coli, Yeast, Human) verified against Kazusa Codon Usage Database raw data (per-thousand frequencies → relative fractions per amino acid).
- **CAI formula**: Matches Sharp & Li (1987) definition: w_i = f_i / max(f_j), CAI = exp((1/L)·Σ ln(w_i)). Zero-frequency codons clamped to 1e-6 per original prescription.
- **Standard genetic code**: All 64 codons verified correct.
- **Optimization strategies**: All thresholds exposed as configurable parameters (`rareCodonThreshold`, `gcTargetMin`, `gcTargetMax`); no hardcoded assumptions.
- **MinimizeSecondary**: Falls through to BalancedOptimization in `SelectOptimalCodon`; separate `ReduceSecondaryStructure` method exists for dedicated secondary structure reduction.

## Date
2026-03-10

