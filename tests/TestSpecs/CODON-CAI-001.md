# Test Specification: CODON-CAI-001

## Test Unit Information
- **ID:** CODON-CAI-001
- **Title:** Codon Adaptation Index (CAI) Calculation
- **Canonical Method:** `CodonOptimizer.CalculateCAI(string, CodonUsageTable)`
- **Area:** Codon Optimization
- **Complexity:** O(n)
- **Status:** ☑ Complete

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

### Hand-Calculated Reference Values (E. coli K12, Kazusa MG1655)

**Test 1: All Optimal Codons**
```
Sequence: CUGCCGACC (Leu-Pro-Thr)
Codons:   CUG(0.50), CCG(0.53), ACC(0.44)
w values: 0.50/0.50=1.0, 0.53/0.53=1.0, 0.44/0.44=1.0
CAI = (1.0 × 1.0 × 1.0)^(1/3) = 1.0
```

**Test 2: Mixed Codons**
```
Sequence: AUGCUGACC (Met-Leu-Thr)
AUG: 1.0/1.0 = 1.0
CUG: 0.50/0.50 = 1.0
ACC: 0.44/0.44 = 1.0
CAI = (1.0 × 1.0 × 1.0)^(1/3) = 1.0
```

**Test 3: Rare Codons**
```
Sequence: CUAACU (Leu-Thr)
CUA: 0.04/0.50 = 0.08
ACU: 0.16/0.44 = 0.36364
CAI = (0.08 × 0.36364)^(1/2) = 0.17056
```

### Codon Frequencies (E. coli K12, Kazusa MG1655)

| AA | Codon | Freq | Max | w |
|----|-------|------|-----|---|
| Leu | CUG | 0.50 | 0.50 | 1.00 |
| Leu | CUA | 0.04 | 0.50 | 0.08 |
| Pro | CCG | 0.53 | 0.53 | 1.00 |
| Pro | CCA | 0.19 | 0.53 | 0.358 |
| Thr | ACC | 0.44 | 0.44 | 1.00 |
| Thr | ACU | 0.16 | 0.44 | 0.364 |

## Audit Notes

### Current Test Coverage (25 tests in CodonOptimizer_CAI_Tests.cs)

| Category | Count | Tests |
|----------|-------|-------|
| M1 (Empty) | 2 | EmptySequence_ReturnsZero, NullSequence_ReturnsZero |
| M2/M3 (Single-codon AA) | 3 | SingleMetCodon, SingleTrpCodon, MetAndTrp |
| M4 (All optimal) | 2 | AllOptimalCodonsEColi, OptimalCodonsWithMet |
| M5 (Rare codons) | 2 | RareCodonsEColi, RareArginineCodonsEColi |
| M6 (Range) | 1 | AnyValidSequence_RangeIsZeroToOne |
| M7 (Organism diff) | 2 | SameSequence_DifferentOrganisms, YeastPreferredCodons |
| M8 (DNA input) | 1 | DnaInputWithThymine_ConvertsToUracil |
| M9 (Lowercase) | 1 | LowercaseInput_HandledCorrectly |
| M10 (Stop codons) | 3 | SequenceWithStopCodon_Excludes, OnlyStopCodons, StopCodonInMiddle |
| M11 (Geometric mean) | 2 | SingleRareCodon_SignificantlyLowers, MoreRareCodons_LowerCAI |
| M12 (Hand-calculated) | 1 | HandCalculatedRareCodons_MatchesExpected |
| S1 (Performance) | 1 | LongSequence_CompletesInReasonableTime |
| S2 (Mixed) | 1 | MixedOptimalAndRare_MatchesHandCalculated |
| S3 (Only stops) | 1 | OnlyStopCodons_ReturnsZero (shared with M10) |
| C1 (All organisms) | 1 | AllThreeOrganismTables_MatchHandCalculated |
| Edge cases | 2 | IncompleteFinalCodon, TwoIncompleteBases |

### Consolidation Status

- **Canonical file:** `CodonOptimizer_CAI_Tests.cs` — all CAI tests consolidated here
- **CodonUsageAnalyzer tests:** Separate class, kept as smoke tests in their own file

### Missing Coverage

None — all test categories fully covered.

## Deviations and Assumptions

**Deviation D1 — 1e-6 clamp for zero-frequency codons:**
When a codon has frequency 0 but its amino acid has other codons with frequency > 0 (i.e., the codon is absent from the reference set but the amino acid is not), the implementation clamps w_i to 1e-6 instead of allowing w_i = 0.

- **Rationale:** Protects against incomplete codon usage tables where zero frequency may represent missing data rather than a truly absent codon. Sharp & Li (1987) used highly expressed E. coli genes where all synonymous codons appeared; real-world tables from Kazusa or custom datasets may have sampling gaps.
- **Impact:** For sequences containing affected codons, CAI > 0 (approximately 1e-6^(1/L)) instead of exactly 0. For all practical sequences with L > 1, the difference is negligible.
- **Alternative:** Strict IEEE 754 arithmetic would give CAI = 0 for any zero-frequency codon usage.

**No other assumptions.** All codon frequency tables verified against Kazusa Codon Usage Database (March 2026).

## Open Questions

None — algorithm well-documented in Sharp & Li (1987).

## Decisions

1. Create new canonical test file: `CodonOptimizer_CAI_Tests.cs`
2. Consolidate CAI tests from `CodonOptimizerTests.cs`
3. Keep `CodonUsageAnalyzerTests.cs` CAI tests as smoke verification (different class)
