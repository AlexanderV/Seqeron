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
| C1 | OptimizeSequence_HarmonizeExpression_MaintainsDistribution | Codon distribution matches host pattern | DNA Chisel `match_codon_usage` |

### Reference cross-checks (added review 2026-09)

| # | Test Name | Description | Evidence |
|---|-----------|-------------|----------|
| R1 | OptimizeSequence_MaximizeCAI_MatchesDnaChiselUseBestCodon | MaximizeCAI = DNA Chisel `CodonOptimize(method="use_best_codon")`, codon for codon | DNA Chisel 3.2.16 |
| R2 | OptimizeSequence_HarmonizeExpression_MatchesTargetUsageDeterministically | Deterministic largest-remainder allocation; same `match_codon_usage` score as DNA Chisel's optimizer (−11.77) | DNA Chisel 3.2.16 |
| R3 | OptimizeSequence_BalancedOptimization_MakesNoNeutralSwaps | GC pass only swaps codons that move GC toward the window and stops on entry | DNA Chisel `EnforceGCContent` |
| R4 | OptimizeSequence_AvoidRareCodons_NoSynonymAboveThreshold_UsesBestCodon | With no synonym above the threshold the best codon is still used | DNA Chisel `use_best_codon` |
| R5 | OptimizeSequence_AmbiguousCodons_HandledPerGeneticCode | GCN → best Ala codon; NNN and non-IUPAC triplets untouched, protein 'X' | GeneticCode / Biopython ambiguity |
| R6 | RemoveRestrictionSites_* (5 tests) | One codon per site, best-frequency substitution, both strands, IUPAC sites, unremovable site kept | REBASE/IUPAC, DNA Chisel `AvoidPattern` |
| R7 | CreateCodonTableFromSequence_MatchesBiopythonRelativeAdaptiveness | w from the built table = Biopython `CodonAdaptationIndex` (0.5 pseudo-count) | Biopython 1.88, Sharp & Li 1987 |

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
- **CAI formula**: Matches Sharp & Li (1987) definition: w_i = f_i / max(f_j), CAI = exp((1/L)·Σ ln(w_i)). Met/Trp and stops are not scored; w < 0.0001 → 0.01 (CodonW `cai_out`) — review 2026-09, CODON-CAI-001 F12/F13 (formerly Met/Trp scored with w = 1 and a 1e-6 clamp).
- **Standard genetic code**: All 64 codons verified correct.
- **Optimization strategies**: All thresholds exposed as configurable parameters (`rareCodonThreshold`, `gcTargetMin`, `gcTargetMax`); no hardcoded assumptions.
- **MinimizeSecondary**: Uses the same codon selection (and GC pass) as BalancedOptimization; the dedicated `ReduceSecondaryStructure` method handles structure reduction.

### Review 2026-09 (CODON-OPT-001, F21–F25)

- **Genetic code**: the private RNA-keyed `StandardGeneticCode` / `AminoAcidToCodons` copies were removed; amino acids and synonymous families now come from `GeneticCode.Standard` (NCBI table 1). Codon families are in NCBI order, which is the documented, deterministic tie-break for every "most frequent synonymous codon" choice. IUPAC-ambiguous triplets resolve through `GeneticCode.Translate` (GCN → Ala → best Ala codon); a triplet with non-IUPAC symbols is left untouched and contributes `X` to the protein.
- **Strategies**: `MaximizeCAI` = DNA Chisel `use_best_codon` (verified identical output). `AvoidRareCodeons` / `BalancedOptimization` no longer keep a rare codon when *no* synonym reaches the threshold — the most frequent synonym is used. `HarmonizeExpression` is no longer weighted-random (`new Random()` per call): it is the deterministic largest-remainder match of the target codon-usage profile (DNA Chisel `match_codon_usage`), which reaches the same DNA Chisel objective score as that library's own optimizer.
- **GC balancing**: only swaps that move GC toward `[gcTargetMin, gcTargetMax]` are applied (no neutral swaps), a candidate that lands inside the window is preferred over one that overshoots, the direction is re-evaluated each step and the pass stops as soon as the sequence is inside the window; the replacement frequency floor is the caller's `rareCodonThreshold` instead of a hard-coded 0.1; the GC count is maintained incrementally (was an O(n²) full rescan).
- **RemoveRestrictionSites**: exactly one codon changes per removed occurrence (the old loop kept rewriting the following codons after the site was already gone), the substitution is the most frequent synonymous codon in the supplied table (the table was previously ignored), sites are matched with IUPAC semantics and on both strands, and an unremovable occurrence no longer aborts the remaining occurrences.
- **ReduceSecondaryStructure**: base pairing delegates to `RnaSecondaryStructure.CanPair` (ViennaRNA pair set, so the G·U wobble counts); the window baseline is re-evaluated after each accepted change; output is always normalised RNA.
- **CreateCodonTableFromSequence**: codons absent from the reference set get the Sharp & Li (1987) count of 0.5 ("following the description in the original paper", Biopython `CodonAdaptationIndex` 1.88), so all 64 codons are present and the derived relative adaptiveness equals Biopython's index exactly. Previously an absent codon was missing from the table and was scored by `CalculateCAI` with the CodonW zero substitute 0.01.

## Date
2026-03-10 (reviewed 2026-09-28, campaign review-2026-09 / B02)

