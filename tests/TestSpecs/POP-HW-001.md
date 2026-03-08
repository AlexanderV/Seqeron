# Test Specification: POP-HW-001

**Test Unit ID:** POP-HW-001
**Title:** Hardy-Weinberg Equilibrium Test
**Algorithm Group:** Population Genetics
**Status:** Complete — Coverage Classification Done
**Last Updated:** 2026-03-08

---

## Scope

### Canonical Methods (Deep Testing Required)

| Method | Class | Complexity |
|--------|-------|------------|
| `TestHardyWeinberg(variantId, AA, Aa, aa, significance)` | PopulationGeneticsAnalyzer | O(1) |

### Associated Types

| Type | Description |
|------|-------------|
| `HardyWeinbergResult` | Record struct containing test results |

---

## Test Cases

### MUST Tests (Required - Evidence-Based)

#### Published Dataset Tests

| ID | Test Case | Evidence Source | Expected |
|----|-----------|-----------------|----------|
| HW-M01 | Ford's moth data (1469, 138, 5) | Wikipedia (Ford 1971) | χ² ≈ 0.83, InEquilibrium = true |
| HW-M02 | Perfect HWE with p=0.5 (25, 50, 25) | Mathematical definition | χ² ≈ 0, InEquilibrium = true |

#### Deviation Detection Tests

| ID | Test Case | Evidence Source | Expected |
|----|-----------|-----------------|----------|
| HW-M03 | Excess heterozygotes (10, 80, 10) | Wikipedia deviation example | χ² >> 3.84, InEquilibrium = false |
| HW-M04 | All heterozygotes (0, 100, 0) | Extreme deviation | χ² very high, InEquilibrium = false |
| HW-M05 | Deficit heterozygotes (45, 10, 45) | Inbreeding pattern | χ² > 3.84, InEquilibrium = false |

#### Expected Count Verification

| ID | Test Case | Evidence Source | Expected |
|----|-----------|-----------------|----------|
| HW-M06 | Fixed major allele (100, 0, 0) | Definition: p=1 | E(AA)=100, E(Aa)=0, E(aa)=0, InEquilibrium = true |
| HW-M07 | Fixed minor allele (0, 0, 100) | Definition: q=1 | E(AA)=0, E(Aa)=0, E(aa)=100, InEquilibrium = true |
| HW-M08 | Equal proportions p=0.5 (36, 48, 16) | Definition | E(AA)≈25, E(Aa)≈50, E(aa)≈25 |

#### Edge Cases

| ID | Test Case | Evidence Source | Expected |
|----|-----------|-----------------|----------|
| HW-M09 | Zero samples (0, 0, 0) | Implementation guard | InEquilibrium = true, PValue = 1 |
| HW-M10 | Single sample (1, 0, 0) | Minimum valid | Result returned, calculations valid |
| HW-M11 | VariantId preserved in result | API contract | result.VariantId == input |

#### Invariant Tests

| ID | Test Case | Evidence Source | Expected |
|----|-----------|-----------------|----------|
| HW-M12 | ChiSquare ≥ 0 | χ² distribution property | Always non-negative |
| HW-M13 | PValue ∈ [0, 1] | Probability range | 0 ≤ PValue ≤ 1 |
| HW-M14 | InEquilibrium consistent with PValue/significance | Definition | InEquilibrium ⟺ (PValue ≥ significance) |
| HW-M15 | Expected counts sum to n | Conservation | E(AA) + E(Aa) + E(aa) = n |

#### Significance Level Tests

| ID | Test Case | Evidence Source | Expected |
|----|-----------|-----------------|----------|
| HW-M16 | Custom significance (α=0.01) | API flexibility | Uses provided significance level |
| HW-M17 | Default significance (α=0.05) | Default parameter | Uses 0.05 when not specified |

### SHOULD Tests (Recommended)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| HW-S01 | Large sample (n=10000) exact match | Numerical stability |
| HW-S02 | Various p values (0.1, 0.3, 0.7, 0.9) | Coverage of frequency spectrum |
| HW-S03 | Borderline significance (AA=46, Aa=49, aa=5; χ²≈3.17, pval≈0.075) | Decision boundary: equilibrium at α=0.05, not at α=0.10 |

### COULD Tests (Optional)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| HW-C01 | Property: same result for symmetric genotypes | Symmetry under allele relabeling |
| HW-C02 | Extreme frequency (p=0.99) with large n | Rare variant scenario |

---

## Invariants to Verify

### Chi-Square Invariants

1. ∀ inputs: χ² ≥ 0
2. When observed = expected: χ² = 0
3. χ² = Σ((O-E)²/E) computed only for E > 0

### P-Value Invariants

1. ∀ inputs: 0 ≤ PValue ≤ 1
2. χ² = 0 → PValue = 1
3. χ² → ∞ → PValue → 0
4. PValue from chi-square CDF with df = 1

### Equilibrium Decision Invariants

1. InEquilibrium ⟺ (PValue ≥ significanceLevel)
2. Default significanceLevel = 0.05
3. Critical χ² at α=0.05, df=1 is 3.841

### Expected Count Invariants

1. E(AA) + E(Aa) + E(aa) = n (within floating-point tolerance)
2. E(AA) = p² × n
3. E(Aa) = 2pq × n
4. E(aa) = q² × n
5. p + q = 1

---

## Audit Results

### Canonical Test File: `PopulationGeneticsAnalyzer_HardyWeinberg_Tests.cs`

22 test methods, 45 test runs total (17 MUST, 3 SHOULD, 2 COULD, plus parameterized invariant cases).

### Supporting Test Files

| File | Tests | Role |
|------|-------|------|
| `PopulationGeneticsProperties.cs` | 2 | Property-based (FsCheck): χ²≥0, PValue∈[0,1] |
| `PopulationSnapshotTests.cs` | 1 | Snapshot regression guard |

---

## Coverage Classification (2026-03-08)

### Canonical (`PopulationGeneticsAnalyzer_HardyWeinberg_Tests.cs`) — 22 test methods

| # | Test Method | Spec ID | Status |
|---|-------------|---------|--------|
| 1 | `TestHardyWeinberg_FordsMothData_IsInEquilibrium` | HW-M01 | ✅ |
| 2 | `TestHardyWeinberg_PerfectEquilibrium_ChiSquareNearZero` | HW-M02 | ✅ |
| 3 | `TestHardyWeinberg_ExcessHeterozygotes_DeviatesFromEquilibrium` | HW-M03 | ✅ |
| 4 | `TestHardyWeinberg_AllHeterozygotes_SignificantDeviation` | HW-M04 | ✅ |
| 5 | `TestHardyWeinberg_DeficitHeterozygotes_InbreedingPattern` | HW-M05 | ✅ |
| 6 | `TestHardyWeinberg_FixedMajorAllele_ExpectedCountsCorrect` | HW-M06 | ✅ |
| 7 | `TestHardyWeinberg_FixedMinorAllele_ExpectedCountsCorrect` | HW-M07 | ✅ |
| 8 | `TestHardyWeinberg_VerifyExpectedCountFormulas` | HW-M08 | ✅ |
| 9 | `TestHardyWeinberg_ZeroSamples_ReturnsEquilibrium` | HW-M09 | ✅ |
| 10 | `TestHardyWeinberg_SingleSample_ReturnsValidResult` | HW-M10 | ✅ |
| 11 | `TestHardyWeinberg_VariantIdPreserved` | HW-M11 | ✅ |
| 12 | `TestHardyWeinberg_ChiSquareNonNegative` (×6 cases) | HW-M12 | ✅ |
| 13 | `TestHardyWeinberg_PValueInValidRange` (×5 cases) | HW-M13 | ✅ |
| 14 | `TestHardyWeinberg_EquilibriumConsistentWithPValue` (×5 cases) | HW-M14 | ✅ |
| 15 | `TestHardyWeinberg_ExpectedCountsSumToN` (×3 cases) | HW-M15 | ✅ |
| 16 | `TestHardyWeinberg_CustomSignificanceLevel` | HW-M16 | ✅ |
| 17 | `TestHardyWeinberg_DefaultSignificanceIs005` | HW-M17 | ✅ |
| 18 | `TestHardyWeinberg_LargeSample_NumericallyStable` | HW-S01 | ✅ |
| 19 | `TestHardyWeinberg_VariousAlleleFrequencies_CorrectExpected` (×4 cases) | HW-S02 | ✅ |
| 20 | `TestHardyWeinberg_BorderlineSignificance_DecisionDependsOnAlpha` | HW-S03 | ✅ |
| 21 | `TestHardyWeinberg_SymmetricGenotypes_SameChiSquareAndPValue` (×3 cases) | HW-C01 | ✅ |
| 22 | `TestHardyWeinberg_ExtremeFrequency_RareVariantStable` | HW-C02 | ✅ |

### Property Tests (`PopulationGeneticsProperties.cs`) — 2 tests

| # | Test Method | Type | Status |
|---|-------------|------|--------|
| 1 | `HardyWeinberg_ChiSquare_IsNonNegative` | FsCheck property | ✅ — complements HW-M12 with random inputs |
| 2 | `HardyWeinberg_PValue_InRange` | FsCheck property | ✅ — complements HW-M13 with random inputs |

### Snapshot Tests (`PopulationSnapshotTests.cs`) — 1 test

| # | Test Method | Type | Status |
|---|-------------|------|--------|
| 1 | `HardyWeinberg_KnownGenotypes_MatchesSnapshot` | Regression guard | ✅ — different purpose from functional tests |

### Classification Summary

- ✅ Covered: 22 canonical + 2 property + 1 snapshot = 25 total (45 test runs)
- ❌ Missing: 0
- ⚠ Weak: 0
- 🔁 Duplicate: 0

### Changes Applied (2026-03-08)

| Action | Test | Detail |
|--------|------|--------|
| ⚠→✅ | HW-M01 | Tightened χ² tolerance: Within(0.1)→Within(0.01); expected counts: Within(1.0)→Within(0.1) |
| ⚠→✅ | HW-M03 | Tightened χ² tolerance: Within(1.0)→Within(0.01) |
| ⚠→✅ | HW-M10 | Added PValue=1, InEquilibrium=true, all expected count assertions |
| ⚠→✅ | HW-M14 | Removed unused `expectedEquilibrium` parameter; added 2 borderline test cases |
| ⚠→✅ | HW-M16 | Replaced conditional assertion with guaranteed borderline data (46,49,5) |
| ⚠→✅ | HW-M17 | Replaced trivial data (25,50,25) with non-trivial borderline data (46,49,5) |
| ⚠→✅ | HW-S02 | Tightened tolerance: Within(1.0)→Within(0.01); added χ²=0 and InEquilibrium assertions |
| 🔁→∅ | Properties `HardyWeinberg_ExpectedFrequencies_SumToTotal` | Removed — duplicate of HW-M15 |
| ❌→✅ | HW-C01 | Added symmetry property test (3 parameterized cases) |
| ❌→✅ | HW-C02 | Added extreme frequency test (p=0.99, n=10000) |

---

## Deviations and Assumptions

None. All formulas, behaviors, and edge cases are sourced from external references:

- Wikipedia: Hardy–Weinberg principle (accessed 2026-02-01, verified 2026-03-08)
- Wikipedia: Chi-squared distribution — CDF table, critical values (accessed 2026-02-01, verified 2026-03-08)
- Ford (1971) Ecological Genetics — Scarlet tiger moth dataset (via Wikipedia)
- Hardy (1908) "Mendelian Proportions in a Mixed Population", Science 28(706):49–50
- Weinberg (1908) "Über den Nachweis der Vererbung beim Menschen", Jahreshefte 64:368–382
- Emigh (1980) "A Comparison of Tests for Hardy–Weinberg Equilibrium", Biometrics 36(4):627–642

### Edge Case Rationale (Not Assumptions)

| Case | Behavior | Justification |
|------|----------|---------------|
| n = 0 | χ²=0, PValue=1, InEquilibrium=true | No data = no evidence against H₀ = fail to reject. Consistent with hypothesis testing: PValue=1 means maximum compatibility with the null hypothesis. |
| E(genotype) = 0 | Skip that term in χ² sum | Standard practice per Wikipedia: "χ² = Σ((O−E)²/E)" is computed only for categories where E > 0, to avoid division by zero. |
| Fixed allele (monomorphic) | χ²=0, InEquilibrium=true | When p=1 (or q=1), all expected counts match observed. Observed = Expected ⇒ χ²=0. |

### Verified Numerical Values (Source: Wikipedia + Independent Python Verification)

| Quantity | Wikipedia Value | Exact Computation | Status |
|----------|----------------|-------------------|--------|
| Ford χ² | 0.83 (rounded) | 0.8309 | ✓ Test uses Within(0.01) |
| Ford E(AA) | 1467.4 | 1467.40 | ✓ Test uses Within(0.1) |
| Ford E(Aa) | 141.2 | 141.21 | ✓ Test uses Within(0.1) |
| Ford E(aa) | 3.4 | 3.40 | ✓ Test uses Within(0.1) |
| Critical χ² (α=0.05, df=1) | 3.84 (2 d.p.) | 3.8415 | ✓ Constant = 3.841 |
| Excess het χ² | 36 | 36.0 (exact) | ✓ Test uses Within(0.01) |
| All het χ² | — | 100.0 (exact) | ✓ Test uses Within(0.01) |
| Deficit het χ² | — | 64.0 (exact) | ✓ Test uses Within(0.01) |
| Borderline χ² | — | 3.1693 | ✓ Test asserts ∈ (2.706, 3.841) |
| df for biallelic HWE | 1 | 3 genotypes − 2 alleles = 1 | ✓ |

---

## Open Questions / Decisions

| # | Question | Decision |
|---|----------|----------|
| 1 | What to return for n=0? | Return equilibrium (χ²=0, p=1) — no evidence against H₀ means fail to reject. Not an assumption; follows from hypothesis testing framework. |
| 2 | Skip terms in χ² when E=0? | Yes — standard practice from Wikipedia formula definition. |

---

## Validation Checklist

- [x] All MUST tests have evidence source (Wikipedia, Ford 1971, mathematical definition)
- [x] All SHOULD tests covered with tightened tolerances
- [x] All COULD tests implemented (symmetry property, extreme frequency)
- [x] Invariants are mathematically verifiable from Wikipedia / chi-square theory
- [x] Edge cases documented (zero samples, single sample, fixed allele, monomorphic)
- [x] All assertions use exact values with tight tolerances (≤0.01 for χ², ≤0.1 for expected counts)
- [x] No conditional assertions — all paths guaranteed to execute
- [x] Cross-verified against independent Python computation (scipy-equivalent math.erf)
- [x] No assumptions — all behaviors sourced from Wikipedia or standard statistics
- [x] No duplicates — removed `HardyWeinberg_ExpectedFrequencies_SumToTotal` from Properties
- [x] Coverage classification complete: 0 missing, 0 weak, 0 duplicate
