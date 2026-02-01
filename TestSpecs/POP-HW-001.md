# Test Specification: POP-HW-001

**Test Unit ID:** POP-HW-001  
**Title:** Hardy-Weinberg Equilibrium Test  
**Algorithm Group:** Population Genetics  
**Status:** Complete  
**Last Updated:** 2026-02-01  

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
| HW-S03 | Borderline significance (χ² ≈ 3.84) | Decision boundary testing |

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

### Existing Tests (PopulationGeneticsAnalyzerTests.cs)

| Test | Status | Action |
|------|--------|--------|
| `TestHardyWeinberg_InEquilibrium_PassesTest` | Weak | Strengthen assertions |
| `TestHardyWeinberg_ExcessHeterozygotes_FailsTest` | Incomplete | Add InEquilibrium assertion |
| `TestHardyWeinberg_CalculatesExpectedCorrectly` | Weak | Add more expected value checks |
| `TestHardyWeinberg_ZeroSamples_HandlesGracefully` | OK | Keep |

### Consolidation Plan

1. **Create dedicated test file:** `PopulationGeneticsAnalyzer_HardyWeinberg_Tests.cs`
2. **Remove Hardy-Weinberg tests from `PopulationGeneticsAnalyzerTests.cs`**
3. **Add missing tests:** Published datasets, invariants, significance levels
4. **Strengthen existing tests:** Use Assert.Multiple for comprehensive checks

---

## Open Questions / Decisions

| # | Question | Decision |
|---|----------|----------|
| 1 | What to return for n=0? | Return equilibrium (χ²=0, p=1) - consistent with current implementation |
| 2 | Skip terms in χ² when E=0? | Yes - standard practice to avoid division by zero |
