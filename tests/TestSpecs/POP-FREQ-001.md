# Test Specification: POP-FREQ-001

**Test Unit ID:** POP-FREQ-001
**Title:** Allele Frequencies
**Algorithm Group:** Population Genetics
**Status:** Complete
**Last Updated:** 2026-03-08

---

## Scope

### Canonical Methods (Deep Testing Required)

| Method | Class | Complexity |
|--------|-------|------------|
| `CalculateAlleleFrequencies(int, int, int)` | PopulationGeneticsAnalyzer | O(1) |
| `CalculateMAF(IEnumerable<int>)` | PopulationGeneticsAnalyzer | O(n) |
| `FilterByMAF(IEnumerable<Variant>, double, double)` | PopulationGeneticsAnalyzer | O(n) |

---

## Test Cases

### MUST Tests (Required - Evidence-Based)

#### CalculateAlleleFrequencies

| ID | Test Case | Evidence Source | Expected |
|----|-----------|-----------------|----------|
| AF-M01 | Frequencies sum to 1.0 | Wikipedia: "p + q = 1" | major + minor ≈ 1.0 |
| AF-M02 | Wikipedia flower example (49-42-9) | Wikipedia Genotype Frequency | p=0.70, q=0.30 |
| AF-M03 | All homozygous major → (1.0, 0.0) | Derived from formula | (1.0, 0.0) |
| AF-M04 | All homozygous minor → (0.0, 1.0) | Derived from formula | (0.0, 1.0) |
| AF-M05 | All heterozygous → (0.5, 0.5) | Derived from formula | (0.5, 0.5) |
| AF-M06 | Equal genotypes 25-50-25 → (0.5, 0.5) | HWE expectation | (0.5, 0.5) |
| AF-M07 | Zero samples → (0, 0) | Edge case handling | (0, 0) |
| AF-M08 | Frequencies are non-negative | Range constraint | major ≥ 0, minor ≥ 0 |
| AF-M09 | Frequencies do not exceed 1.0 | Range constraint | major ≤ 1, minor ≤ 1 |
| AF-M10 | Wikipedia diploid example (6-3-1) | Wikipedia Allele Frequency § Diploids → Example | p=0.75, q=0.25 |
| AF-M11 | NDSU blood type example (1787-3039-1303) | NDSU McClean (1998) via Wayback | p≈0.5394, q≈0.4606 |
| AF-M12 | NDSU molecular example (30-50-20) | NDSU McClean (1998) via Wayback | p=0.55, q=0.45 |
| AF-M13 | Negative hom_maj count → throws | Input validation | ArgumentOutOfRangeException |
| AF-M14 | Negative heterozygous count → throws | Input validation | ArgumentOutOfRangeException |
| AF-M15 | Negative hom_min count → throws | Input validation | ArgumentOutOfRangeException |

#### CalculateMAF

| ID | Test Case | Evidence Source | Expected |
|----|-----------|-----------------|----------|
| MAF-M01 | Alt frequency < 0.5 returns alt freq | MAF definition | MAF = alt_freq |
| MAF-M02 | Alt frequency > 0.5 returns 1-alt_freq | MAF definition | MAF = 1 - alt_freq |
| MAF-M03 | MAF is always ≤ 0.5 | MAF invariant | MAF ≤ 0.5 |
| MAF-M04 | MAF is always ≥ 0 | MAF invariant | MAF ≥ 0 |
| MAF-M05 | Monomorphic ref (all 0) → MAF = 0 | Fixed allele | MAF = 0 |
| MAF-M06 | Monomorphic alt (all 2) → MAF = 0 | Fixed allele | MAF = 0 |
| MAF-M07 | Empty genotypes → MAF = 0 | Edge case | MAF = 0 |
| MAF-M08 | Perfect 50/50 split → MAF = 0.5 | Boundary | MAF = 0.5 |

#### FilterByMAF

| ID | Test Case | Evidence Source | Expected |
|----|-----------|-----------------|----------|
| FLT-M01 | Filters variants below minMAF | Filter logic | Excluded |
| FLT-M02 | Keeps variants at/above minMAF | Filter logic | Included |
| FLT-M03 | Filters variants above maxMAF | Filter logic | Excluded |
| FLT-M04 | Keeps variants at/below maxMAF | Filter logic | Included |
| FLT-M05 | Empty input → empty output | Edge case | Empty enumerable |
| FLT-M06 | All filtered → empty output | Edge case | Empty enumerable |
| FLT-M07 | None filtered → all pass | Edge case | All variants |

### SHOULD Tests (Recommended)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| AF-S01 | Large population (10000 samples) | No overflow |
| AF-S02 | Single sample homozygous | Minimum valid input |
| MAF-S01 | Various genotype patterns | Algorithm verification |
| FLT-S01 | Boundary MAF values (at exact threshold) | Precision handling |
| FLT-S02 | Filter preserves input order | Ordering guarantee |

### COULD Tests (Optional)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| AF-C01 | Property test: sum invariant | Randomized testing |
| AF-C02 | Property test: range invariant [0,1] | Randomized testing |
| MAF-C01 | Property test: MAF range invariant | Randomized testing |

---

## Invariants to Verify

### Allele Frequency Invariants

1. ∀ (hom_maj, het, hom_min): major + minor = 1.0 (when total > 0)
2. ∀ result: 0 ≤ major ≤ 1 ∧ 0 ≤ minor ≤ 1
3. major = (2×hom_maj + het) / (2×total)
4. minor = (2×hom_min + het) / (2×total)

### MAF Invariants

1. ∀ genotypes: 0 ≤ MAF ≤ 0.5
2. ∀ genotypes: MAF = min(alt_freq, 1 - alt_freq)
3. MAF of monomorphic locus = 0

### Filter Invariants

1. ∀ filtered: minMAF ≤ MAF(v) ≤ maxMAF
2. FilterByMAF preserves order
3. FilterByMAF is lazy (IEnumerable)

---

## Deviations and Assumptions

None. All formulas, behaviors, and edge cases are sourced from external references:

- Wikipedia: Allele frequency (accessed 2026-02-01, verified 2026-03-08)
- Wikipedia: Minor allele frequency (accessed 2026-02-01, verified 2026-03-08)
- Wikipedia: Genotype frequency (accessed 2026-02-01, verified 2026-03-08)
- Gillespie (2004) Population Genetics: A Concise Guide, ISBN 978-0-8018-8008-7
- NDSU Population Genetics (McClean, 1998) via Wayback Machine (archived 2024-05-12)

Negative genotype counts are rejected with `ArgumentOutOfRangeException` — genotype counts are non-negative by definition.

**Canonical file:** `PopulationGeneticsAnalyzer_AlleleFrequency_Tests.cs` — 49 test cases (30 MUST, 11 SHOULD, 8 invariant/parameterized).
**Property file:** `Properties/PopulationGeneticsProperties.cs` — AF-C01 (sum, 1e-10), AF-C02 (range [0,1]), MAF-C01 (MAF range [0, 0.5]).

### Coverage Classification (2026-03-08)

| ID | Status | Notes |
|----|--------|-------|
| AF-M01–M15 | ✅ Covered | All MUST tests present with exact assertions (Within(1e-10) or exact fractions) |
| AF-S01 | ✅ Covered | Large population — exact frequency + sum + NaN checks |
| AF-S02 | ✅ Covered | Single sample — 3 parameterized cases with exact values |
| MAF-M01–M08 | ✅ Covered | Exact fractions or Within(1e-10) |
| MAF-S01 | ✅ Covered | 3 patterns with Within(1e-10) |
| FLT-M01–M07 | ✅ Covered | All filter edge cases |
| FLT-S01 | ✅ Covered | High AF MAF calculation |
| FLT-S02 | ✅ Covered | Order preservation |
| AF-C01 | ✅ Covered | Property: sum invariant — Within(1e-10) (PopulationGeneticsProperties) |
| AF-C02 | ✅ Covered | Property: range invariant [0,1] (PopulationGeneticsProperties) |
| MAF-C01 | ✅ Covered | Property: MAF range invariant (PopulationGeneticsProperties) |
