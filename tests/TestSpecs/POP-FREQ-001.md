# Test Specification: POP-FREQ-001

**Test Unit ID:** POP-FREQ-001
**Title:** Allele Frequencies
**Algorithm Group:** Population Genetics
**Status:** Complete
**Last Updated:** 2026-02-01

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

### COULD Tests (Optional)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| AF-C01 | Property test: sum invariant | Randomized testing |
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

## Audit of Existing Tests

### PopulationGeneticsAnalyzerTests.cs (Lines 12-86)

| Existing Test | Classification | Action |
|---------------|----------------|--------|
| CalculateAlleleFrequencies_EqualGenotypes_Returns50Percent | Covered (AF-M06) | Keep, enhance |
| CalculateAlleleFrequencies_AllHomozygousMajor_Returns100Percent | Covered (AF-M03) | Keep |
| CalculateAlleleFrequencies_ZeroSamples_ReturnsZero | Covered (AF-M07) | Keep |
| CalculateMAF_FromGenotypes_CalculatesCorrectly | Covered (MAF-M01) | Keep |
| CalculateMAF_HighAltFreq_ReturnsMinorAllele | Covered (MAF-M06) | Keep |
| FilterByMAF_FiltersCorrectly | Weak (FLT-M01/M02 combined) | Enhance |

### Missing Tests (All Closed)

- ~~AF-M02 (Wikipedia example)~~ ✅ Covered
- ~~AF-M04 (All homozygous minor)~~ ✅ Covered
- ~~AF-M05 (All heterozygous)~~ ✅ Covered
- ~~AF-M01 (Sum to 1 invariant)~~ ✅ Covered
- ~~AF-M08, AF-M09 (Range constraints)~~ ✅ Covered
- ~~MAF-M02 through MAF-M05, MAF-M07, MAF-M08~~ ✅ Covered
- ~~FLT-M03 through FLT-M07~~ ✅ Covered

### Consolidation Plan

1. Create dedicated test file: `PopulationGeneticsAnalyzer_AlleleFrequency_Tests.cs`
2. Move and enhance existing allele frequency tests
3. Add all missing MUST tests
4. Remove allele frequency tests from original file (to avoid duplication)

---

## Open Questions / Decisions

None - algorithm behavior is well-defined.
