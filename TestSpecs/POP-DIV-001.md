# Test Specification: POP-DIV-001

**Test Unit ID:** POP-DIV-001  
**Title:** Diversity Statistics  
**Algorithm Group:** Population Genetics  
**Status:** Complete  
**Last Updated:** 2026-02-01  

---

## Scope

### Canonical Methods (Deep Testing Required)

| Method | Class | Complexity |
|--------|-------|------------|
| `CalculateNucleotideDiversity(seqs)` | PopulationGeneticsAnalyzer | O(n² × m) |
| `CalculateWattersonTheta(S, n, L)` | PopulationGeneticsAnalyzer | O(n) |
| `CalculateTajimasD(π, θ, S, n)` | PopulationGeneticsAnalyzer | O(n) |
| `CalculateDiversityStatistics(seqs)` | PopulationGeneticsAnalyzer | O(n² × m) |

---

## Test Cases

### MUST Tests (Required - Evidence-Based)

#### CalculateNucleotideDiversity

| ID | Test Case | Evidence Source | Expected |
|----|-----------|-----------------|----------|
| ND-M01 | Identical sequences → π = 0 | Definition (Nei & Li 1979) | 0 |
| ND-M02 | All different (2 seqs, all positions) → π = 1.0 | Definition | 1.0 |
| ND-M03 | Single sequence → π = 0 | Edge case (n < 2) | 0 |
| ND-M04 | Wikipedia example (n=5, 20 pairwise diffs, 10 comps, L=20) | Wikipedia Tajima's D | k̂ = 2.0, π = 0.1 |
| ND-M05 | Empty input → π = 0 | Edge case | 0 |
| ND-M06 | π is always non-negative | Range invariant | π ≥ 0 |

#### CalculateWattersonTheta

| ID | Test Case | Evidence Source | Expected |
|----|-----------|-----------------|----------|
| WT-M01 | S=10, n=10, L=1000 → θ ≈ 0.00353 | Wikipedia Watterson | 0.00353 (±0.0005) |
| WT-M02 | S=0 → θ = 0 | Definition | 0 |
| WT-M03 | n < 2 → θ = 0 | Edge case | 0 |
| WT-M04 | L ≤ 0 → θ = 0 | Edge case | 0 |
| WT-M05 | n=2 (minimum valid) → θ = S / (1 × L) | a₁ = 1 for n=2 | S/L |
| WT-M06 | θ is always non-negative | Range invariant | θ ≥ 0 |

#### CalculateTajimasD

| ID | Test Case | Evidence Source | Expected |
|----|-----------|-----------------|----------|
| TD-M01 | π = θ → D ≈ 0 (near zero) | Wikipedia: neutral evolution | |D| < 1 |
| TD-M02 | π << θ → D < 0 (negative) | Wikipedia: positive selection | D < 0 |
| TD-M03 | π >> θ → D > 0 (positive) | Wikipedia: balancing selection | D > 0 |
| TD-M04 | S = 0 → D = 0 | Edge case | 0 |
| TD-M05 | n < 3 → D = 0 | Edge case (minimum n=3) | 0 |
| TD-M06 | Variance ≤ 0 → D = 0 | Implementation guard | 0 |

#### CalculateDiversityStatistics

| ID | Test Case | Evidence Source | Expected |
|----|-----------|-----------------|----------|
| DS-M01 | Returns all metrics (SampleSize, S, π, θ, D) | Integration | All fields populated |
| DS-M02 | Single sequence → all zeros | Edge case | π=0, θ=0, D=0, S=0 |
| DS-M03 | Empty input → all zeros | Edge case | All zeros |
| DS-M04 | SampleSize matches input count | Consistency | n = seqList.Count |
| DS-M05 | SegregratingSites correctly counted | Definition | Matches manual count |

### SHOULD Tests (Recommended)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| ND-S01 | 3+ sequences with partial polymorphism | Real-world scenario |
| WT-S01 | Various sample sizes (n=3,5,10,20) | Harmonic number accuracy |
| TD-S01 | Moderate D values (selection signal) | Algorithm verification |
| DS-S01 | Heterozygosity values are in [0,1] | Range validation |
| DS-S02 | H_expected ≥ H_observed for polymorphic data | Theoretical relationship |

### COULD Tests (Optional)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| ND-C01 | Property test: π invariant under sequence reordering | Symmetry |
| WT-C01 | Large sample performance (n=100) | Performance baseline |
| TD-C01 | Wikipedia example full calculation | Reference validation |

---

## Invariants to Verify

### Nucleotide Diversity Invariants

1. ∀ sequences: π ≥ 0
2. ∀ identical sequences: π = 0
3. ∀ n < 2: π = 0 (undefined, return 0)
4. π = (total pairwise differences) / (comparisons × length)

### Watterson's Theta Invariants

1. ∀ parameters: θ ≥ 0
2. ∀ S = 0: θ = 0
3. ∀ n < 2 or L ≤ 0: θ = 0
4. θ = S / (a_n × L) where a_n = Σ(1/i) for i=1 to n-1

### Tajima's D Invariants

1. ∀ S = 0: D = 0
2. ∀ n < 3: D = 0 (undefined)
3. When π ≈ θ: D ≈ 0 (within statistical variance)
4. When π < θ: D < 0 (typically)
5. When π > θ: D > 0 (typically)

### DiversityStatistics Invariants

1. SampleSize = number of input sequences
2. SegregratingSites = count of polymorphic positions
3. All numeric values ≥ 0
4. 0 ≤ HeterozygosityObserved ≤ 1
5. 0 ≤ HeterozygosityExpected ≤ 1

---

## Audit Summary

### Existing Tests (PopulationGeneticsAnalyzerTests.cs)

| Test | Status | Notes |
|------|--------|-------|
| CalculateNucleotideDiversity_IdenticalSequences_ReturnsZero | Keep | ND-M01 |
| CalculateNucleotideDiversity_AllDifferent_ReturnsPositive | Keep | ND-M02 |
| CalculateNucleotideDiversity_SingleSequence_ReturnsZero | Duplicate | Remove, redundant |
| CalculateWattersonTheta_KnownValues_CalculatesCorrectly | Keep | WT-M01 |
| CalculateWattersonTheta_SmallSample_HandlesCorrectly | Keep | WT-M05 variant |
| CalculateTajimasD_NeutralEvolution_NearZero | Keep | TD-M01 |
| CalculateTajimasD_PositiveSelection_Negative | Keep | TD-M02 |
| CalculateTajimasD_NoSegregratingSites_ReturnsZero | Keep | TD-M04 |
| CalculateDiversityStatistics_ReturnsAllMetrics | Keep | DS-M01 |
| CalculateDiversityStatistics_SingleSequence_ReturnsZeroDiversity | Keep | DS-M02 |

### Consolidation Plan

1. **Create new file:** `PopulationGeneticsAnalyzer_Diversity_Tests.cs`
2. **Move and enhance:** All diversity-related tests from `PopulationGeneticsAnalyzerTests.cs`
3. **Remove duplicates:** Merge redundant single-sequence tests
4. **Add missing:** Empty input, edge cases, invariant tests
5. **Organize:** Group by method with clear regions

---

## Open Questions / Decisions

1. **Heterozygosity adaptation:** Implementation uses polymorphic site fraction for observed heterozygosity (sequence data adaptation). Documented as ASSUMPTION - acceptable for haploid sequence data use case.

2. **Typo in record:** `SegregratingSites` has a typo (should be `SegregatingSites`). Not fixing to maintain API compatibility.

---

## Test File

`Seqeron.Genomics.Tests/PopulationGeneticsAnalyzer_Diversity_Tests.cs`
