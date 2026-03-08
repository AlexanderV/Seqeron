# Test Specification: POP-DIV-001

**Test Unit ID:** POP-DIV-001
**Title:** Diversity Statistics
**Algorithm Group:** Population Genetics
**Status:** Complete
**Last Updated:** 2026-03-08

---

## Scope

### Canonical Methods (Deep Testing Required)

| Method | Class | Complexity |
|--------|-------|------------|
| `CalculateNucleotideDiversity(seqs)` | PopulationGeneticsAnalyzer | O(n² × m) |
| `CalculateWattersonTheta(S, n, L)` | PopulationGeneticsAnalyzer | O(n) |
| `CalculateTajimasD(k̂, S, n)` | PopulationGeneticsAnalyzer | O(n) |
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
| ND-M04 | Wikipedia Tajima's D example (n=5, L=20, total pairwise diffs=20, comparisons=10) | Wikipedia Tajima's D | k̂ = 2.0, π = 0.1 |
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
| TD-M01 | k̂ = S/a₁ (neutral evolution) → D = 0 | Wikipedia: neutral evolution | D = 0 |
| TD-M02 | k̂ << S/a₁ → D < 0 (negative) | Wikipedia: positive selection | D < 0 |
| TD-M03 | k̂ >> S/a₁ → D > 0 (positive) | Wikipedia: balancing selection | D > 0 |
| TD-M04 | S = 0 → D = 0 | Edge case | 0 |
| TD-M05 | n < 3 → D = 0 | Edge case (minimum n=3) | 0 |
| TD-M06 | n = 3 → produces valid (non-NaN, non-Infinity) D | Minimum valid sample | Valid float |

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
| DS-S01 | Heterozygosity values are in [0,1] | Range validation |
| DS-S02 | Identical sequences → zero segregating sites | Zero-diversity invariant |

### COULD Tests (Reference Validation)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| TD-C01 | Wikipedia Tajima's D full hand-calculation (n=5, S=4, k̂=2.0 → D ≈ 0.273) | Exact numerical verification against Wikipedia |
| TD-C02 | End-to-end: Wikipedia sequences → CalculateDiversityStatistics → D ≈ 0.273, π = 0.1 | Full pipeline verification |

---

## Invariants to Verify

### Nucleotide Diversity Invariants

1. ∀ sequences: π ≥ 0
2. ∀ identical sequences: π = 0
3. ∀ n < 2: π = 0 (undefined, return 0)
4. π = Σ d_ij / (C(n,2) × L)

### Watterson's Theta Invariants

1. ∀ parameters: θ ≥ 0
2. ∀ S = 0: θ = 0
3. ∀ n < 2 or L ≤ 0: θ = 0
4. θ = S / (a_n × L) where a_n = Σ(1/i) for i=1 to n-1

### Tajima's D Invariants

1. ∀ S = 0: D = 0
2. ∀ n < 3: D = 0 (undefined)
3. When k̂ = S/a₁: D = 0
4. When k̂ < S/a₁: D < 0
5. When k̂ > S/a₁: D > 0

### DiversityStatistics Invariants

1. SampleSize = number of input sequences
2. SegregratingSites = count of polymorphic positions
3. All numeric values ≥ 0
4. 0 ≤ HeterozygosityObserved ≤ 1
5. 0 ≤ HeterozygosityExpected ≤ 1
6. HeterozygosityObserved = n/(n−1) × HeterozygosityExpected (Nei bias correction)

---

## Notes

1. **Tajima's D API:** `CalculateTajimasD(averagePairwiseDifferences, segregatingSites, sampleSize)` takes unnormalized k̂ (average number of pairwise differences per pair of sequences, NOT per-site π). The Watterson estimate S/a₁ is computed internally. This matches the Wikipedia formula: D = (k̂ − S/a₁) / √(e₁S + e₂S(S−1)).

2. **Typo in record:** `SegregratingSites` has a typo (should be `SegregatingSites`). Not fixing to maintain API compatibility.

3. **Heterozygosity for haploid sequence data:** `HeterozygosityObserved` = Nei's (1978) unbiased gene diversity: n/(n−1) × (1−Σp²) averaged per site. `HeterozygosityExpected` = basic gene diversity: (1−Σp²) averaged per site (Wikipedia: Zygosity). Both are standard formulas with well-defined semantics for haploid data. For haploid sequences, Nei's unbiased gene diversity is mathematically equivalent to nucleotide diversity π.

---

## Coverage Classification

| Test Method | Spec ID | Status | Notes |
|-------------|---------|--------|-------|
| `CalculateNucleotideDiversity_IdenticalSequences_ReturnsZero` | ND-M01 | ✅ Covered | Exact value: 0 |
| `CalculateNucleotideDiversity_AllDifferent_ReturnsOne` | ND-M02 | ✅ Covered | Exact value: 1.0 |
| `CalculateNucleotideDiversity_SingleSequence_ReturnsZero` | ND-M03 | ✅ Covered | |
| `CalculateNucleotideDiversity_WikipediaExample_CalculatesCorrectly` | ND-M04 | ✅ Covered | π = 0.1 (±0.0001) |
| `CalculateNucleotideDiversity_EmptyInput_ReturnsZero` | ND-M05 | ✅ Covered | |
| `CalculateNucleotideDiversity_VariousInputs_AlwaysNonNegative` | ND-M06 | ✅ Covered | Invariant check |
| `CalculateNucleotideDiversity_PartialPolymorphism_CalculatesCorrectly` | ND-S01 | ✅ Covered | Exact: 2/12 |
| `CalculateWattersonTheta_WikipediaExample_CalculatesCorrectly` | WT-M01 | ✅ Covered | θ ≈ 0.00353 |
| `CalculateWattersonTheta_ZeroSegregatingSites_ReturnsZero` | WT-M02 | ✅ Covered | |
| `CalculateWattersonTheta_SampleSizeLessThanTwo_ReturnsZero` | WT-M03 | ✅ Covered | |
| `CalculateWattersonTheta_InvalidSequenceLength_ReturnsZero` | WT-M04 | ✅ Covered | |
| `CalculateWattersonTheta_MinimumSampleSize_CalculatesCorrectly` | WT-M05 | ✅ Covered | Exact: S/L |
| `CalculateWattersonTheta_VariousInputs_AlwaysNonNegative` | WT-M06 | ✅ Covered | Invariant check |
| `CalculateWattersonTheta_VariousSampleSizes_HarmonicNumberCorrect` | WT-S01 | ✅ Covered | Exact values for n=3,5,10 |
| `CalculateTajimasD_NeutralEvolution_NearZero` | TD-M01 | ✅ Covered | D = 0 when k̂ = S/a₁ |
| `CalculateTajimasD_PositiveSelection_Negative` | TD-M02 | ✅ Covered | Sign check |
| `CalculateTajimasD_BalancingSelection_Positive` | TD-M03 | ✅ Covered | Sign check |
| `CalculateTajimasD_NoSegregatingSites_ReturnsZero` | TD-M04 | ✅ Covered | |
| `CalculateTajimasD_SampleSizeLessThanThree_ReturnsZero` | TD-M05 | ✅ Covered | |
| `CalculateTajimasD_MinimumValidSampleSize_CalculatesValue` | TD-M06 | ✅ Covered | Non-NaN, non-Infinity |
| `CalculateTajimasD_WikipediaExample_ExactValue` | TD-C01 | ✅ Covered | D ≈ 0.273 (±0.005) |
| `CalculateDiversityStatistics_WikipediaExample_CorrectTajimasD` | TD-C02 | ✅ Covered | Full pipeline |
| `CalculateDiversityStatistics_ReturnsAllMetrics` | DS-M01 | ✅ Covered | All 7 fields exact |
| `CalculateDiversityStatistics_SingleSequence_ReturnsZeros` | DS-M02 | ✅ Covered | |
| `CalculateDiversityStatistics_EmptyInput_ReturnsZeros` | DS-M03 | ✅ Covered | |
| `CalculateDiversityStatistics_SampleSizeMatchesInputCount` | DS-M04 | ✅ Covered | |
| `CalculateDiversityStatistics_SegregatingSitesCountedCorrectly` | DS-M05 | ✅ Covered | Exact S = 2 |
| `CalculateDiversityStatistics_HeterozygosityInValidRange` | DS-S01 | ✅ Covered | Range [0,1] |
| `CalculateDiversityStatistics_IdenticalSequences_ZeroSegregatingSites` | DS-S02 | ✅ Covered | |

**Total:** 29 tests, 29 ✅ Covered, 0 ❌ Missing, 0 ⚠ Weak, 0 🔁 Duplicate

### Changes Made

| Action | Test | Details |
|--------|------|---------|
| ⚠→✅ Strengthened | DS-M01 | Replaced permissive `GreaterThanOrEqualTo(0)` with exact hand-computed values for all 7 fields; added TajimasD assertion |
| ⚠→✅ Strengthened | WT-S01 | Added exact θ values for n=5 (0.048) and n=10 (0.03535) alongside ordering assertions |
| 🔁 Removed | `DiversePopulation_AllMetricsReasonable` | Redundant with strengthened DS-M01 + TD-C02; had only weak `GreaterThan(0)` assertions |

---

## Test File

`Seqeron.Genomics.Tests/PopulationGeneticsAnalyzer_Diversity_Tests.cs`
