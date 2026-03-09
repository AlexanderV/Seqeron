# META-BETA-001: Beta Diversity Test Specification

**Test Unit ID:** META-BETA-001  
**Algorithm:** Beta Diversity (CalculateBetaDiversity)  
**Class:** MetagenomicsAnalyzer  
**Method:** `CalculateBetaDiversity(string, IReadOnlyDictionary<string, double>, string, IReadOnlyDictionary<string, double>)`  
**Status:** ☑ Complete  
**Created:** 2026-02-04  
**Last Updated:** 2026-03-09

## Test Categories

### MUST Tests (Evidence-Based)

#### M1: Published Example Verification
**Source:** Wikipedia Bray–Curtis dissimilarity § Example

1. **Wikipedia Aquarium Example**
   - Tank 1: Goldfish=6, Guppy=7, Rainbow fish=4; Tank 2: Goldfish=10, Guppy=0, Rainbow fish=6
   - Bray-Curtis = 13/33 ≈ 0.3939 (Wikipedia published result)
   - Jaccard Distance = 1/3 ≈ 0.3333 (computed from definition)
   - SharedSpecies=2, UniqueToSample1=1, UniqueToSample2=0

#### M2: Mathematical Correctness
**Source:** Bray & Curtis (1957), Jaccard (1901)

2. **Identical Samples Return Zero Distance**
   - Bray-Curtis = 0 when samples identical
   - Jaccard Distance = 0 when samples identical
   - Both sample names preserved in result

3. **Completely Disjoint Samples Return Maximum Distance**
   - Bray-Curtis = 1 when no shared species
   - Jaccard Distance = 1 when no shared species
   - Correct species count statistics

4. **Symmetry Property**
   - Distance(A, B) = Distance(B, A) for both metrics
   - Sample name order preserved correctly

#### M3: Edge Case Handling
**Source:** Mathematical analysis of formulas; Wikipedia Jaccard index (domain: non-empty sets)

5. **Empty Sample Handling**
   - Both samples empty → returns 0 (convention: identical inputs, see Design Decisions in Evidence)
   - One sample empty → maximum distance = 1
   - No runtime exceptions

6. **Single Species Samples**
   - Same single species, same abundance → distance = 0
   - Same single species, different abundance → BC = 5/11, Jaccard = 0
   - Different single species → distance = 1

7. **Zero Abundance Handling**
   - Species with 0 abundance treated as absent (per Jaccard presence/absence definition)
   - BC and Jaccard computed only over species with nonzero abundance
   - Exact values verified: BC=0.7, Jaccard=2/3 for test scenario

#### M4: Range and Count Validation
**Source:** Mathematical definitions

8. **Distance Range Constraints**
   - Bray-Curtis ∈ [0, 1] for all valid inputs
   - Jaccard Distance ∈ [0, 1] for all valid inputs

9. **Species Count Accuracy**
   - SharedSpecies + UniqueToSample1 + UniqueToSample2 = total unique species
   - Exact BC and Jaccard values verified for multi-species scenarios

### SHOULD Tests (Quality Assurance)

#### S1: Abundance vs. Presence Distinction

10. **Bray-Curtis Responds to Abundance, Jaccard Ignores It**
    - Same species present in both samples, different abundance distributions
    - Balanced vs sample2: BC=0 (identical); Skewed vs sample2: BC=0.4
    - Jaccard identical in both comparisons (same species sets)

#### S2: Sample Name Preservation

11. **Complex Sample Names Preserved**
    - Unicode characters, special symbols in names
    - Names preserved exactly in result

## Current Test Coverage

**Canonical Test File:** `MetagenomicsAnalyzer_BetaDiversity_Tests.cs`

| Test | Spec | Status |
|------|------|--------|
| WikipediaAquariumExample_MatchesPublishedResult | M1.1 | ✅ |
| IdenticalSamples_ReturnsZeroDistance | M2.2 | ✅ |
| CompletelyDisjointSamples_ReturnsMaximumDistance | M2.3 | ✅ |
| SymmetryProperty_ReturnsIdenticalResults | M2.4 | ✅ |
| BothSamplesEmpty_HandlesGracefully | M3.5 | ✅ |
| OneSampleEmpty_ReturnsMaximumDistance | M3.5 | ✅ |
| SameSingleSpeciesSameAbundance_ReturnsZeroDistance | M3.6 | ✅ |
| SameSingleSpeciesDifferentAbundance_BrayCurtisNonZero | M3.6 | ✅ |
| DifferentSingleSpecies_ReturnsMaximumDistance | M3.6 | ✅ |
| ZeroAbundanceHandling_TreatsAsAbsent | M3.7 | ✅ |
| SpeciesCountAccuracy_MatchesExpected | M4.9 | ✅ |
| BrayCurtisRespondsToAbundance_JaccardIgnoresIt | S1.10 | ✅ |
| ComplexSampleNames_PreservedCorrectly | S2.11 | ✅ |

## Coverage Classification (2026-03-09)

| Category | Count | Action |
|----------|-------|--------|
| ✅ Covered | 11 | No changes |
| ⚠ Weak → Strengthened | 2 | BothSamplesEmpty: added BC=0, Jaccard=0 assertions; SymmetryProperty: added exact BC=3/5, Jaccard=2/3 + species counts + name checks |
| 🔁 Duplicate → Removed | 1 | DistanceRangeConstraints (all 4 scenarios subsumed by other exact-value tests; range [0,1] implicitly verified) |
| ❌ Missing | 0 | — |

**M4.8 (Range Constraints):** Verified implicitly by all 12 exact-value tests that assert specific values within [0, 1]. Dedicated range-only test removed as redundant.

## Deviations and Assumptions

**None.** All test expected values are derived directly from:
- The published Wikipedia Bray-Curtis worked example (Tank 1/Tank 2)
- Hand computation from the source formulas (Bray-Curtis: `1 - 2*C/(S_j+S_k)`; Jaccard: `1 - |A∩B|/|A∪B|`)

Empty sample handling is a documented design decision (see Evidence), not an assumption.