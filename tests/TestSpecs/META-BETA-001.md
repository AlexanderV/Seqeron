# META-BETA-001: Beta Diversity Test Specification

**Test Unit ID:** META-BETA-001  
**Algorithm:** Beta Diversity (CalculateBetaDiversity)  
**Class:** MetagenomicsAnalyzer  
**Method:** `CalculateBetaDiversity(string, IReadOnlyDictionary<string, double>, string, IReadOnlyDictionary<string, double>)`  
**Created:** 2026-02-04  

## Test Categories

### MUST Tests (Evidence-Based)

#### M1: Mathematical Correctness
**Source:** Bray & Curtis (1957), Jaccard (1901)

1. **Identical Samples Return Zero Distance**
   - Bray-Curtis = 0 when samples identical
   - Jaccard Distance = 0 when samples identical
   - Both sample names preserved in result

2. **Completely Disjoint Samples Return Maximum Distance**
   - Bray-Curtis = 1 when no shared species
   - Jaccard Distance = 1 when no shared species
   - Correct species count statistics

3. **Symmetry Property**
   - Distance(A, B) = Distance(B, A) for both metrics
   - Sample name order preserved correctly

#### M2: Edge Case Handling
**Source:** Mathematical analysis of formulas

4. **Empty Sample Handling**
   - Both samples empty → defined behavior (implementation choice)
   - One sample empty → maximum distance = 1
   - No runtime exceptions

5. **Single Species Samples**
   - Same single species → distance = 0
   - Different single species → distance = 1

#### M3: Range Validation  
**Source:** Mathematical definitions

6. **Distance Range Constraints**
   - Bray-Curtis ∈ [0, 1] for all valid inputs
   - Jaccard Distance ∈ [0, 1] for all valid inputs

7. **Species Count Accuracy**
   - SharedSpecies = intersection count
   - UniqueToSample1 + UniqueToSample2 + SharedSpecies = total unique species
   - Counts match actual species overlap

### SHOULD Tests (Quality Assurance)

#### S1: Realistic Ecological Scenarios

8. **Partial Species Overlap**
   - Mixed composition with known overlap percentage
   - Verify intermediate distance values are reasonable

9. **Abundance vs. Presence Effects**
   - Same species, different abundances
   - Verify Bray-Curtis responds to abundance differences
   - Verify Jaccard ignores abundance (presence/absence only)

#### S2: Boundary Conditions

10. **Zero Abundance Handling**
    - Species with 0 abundance treated as absent
    - Mixed zero/non-zero abundances processed correctly

11. **Very Small/Large Abundances**
    - Numerical stability with extreme values
    - Proportional calculations maintain precision

### COULD Tests (Extended Coverage)

#### C1: Performance

12. **Large Sample Processing**
    - Many species (1000+) processed efficiently
    - Memory usage reasonable

13. **Complex Sample Names**
    - Unicode characters, special symbols
    - Long sample names preserved correctly

## Current Test Audit

### Existing Coverage Analysis
**Source:** MetagenomicsAnalyzerTests.cs

**Covered:**
- ✅ M1.1: Identical samples → zero distance
- ✅ M2.1: Completely different samples → max distance  
- ✅ S1.1: Partial overlap with intermediate distance
- ✅ Basic sample name preservation

**Missing/Weak:**
- ❌ M1.3: Symmetry property not explicitly tested
- ❌ M2.4: Empty sample handling not tested
- ❌ M2.5: Single species scenarios not tested
- ❌ M3.6: Range constraints not validated with Assert.That(..., Is.InRange(0, 1))
- ❌ M3.7: Species count accuracy not validated 
- ❌ S1.2: Abundance vs presence distinction not tested
- ❌ S2.1: Zero abundance handling not tested

### Test Integration Plan

**Canonical Test File:** `MetagenomicsAnalyzer_BetaDiversity_Tests.cs` (new dedicated file)  
**Current Tests:** Migrate from MetagenomicsAnalyzerTests.cs #region Beta Diversity Tests  
**Wrapper Tests:** Keep minimal smoke test in original location

## Implementation Notes

- Use mathematical examples with known outcomes for validation
- Focus on algorithmic correctness over biological realism
- Ensure deterministic test results
- Use Assert.Multiple for related invariants
- Test both metrics independently and in combination

## Open Questions

1. **Empty Sample Behavior**: Current implementation behavior undefined - should be specified
2. **UniFrac Distance**: Always returns 0 - should this be tested or documented as not implemented?
3. **Precision Requirements**: No documented precision requirements for floating-point comparisons