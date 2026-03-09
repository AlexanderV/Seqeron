# META-BETA-001 Evidence

**Test Unit:** META-BETA-001 (Beta Diversity)  
**Algorithm:** CalculateBetaDiversity  
**Date Collected:** 2026-02-04  
**Last Updated:** 2026-03-09

## Authoritative Sources

### Primary Sources

1. **Whittaker (1960)** — "Vegetation of the Siskiyou Mountains, Oregon and California"
   - Ecological Monographs 30(3): 279–338
   - Original definition of beta diversity concept
   - Source: [Wikipedia Beta Diversity](https://en.wikipedia.org/wiki/Beta_diversity)

2. **Bray & Curtis (1957)** — "An Ordination of the Upland Forest Communities of Southern Wisconsin"
   - Ecological Monographs 27(4): 325–349
   - Original Bray-Curtis dissimilarity formula
   - Source: [Wikipedia Bray–Curtis dissimilarity](https://en.wikipedia.org/wiki/Bray%E2%80%93Curtis_dissimilarity)

3. **Jaccard (1901)** — "Étude comparative de la distribution florale dans une portion des Alpes et des Jura"
   - Bulletin de la Société vaudoise des sciences naturelles 37(142): 547–579
   - Original Jaccard index definition
   - Source: [Wikipedia Jaccard index](https://en.wikipedia.org/wiki/Jaccard_index)

### Key Algorithmic Definitions

#### Bray-Curtis Dissimilarity
- **Formula**: `BC_jk = 1 - 2 * C_jk / (S_j + S_k)`
- **Expanded**: `BC_jk = 1 - 2 * Σ min(N_ij, N_ik) / Σ (N_ij + N_ik)`
- **Where**: N_ij = count of species i at site j; C_jk = sum of lesser counts; S_j, S_k = total specimens per site
- **Range**: [0, 1] where 0 = identical composition, 1 = no shared species
- **Properties**: NOT a true distance metric (violates triangle inequality)
- **Source**: Wikipedia Bray–Curtis dissimilarity, formula + properties

#### Jaccard Distance
- **Formula**: `J_distance = 1 - |A ∩ B| / |A ∪ B|`
- **For binary attributes**: `d_J = (M01 + M10) / (M01 + M10 + M11)` (equivalent)
- **Range**: [0, 1] where 0 = identical sets, 1 = no shared elements
- **Properties**: True distance metric (satisfies triangle inequality)
- **Domain**: Defined for finite non-empty sample sets
- **Source**: Wikipedia Jaccard index, Overview section

## Published Verification Data

### Wikipedia Bray-Curtis Worked Example
**Source**: [Wikipedia Bray–Curtis dissimilarity § Example](https://en.wikipedia.org/wiki/Bray%E2%80%93Curtis_dissimilarity#Example)

| Species | Tank 1 | Tank 2 | Min |
|---|---|---|---|
| Goldfish | 6 | 10 | 6 |
| Guppy | 7 | 0 | 0 |
| Rainbow fish | 4 | 6 | 4 |
| **Total** | **17** | **16** | **10** |

**Bray-Curtis**: C_jk = 6+0+4 = 10; S_j = 17, S_k = 16  
BC = 1 − (2 × 10)/(17 + 16) = 1 − 20/33 = **13/33 ≈ 0.3939**

**Jaccard for same data**: shared = {Goldfish, Rainbow fish} = 2; union = {Goldfish, Guppy, Rainbow fish} = 3  
J_dist = 1 − 2/3 = **1/3 ≈ 0.3333** (Guppy absent in Tank 2 since count = 0)

## Algorithmic Invariants (Source-Verified)

1. **Range Constraints**: Bray-Curtis ∈ [0, 1]; Jaccard Distance ∈ [0, 1]
2. **Symmetry**: Distance(A, B) = Distance(B, A) for both metrics
3. **Identity**: Distance(A, A) = 0 for any valid sample A
4. **Disjoint Maximum**: Completely disjoint samples → distance = 1 for both metrics
5. **Abundance Sensitivity**: Bray-Curtis uses abundance values; Jaccard uses presence/absence only

All derived directly from mathematical definitions in Wikipedia sources.

## Design Decisions

1. **Empty Sample Convention**: Both Bray-Curtis and Jaccard formulas involve division by total counts. When both samples are empty, the denominator is 0, making the formula undefined. Wikipedia Jaccard article states the index is defined for "finite non-empty sample sets." The implementation returns 0 for two identical empty samples (zero distance for identical inputs) and 1.0 when one sample is empty (maximum distance), consistent with scipy.spatial.distance convention.

2. **UniFrac Distance**: Set to 0 as it requires a phylogenetic tree, which is outside the scope of this analysis. The field exists in the `BetaDiversity` record for future extension.

## References

1. Whittaker, R.H. (1960). Vegetation of the Siskiyou Mountains, Oregon and California. Ecological Monographs, 30(3), 279–338.
2. Bray, J.R. & Curtis, J.T. (1957). An ordination of the upland forest communities of southern Wisconsin. Ecological Monographs, 27(4), 325–349.
3. Jaccard, P. (1901). Étude comparative de la distribution florale dans une portion des Alpes et des Jura. Bulletin de la Société vaudoise des sciences naturelles, 37(142), 547–579.
4. Wikipedia contributors. Beta diversity. Wikipedia. Retrieved 2026-03-09.
5. Wikipedia contributors. Bray–Curtis dissimilarity. Wikipedia. Retrieved 2026-03-09.
6. Wikipedia contributors. Jaccard index. Wikipedia. Retrieved 2026-03-09.

## Coverage Classification Changes (2026-03-09)

| Category | Count | Details |
|----------|-------|---------|
| ✅ Covered | 11 | No changes needed |
| ⚠ Weak → Strengthened | 2 | BothSamplesEmpty: added BC=0 and Jaccard=0 assertions (convention per Design Decisions); SymmetryProperty: added exact hand-computed values BC=3/5, Jaccard=2/3, species counts, sample name preservation |
| 🔁 Duplicate → Removed | 1 | DistanceRangeConstraints: 4 scenarios ({A:1}vs{A:1}, {A:1}vs{B:1}, {A:0.3,B:0.7}vs{B:0.4,C:0.6}, {}vs{A:1}) all subsumed by IdenticalSamples, DifferentSingleSpecies, SymmetryProperty, OneSampleEmpty respectively; range [0,1] implicitly verified by 12 exact-value tests |
| ❌ Missing | 0 | All spec items covered |