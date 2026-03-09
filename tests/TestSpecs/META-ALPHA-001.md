# META-ALPHA-001: Alpha Diversity Test Specification

## Test Unit Information

| Field | Value |
|-------|-------|
| **ID** | META-ALPHA-001 |
| **Area** | Metagenomics |
| **Canonical Methods** | `MetagenomicsAnalyzer.CalculateAlphaDiversity` |
| **Complexity** | O(n) where n = number of taxa |
| **Invariants** | H ≥ 0; λ ∈ [0,1]; J ∈ [0,1]; S ≥ 0 |
| **Status** | ☑ Complete |

## Methods Under Test

| Method | Class | Type | Test Depth |
|--------|-------|------|------------|
| `CalculateAlphaDiversity(abundances)` | MetagenomicsAnalyzer | Canonical | Deep |

## Evidence Sources

1. **Wikipedia — Diversity Index:** Shannon, Simpson, Inverse Simpson formulas
2. **Wikipedia — Alpha Diversity:** Concept definition
3. **Wikipedia — Species Richness:** Observed species count
4. **Wikipedia — Species Evenness:** Pielou evenness formula
5. **Shannon (1948):** H = −Σ pᵢ ln(pᵢ)
6. **Simpson (1949):** λ = Σ pᵢ²
7. **Hill (1973):** Unified diversity framework
8. **Chao (1984):** Chao1 estimator

## Test Categories

### MUST Tests (Evidence-Backed)

| ID | Test | Evidence | Justification |
|----|------|----------|---------------|
| M1 | Empty abundances → all metrics = 0 | Standard robustness | No data → no diversity |
| M2 | Null abundances → all metrics = 0 | Standard robustness | Null safety |
| M3 | Single species → Shannon = 0 | Shannon theory | −1·ln(1) = 0 |
| M4 | Single species → Simpson = 1.0 | Simpson theory | 1² = 1 |
| M5 | Single species → InverseSimpson = 1.0 | Hill (1973) | 1/1 = 1 |
| M6 | Single species → Pielou = 0 | Convention | ln(1) = 0 → undefined |
| M7 | Single species → ObservedSpecies = 1 | Definition | Count of species |
| M8 | Two equal species → Shannon = ln(2) | Shannon formula | −2×(0.5×ln(0.5)) = ln(2) |
| M9 | Two equal species → Simpson = 0.5 | Simpson formula | 0.5² + 0.5² = 0.5 |
| M10 | Two equal species → Pielou = 1.0 | Evenness theory | Perfect evenness |
| M11 | Four equal species → Shannon = ln(4) | Shannon formula | −4×(0.25×ln(0.25)) |
| M12 | Four equal species → Simpson = 0.25 | Simpson formula | 4×0.25² = 0.25 |
| M13 | Four equal species → InverseSimpson = 4.0 | Hill (1973) | 1/0.25 = 4 |
| M14 | Zero abundances filtered out | Implementation | ln(0) undefined |
| M15 | Abundances normalized before calculation | Implementation | Sum may not be 1 |
| M16 | Shannon ≥ 0 always | Shannon theory | −Σ pᵢ ln(pᵢ) ≥ 0 for pᵢ ∈ (0,1] |
| M17 | Simpson ∈ [0, 1] | Simpson formula | pᵢ ∈ [0,1] → Σpᵢ² ∈ [0,1] |
| M18 | Pielou ∈ [0, 1] for S > 1 | Evenness bounds | H ≤ ln(S) → J ≤ 1 |
| M19 | Chao1 ≥ ObservedSpecies | Chao (1984) | S_Chao1 = S_obs + f₁²/(2·f₂) ≥ S_obs |
| M20 | Chao1 with count data: singletons + doubletons | Chao (1984) | f₁=2, f₂=1, S_obs=5 → 7 |
| M21 | Chao1 bias-corrected when f₂=0 | Chao (1984) | S_obs + f₁·(f₁−1)/2 |
| M22 | Chao1 = S_obs for proportional data | Chao (1984) | No integer counts → no singletons |

### SHOULD Tests

| ID | Test | Rationale |
|----|------|-----------|
| S1 | Highly uneven distribution → low Shannon, high Simpson | Diversity theory |
| S2 | Shannon increases with increasing richness (even distribution) | Diversity property |
| S3 | Simpson decreases with increasing richness (even distribution) | Dominance decreases |
| S4 | InverseSimpson = S for perfectly even distribution | Effective species count |

### COULD Tests

| ID | Test | Rationale |
|----|------|-----------|
| C1 | Large input (10000 taxa) processed efficiently | Scalability |
| C2 | Numerical stability with very small abundances | Floating-point edge case |

## Coverage Classification

### Test File

| File | Tests | Coverage |
|------|-------|----------|
| MetagenomicsAnalyzer_AlphaDiversity_Tests.cs | 29 tests (20 [Test] + 4 [TestCase] × S4 + 6 × invariants) | ✅ Complete (M1-M22, S1-S4, C1-C2) |

### Test → Spec Mapping

| Test | Covers | Assessment |
|------|--------|------------|
| `EmptyAbundances_ReturnsAllZeros` | M1 | ✅ Covered — Assert.Multiple for all 6 metrics |
| `NullAbundances_ReturnsAllZeros` | M2 | ✅ Covered — Assert.Multiple for all 6 metrics |
| `SingleSpecies_ReturnsCorrectMetrics` | M3-M7 | ✅ Covered — exact values with tolerance |
| `TwoEqualSpecies_ReturnsCorrectMetrics` | M8-M10 | ✅ Covered — exact ln(2), 0.5, 2.0, 1.0 |
| `FourEqualSpecies_ReturnsCorrectMetrics` | M11-M13 | ✅ Covered — exact ln(4), 0.25, 4.0, 1.0 |
| `EvenDistribution_InverseSimpsonEqualsSpeciesCount` | S4 | ✅ Covered — parametric (2, 5, 10, 100) |
| `HighlyUneven_LowDiversityHighDominance` | S1 | ✅ Covered — exact computed Shannon, Simpson, InverseSimpson, Pielou, Chao1 |
| `UnevenThreeSpecies_CalculatesCorrectly` | S1 extra | ✅ Covered — exact 3-species computation |
| `ZeroAbundances_FilteredOut` | M14 | ✅ Covered — exact values after filtering |
| `AllZeroAbundances_ReturnsZeros` | M14 | ✅ Covered — all-zero edge case |
| `UnnormalizedAbundances_Normalized` | M15 | ✅ Covered — exact values after normalization |
| `LargeTaxonCount_HandledCorrectly` | C1 | ✅ Covered — 10000 taxa, exact Shannon=ln(S), Simpson=1/S |
| `VerySmallAbundances_NumericallyStable` | C2 | ✅ Covered — exact computed Shannon, Simpson, InverseSimpson, Pielou |
| `CountDataWithSingletons_Chao1ExceedsObserved` | M20 | ✅ Covered — Chao1 = 7 for {50,30,1,1,2} |
| `CountDataNoDoubletons_Chao1BiasCorrected` | M21 | ✅ Covered — Chao1 = 7 (bias-corrected) |
| `CountDataNoSingletons_Chao1EqualsObserved` | M22 | ✅ Covered — Chao1 = S_obs when f₁=0 |
| `ProportionalData_Chao1EqualsObserved` | M22 | ✅ Covered — Chao1 = S_obs for fractional data |
| `VariousDistributions_SatisfyTheoreticalBounds` | M16-M19 | ✅ Covered — 6 distributions × bounds |
| `ShannonIncreasesWithRichness` | S2 | ✅ Covered — monotonic for S=1..16 |
| `SimpsonDecreasesWithRichness` | S3 | ✅ Covered — monotonic for S=1..16 |

### Test Organization

```
MetagenomicsAnalyzer_AlphaDiversity_Tests.cs
├── Empty and Null Input
│   ├── EmptyAbundances_ReturnsAllZeros (M1)
│   └── NullAbundances_ReturnsAllZeros (M2)
├── Single Species
│   └── SingleSpecies_ReturnsCorrectMetrics (M3-M7)
├── Even Distributions
│   ├── TwoEqualSpecies_ReturnsCorrectMetrics (M8-M10)
│   ├── FourEqualSpecies_ReturnsCorrectMetrics (M11-M13)
│   └── EvenDistribution_InverseSimpsonEqualsSpeciesCount (S4)
├── Uneven Distributions
│   ├── HighlyUneven_LowDiversityHighDominance (S1)
│   └── UnevenThreeSpecies_CalculatesCorrectly
├── Edge Cases
│   ├── ZeroAbundances_FilteredOut (M14)
│   ├── AllZeroAbundances_ReturnsZeros (M14)
│   ├── UnnormalizedAbundances_Normalized (M15)
│   ├── LargeTaxonCount_HandledCorrectly (C1)
│   ├── VerySmallAbundances_NumericallyStable (C2)
│   ├── CountDataWithSingletons_Chao1ExceedsObserved (M20)
│   ├── CountDataNoDoubletons_Chao1BiasCorrected (M21)
│   ├── CountDataNoSingletons_Chao1EqualsObserved (M22)
│   └── ProportionalData_Chao1EqualsObserved (M22)
└── Invariants
    ├── VariousDistributions_SatisfyTheoreticalBounds (M16-M19)
    ├── ShannonIncreasesWithRichness (S2)
    └── SimpsonDecreasesWithRichness (S3)
```

## Open Questions

None — all behavior is defined by authoritative sources.

## Decisions

1. Pielou evenness = 0 when S ≤ 1 — per Pielou (1966): ln(1) = 0 makes J = H/ln(S) undefined; J = 0 is the standard ecological convention
2. Natural logarithm (ln / nats) for Shannon index — per Shannon (1948), the most common base in ecology per Wikipedia Diversity Index
3. Chao1 uses the full formula from Chao (1984) for integer count data; returns S_obs for proportional data where singletons/doubletons are undefined
