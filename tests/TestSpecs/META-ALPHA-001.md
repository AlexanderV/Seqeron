# META-ALPHA-001: Alpha Diversity Test Specification

## Test Unit Information

| Field | Value |
|-------|-------|
| **ID** | META-ALPHA-001 |
| **Area** | Metagenomics |
| **Canonical Methods** | `MetagenomicsAnalyzer.CalculateAlphaDiversity` |
| **Complexity** | O(n) where n = number of taxa |
| **Invariants** | H ≥ 0; λ ∈ [0,1]; J ∈ [0,1]; S ≥ 0 |

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
| M19 | Chao1 ≥ ObservedSpecies | Chao theory | Estimates ≥ observed |

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

## Existing Test Audit

### Current Test Pool

| File | Tests | Coverage |
|------|-------|----------|
| MetagenomicsAnalyzerTests.cs | 5 CalculateAlphaDiversity tests | Partial (M1, M3, M4, M8-M10) |

### Existing Test Analysis

| Existing Test | Covers | Assessment | Action |
|---------------|--------|------------|--------|
| `CalculateAlphaDiversity_SingleSpecies_LowDiversity` | M3, M4, M7 | Weak assertions | Strengthen with exact values |
| `CalculateAlphaDiversity_EvenDistribution_HighDiversity` | M8, M9 (partial) | Uses GreaterThan | Replace with exact formula values |
| `CalculateAlphaDiversity_UnevenDistribution_CalculatesCorrectly` | S1 | Weak bounds | Strengthen with exact calculation |
| `CalculateAlphaDiversity_EmptyAbundances_ReturnsZero` | M1 | Missing full invariant check | Add Assert.Multiple |
| `CalculateAlphaDiversity_CalculatesPielouEvenness` | M10 | Good | Keep |

### Missing Tests

| ID | Gap | Priority |
|----|-----|----------|
| M2 | Null input | Must add |
| M5, M6 | Single species InverseSimpson, Pielou | Must add |
| M11-M13 | Four equal species exact values | Must add |
| M14 | Zero abundance filtering | Must add |
| M15 | Normalization | Must add |
| M16-M18 | Invariant bounds | Must add |
| M19 | Chao1 bound | Must add |
| S2-S4 | Richness scaling properties | Should add |

### Duplicate/Weak Tests

| Test | Issue | Resolution |
|------|-------|------------|
| EvenDistribution test | Uses vague GreaterThan | Replace with exact ln(4) |
| UnevenDistribution test | Bounds too wide | Tighten with calculated values |

## Consolidation Plan

### Target Structure

Create `MetagenomicsAnalyzer_AlphaDiversity_Tests.cs`:
- Extract and strengthen existing tests from MetagenomicsAnalyzerTests.cs
- Remove Alpha Diversity tests from MetagenomicsAnalyzerTests.cs
- Add all missing MUST tests
- Group by index type using regions
- Use `Assert.Multiple` for multi-invariant checks

### Test Organization

```
MetagenomicsAnalyzer_AlphaDiversity_Tests.cs
├── Empty and Null Input
│   ├── EmptyAbundances_ReturnsAllZeros
│   └── NullAbundances_ReturnsAllZeros
├── Single Species
│   └── SingleSpecies_ReturnsCorrectMetrics
├── Even Distributions
│   ├── TwoEqualSpecies_ReturnsCorrectMetrics
│   └── FourEqualSpecies_ReturnsCorrectMetrics
├── Uneven Distributions
│   ├── HighlyUneven_LowDiversityHighDominance
│   └── UnevenDistribution_CorrectPielouEvenness
├── Edge Cases
│   ├── ZeroAbundances_FilteredOut
│   ├── UnnormalizedAbundances_Normalized
│   └── LargeTaxonCount_HandledCorrectly
└── Invariants
    └── AllResults_SatisfyTheoreticalBounds
```

## Open Questions

None — all behavior is defined by authoritative sources.

## Decisions

1. Pielou evenness = 0 when S ≤ 1 (implementation convention)
2. Chao1 = ObservedSpecies (simplified, documented as ASSUMPTION)
3. Use natural logarithm (ln) for Shannon index
