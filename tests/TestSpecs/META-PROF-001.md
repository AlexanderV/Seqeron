# META-PROF-001: Taxonomic Profile Test Specification

## Test Unit Information

| Field | Value |
|-------|-------|
| **ID** | META-PROF-001 |
| **Area** | Metagenomics |
| **Canonical Methods** | `MetagenomicsAnalyzer.GenerateTaxonomicProfile` |
| **Complexity** | O(n) where n = number of classifications |
| **Invariants** | 0 ≤ abundance ≤ 1; Σ(abundances) ≈ 1.0; ClassifiedReads ≤ TotalReads |
| **Status** | ☑ Complete |

## Methods Under Test

| Method | Class | Type | Test Depth |
|--------|-------|------|------------|
| `GenerateTaxonomicProfile(classifications)` | MetagenomicsAnalyzer | Canonical | Deep |

## Evidence Sources

1. **Wikipedia — Metagenomics:** Taxonomic profiling, relative abundance normalization
2. **Wikipedia — Relative abundance distribution:** Abundance computation in ecological contexts
3. **Shannon (1948):** Shannon diversity index formula
4. **Simpson (1949):** Simpson concentration index formula
5. **Segata et al. (2012):** MetaPhlAn paper on taxonomic profiling methodology

## Test Categories

### MUST Tests (Evidence-Backed)

| ID | Test | Evidence | Justification |
|----|------|----------|---------------|
| M1 | Empty input returns TotalReads=0, ClassifiedReads=0, empty abundances | Standard robustness | No input → no output |
| M2 | Single classified read produces abundance=1.0 for that taxon | Abundance formula | 1/1 = 1.0 |
| M3 | Unclassified reads excluded from ClassifiedReads count | MetaPhlAn: unclassified excluded | Normalization correctness |
| M4 | Unclassified reads excluded from abundance denominator | MetaPhlAn/Wikipedia | Abundance = count / classified_total |
| M5 | TotalReads = |input classifications| | Counting invariant | Basic consistency |
| M6 | ClassifiedReads = count(non-unclassified) | Definition | Filtering correctness |
| M7 | ClassifiedReads ≤ TotalReads | Bound invariant | Logical constraint |
| M8 | Sum of kingdom abundances ≈ 1.0 | Normalization invariant | Abundances are fractions |
| M9 | Shannon diversity ≥ 0 | Shannon formula (−Σp·ln(p)) | ln(p) ≤ 0 for p ≤ 1 |
| M10 | Simpson diversity ∈ [0, 1] | Simpson formula (Σp²) | p ∈ [0,1] |
| M11 | Single species → Shannon = 0 | Shannon theory | No uncertainty |
| M12 | Single species → Simpson = 1.0 | Simpson theory | Certainty |
| M13 | Empty strings filtered from rank-specific abundances | Implementation note | Clean output |
| M14 | Multiple taxa with equal counts → equal abundances | Uniformity test | Fair distribution |

### SHOULD Tests

| ID | Test | Rationale |
|----|------|-----------|
| S1 | Uniform distribution produces expected Shannon value | Formula verification |
| S2 | Uniform distribution produces expected Simpson value | Formula verification |
| S3 | All ranks have consistent total counts | Cross-rank consistency |
| S4 | High skew produces low Shannon, high Simpson | Diversity theory |

### COULD Tests

| ID | Test | Rationale |
|----|------|-----------|
| C1 | Large input processed efficiently | Scalability |
| C2 | Abundance values are exact (no floating point accumulation errors) | Numerical stability |

## Existing Test Audit

### Current Test Pool

| File | Tests | Coverage |
|------|-------|----------|
| MetagenomicsAnalyzerTests.cs | 3 GenerateTaxonomicProfile tests | M2, M3, M9-M10 (partial) |

### Existing Test Analysis

| Existing Test | Covers | Action |
|---------------|--------|--------|
| `GenerateTaxonomicProfile_CalculatesAbundances` | M2 (partial), M5 (partial) | Strengthen with invariants |
| `GenerateTaxonomicProfile_CalculatesDiversity` | M9, M10 (partial) | Add formula verification |
| `GenerateTaxonomicProfile_WithUnclassified_ExcludesFromAbundance` | M3, M4 | Keep, strengthen assertions |

### Missing Tests (All Closed)

| ID | Gap | Status |
|----|-----|--------|
| M1 | Empty input | ✅ Covered |
| M6 | ClassifiedReads count verification | ✅ Covered |
| M7 | Bound invariant | ✅ Covered |
| M8 | Sum normalization | ✅ Covered |
| M11 | Single species Shannon=0 | ✅ Covered |
| M12 | Single species Simpson=1.0 | ✅ Covered |
| M13 | Empty string filtering | ✅ Covered |
| M14 | Equal counts → equal abundances | ✅ Covered |

## Consolidation Plan

### Target Structure

Create `MetagenomicsAnalyzer_TaxonomicProfile_Tests.cs`:
- Migrate and consolidate existing GenerateTaxonomicProfile tests
- Add all missing MUST tests
- Use `Assert.Multiple` for related invariants
- Group tests by aspect (abundances, diversity, edge cases)

### Test File Structure

```
MetagenomicsAnalyzer_TaxonomicProfile_Tests.cs
├── #region Abundance Tests
│   ├── EmptyInput_ReturnsEmptyProfile [M1]
│   ├── SingleClassification_AbundanceIsOne [M2]
│   ├── MultipleClassifications_AbundancesSumToOne [M8]
│   ├── EqualCounts_ProduceEqualAbundances [M14]
│   └── EmptyRankValues_FilteredFromAbundance [M13]
├── #region Unclassified Handling
│   ├── UnclassifiedReads_ExcludedFromClassifiedCount [M3, M6]
│   └── UnclassifiedReads_ExcludedFromAbundanceDenominator [M4]
├── #region Read Count Invariants
│   ├── TotalReads_EqualsInputCount [M5]
│   └── ClassifiedReads_LessThanOrEqualToTotalReads [M7]
└── #region Diversity Metrics
    ├── SingleSpecies_ShannonIsZero [M11]
    ├── SingleSpecies_SimpsonIsOne [M12]
    ├── ShannonDiversity_NonNegative [M9]
    └── SimpsonDiversity_InZeroOneRange [M10]
```

## Open Questions

None. Algorithm behavior is well-defined by sources.

## Decisions

| Decision | Rationale |
|----------|-----------|
| Shannon uses natural log | Consistent with ecological convention and implementation |
| Simpson is concentration (Σp²) not diversity (1-Σp²) | Implementation uses concentration form |
| Empty input → all zero metrics | Safe default, no division by zero |
