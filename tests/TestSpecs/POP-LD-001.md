# TestSpec: POP-LD-001 - Linkage Disequilibrium

## Test Unit Information
- **ID:** POP-LD-001
- **Area:** Population Genetics (PopGen)
- **Status:** ☑ Complete
- **Evidence:** [POP-LD-001.md](../docs/Evidence/POP-LD-001.md)

---

## Methods Under Test

| Method | Class | Type | Test Depth |
|--------|-------|------|------------|
| `CalculateLD(var1, var2, genotypes, distance)` | PopulationGeneticsAnalyzer | Canonical | Deep |
| `FindHaplotypeBlocks(variants, ldThreshold)` | PopulationGeneticsAnalyzer | Canonical | Deep |

---

## Test Categories

### MUST Tests (Required for Complete Status)

| ID | Test Case | Source | Rationale |
|----|-----------|--------|-----------|
| M1 | CalculateLD returns r² = 0 for empty genotypes | Wikipedia (LD definition) | No data means no association |
| M2 | CalculateLD returns r² in range [0, 1] | Hill & Robertson (1968) | Mathematical invariant |
| M3 | CalculateLD returns D' in range [0, 1] | Lewontin (1964) | Normalized measure invariant |
| M4 | CalculateLD preserves variant IDs | API contract | Data integrity |
| M5 | CalculateLD preserves distance | API contract | Data integrity |
| M6 | CalculateLD with perfect LD returns high r² | Wikipedia (LD definition) | Complete association |
| M7 | CalculateLD with no LD returns low r² | Wikipedia (LD definition) | Random association |
| M8 | CalculateLD handles monomorphic locus (no division by zero) | Mathematical edge case | Denominator = 0 protection |
| M9 | FindHaplotypeBlocks returns empty for single variant | Gabriel (2002) | Block requires ≥2 variants |
| M10 | FindHaplotypeBlocks returns empty for < 2 variants | Gabriel (2002) | Block definition |
| M11 | FindHaplotypeBlocks creates block for high LD variants | Gabriel (2002) | Block detection |
| M12 | FindHaplotypeBlocks returns empty for low LD variants | Gabriel (2002) | No block when LD below threshold |
| M13 | FindHaplotypeBlocks block.Start ≤ block.End | Invariant | Position ordering |
| M14 | FindHaplotypeBlocks block contains ≥2 variants | Gabriel (2002) | Block definition |

### SHOULD Tests (Recommended)

| ID | Test Case | Source | Rationale |
|----|-----------|--------|-----------|
| S1 | CalculateLD with single genotype pair | Statistical validity | Minimum sample |
| S2 | CalculateLD with all homozygous major (0,0) | No variation | Edge case |
| S3 | CalculateLD with all homozygous minor (2,2) | No variation | Edge case |
| S4 | FindHaplotypeBlocks respects ldThreshold parameter | API contract | Threshold customization |
| S5 | FindHaplotypeBlocks orders variants by position | Gabriel (2002) | Correct block boundaries |
| S6 | FindHaplotypeBlocks handles multiple blocks | Gabriel (2002) | Complex genome regions |

### COULD Tests (Optional Enhancements)

| ID | Test Case | Source | Rationale |
|----|-----------|--------|-----------|
| C1 | CalculateLD with large sample size | Performance | Scalability |
| C2 | FindHaplotypeBlocks with many variants | Performance | Scalability |
| C3 | CalculateLD with mixed heterozygotes | Phasing complexity | Algorithm accuracy |

---

## Existing Tests Audit

### Current Location
- File: `Seqeron.Genomics.Tests/PopulationGeneticsAnalyzer_LinkageDisequilibrium_Tests.cs`

### Test Coverage

| Test Name | Covers | Status |
|-----------|--------|--------|
| `CalculateLD_EmptyGenotypes_ReturnsZeroLD` | M1 | ✅ Covered |
| `CalculateLD_PreservesVariantIdsAndDistance` | M4, M5 | ✅ Covered |
| `CalculateLD_PerfectLD_ReturnsExactValues` | M6 | ✅ Covered (r²=1.0, D'=1.0) |
| `CalculateLD_NoLD_ReturnsZeroValues` | M7 | ✅ Covered (r²=0.0, D'=0.0) |
| `CalculateLD_RSquared_AlwaysInValidRange` | M2 | ✅ Covered |
| `CalculateLD_DPrime_AlwaysInValidRange` | M3 | ✅ Covered |
| `CalculateLD_MonomorphicFirstLocus_ReturnsZeroRSquared` | M8 | ✅ Covered |
| `CalculateLD_MonomorphicSecondLocus_ReturnsZeroRSquared` | M8 | ✅ Covered |
| `CalculateLD_SingleGenotypePair_ReturnsValidResult` | S1 | ✅ Covered |
| `CalculateLD_AllHomozygousMajor_ReturnsZeroRSquared` | S2 | ✅ Covered |
| `CalculateLD_AllHomozygousMinor_ReturnsZeroRSquared` | S3 | ✅ Covered |
| `FindHaplotypeBlocks_SingleVariant_ReturnsNoBlocks` | M9 | ✅ Covered |
| `FindHaplotypeBlocks_EmptyVariants_ReturnsNoBlocks` | M10 | ✅ Covered |
| `FindHaplotypeBlocks_HighLD_CreatesBlock` | M11 | ✅ Covered |
| `FindHaplotypeBlocks_LowLD_ReturnsNoBlocks` | M12 | ✅ Covered |
| `FindHaplotypeBlocks_BlockPositions_StartLessThanOrEqualToEnd` | M13 | ✅ Covered |
| `FindHaplotypeBlocks_EachBlock_ContainsAtLeastTwoVariants` | M14 | ✅ Covered |
| `FindHaplotypeBlocks_ThresholdParameter_AffectsBlockFormation` | S4 | ✅ Covered |
| `FindHaplotypeBlocks_UnorderedInput_OrdersByPosition` | S5 | ✅ Covered |
| `FindHaplotypeBlocks_MultipleBlocks_DetectsAll` | S6 | ✅ Covered |

---

## Test Data

### Perfect LD Dataset
```csharp
var genotypes = new List<(int, int)>
{
    (0, 0), (0, 0), (1, 1), (1, 1), (2, 2), (2, 2)
};
// Identical genotype vectors → r² = 1.0, D' = 1.0
// Math: mean=1, Cov=Var=2/3 → r²=1.0; D=1/3, D_max=0.25, clamped → D'=1.0
```

### No LD Dataset (3×3 balanced design)
```csharp
var genotypes = new List<(int, int)>
{
    (0, 0), (0, 1), (0, 2),
    (1, 0), (1, 1), (1, 2),
    (2, 0), (2, 1), (2, 2)
};
// Every X₂ value appears equally within each X₁ level → Cov = 0 → r² = 0, D' = 0
```

### Monomorphic Locus 1
```csharp
var genotypes = new List<(int, int)>
{
    (0, 0), (0, 1), (0, 2), (0, 0)
};
// First locus monomorphic (all 0), Var(X₁) = 0
// Expected: r² = 0, D' = 0 (protected from division by zero)
```

---

## Deviations and Assumptions

**None.** All tests and implementation are grounded in the authoritative sources listed in the Evidence document.

- r² is computed as the squared Pearson correlation of genotype values (0, 1, 2). From Wikipedia (LD for diploid frequencies): the diploid correlation R_AB equals the haplotype-level r_AB.
- D is estimated from the diploid genotype covariance: D = Cov(X₁,X₂)/2, per Wikipedia (LD for diploid frequencies).
- D' is normalized per Lewontin (1964): D' = |D| / D_max, clamped to [0, 1].
- Haplotype block detection uses adjacent-pair r² threshold (simplified Gabriel et al. 2002). Default threshold: 0.7.

---

## Last Updated
2026-03-08
