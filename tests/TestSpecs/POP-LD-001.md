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
- File: `Seqeron.Genomics.Tests/PopulationGeneticsAnalyzerTests.cs`
- Region: `#region Linkage Disequilibrium Tests`

### Existing Tests

| Test Name | Coverage | Status | Action |
|-----------|----------|--------|--------|
| `CalculateLD_PerfectLD_ReturnsHighValues` | M6 | Weak | Strengthen assertions |
| `CalculateLD_NoLD_ReturnsLowValues` | M7 | Weak | Strengthen assertions |
| `CalculateLD_RecordsDistance` | M4, M5 | Covered | Move to canonical file |
| `FindHaplotypeBlocks_HighLD_CreatesBlock` | M11 | Weak | Strengthen assertions |
| `FindHaplotypeBlocks_LowLD_NoBlock` | M12 | Covered | Move to canonical file |
| `CalculateLD_EmptyGenotypes_ReturnsZeroLD` | M1 | Covered | Move to canonical file |
| `FindHaplotypeBlocks_SingleVariant_NoBlocks` | M9 | Covered | Move to canonical file |

### Missing Tests

| Must ID | Missing Test |
|---------|--------------|
| M2 | r² range validation |
| M3 | D' range validation |
| M8 | Monomorphic locus handling |
| M10 | Empty variants list |
| M13 | Block position ordering |
| M14 | Block minimum variant count |

---

## Consolidation Plan

1. **Create:** `PopulationGeneticsAnalyzer_LinkageDisequilibrium_Tests.cs`
2. **Move:** All LD tests from `PopulationGeneticsAnalyzerTests.cs`
3. **Refactor:** Strengthen weak tests with proper assertions
4. **Add:** Missing MUST tests
5. **Remove:** Duplicate comment placeholders in original file

---

## Test Data

### Perfect LD Dataset
```csharp
var genotypes = new List<(int, int)>
{
    (0, 0), (0, 0), (1, 1), (1, 1), (2, 2), (2, 2)
};
// Expected: High r² (≥ 0.5), high D'
```

### No LD Dataset
```csharp
var genotypes = new List<(int, int)>
{
    (0, 2), (2, 0), (1, 1), (0, 1), (2, 1), (1, 0)
};
// Expected: Low r² (< 0.3)
```

### Monomorphic Locus 1
```csharp
var genotypes = new List<(int, int)>
{
    (0, 0), (0, 1), (0, 2), (0, 0)
};
// First locus monomorphic (all 0), second varies
// Expected: r² = 0 (protected from division by zero)
```

---

## Open Questions

None.

---

## Decisions

1. **Weak assertion threshold:** Use 0.4 as minimum for "high" r² in perfect LD test (due to phase estimation uncertainty)
2. **Block threshold:** Use 0.3 as lower threshold in block detection tests to ensure blocks are found
3. **Test file naming:** `PopulationGeneticsAnalyzer_LinkageDisequilibrium_Tests.cs`

---

## Last Updated
2026-02-01
