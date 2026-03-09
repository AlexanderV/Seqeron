# META-BIN-001: Genome Binning - Test Specification

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | META-BIN-001 |
| **Algorithm** | Genome Binning (MAG assembly) |
| **Canonical Method** | `MetagenomicsAnalyzer.BinContigs(...)` |
| **Complexity** | O(n × k × i) |
| **Evidence Document** | [META-BIN-001-Evidence.md](../docs/Evidence/META-BIN-001-Evidence.md) |

## Method Under Test

```csharp
public static IEnumerable<GenomeBin> BinContigs(
    IEnumerable<(string ContigId, string Sequence, double Coverage)> contigs,
    int numBins = 10,
    double minBinSize = 500000,
    double expectedGenomeSize = 4_000_000)
```

### Output Record

```csharp
public readonly record struct GenomeBin(
    string BinId,
    IReadOnlyList<string> ContigIds,
    double TotalLength,
    double GcContent,
    double Coverage,
    double Completeness,
    double Contamination,
    string PredictedTaxonomy);
```

## Test Categories

### Must Tests (Required for DoD)

| ID | Test Name | Rationale | Source |
|----|-----------|-----------|--------|
| M1 | EmptyInput_ReturnsEmpty | Empty contigs must yield empty result | Standard edge case |
| M2 | SingleContig_BelowMinSize_ReturnsEmpty | Contigs below minBinSize excluded | Implementation behavior |
| M3 | ValidContigs_ReturnsNonEmptyBins | Non-trivial input must produce bins | Wikipedia |
| M4 | BinIds_AreUnique | Each bin must have unique ID | Invariant |
| M5 | Completeness_InValidRange | Must be [0, 100]; verified exact formula: min(totalLength/expectedGenomeSize×100, 100) | CheckM standard |
| M6 | Contamination_InValidRange | Must be [0, 100]; verified: uniform GC→0, max-variance GC→100 | CheckM standard |
| M7 | GcContent_InValidRange | Must be [0, 1]; verified exact values for pure high-GC (1.0) and low-GC (0.0) bins | Mathematical invariant |
| M8 | TotalLength_EqualsContigLengthSum | Length accounting correct | Mathematical invariant |
| M9 | MinBinSize_FiltersSmallBins | Bins below threshold excluded | Implementation spec |
| M10 | HighGcVsLowGc_SeparatesIntoDifferentBins | K-means separates extreme GC populations | Wikipedia, TETRA paper |
| M11 | ContigIds_PreservedInBins | All binned contigs traceable; bins disjoint | Data integrity |
| M12 | Coverage_PreservedInBins | Coverage = arithmetic mean of contig coverages; verified exact value | Implementation spec |

### Should Tests (Recommended)

| ID | Test Name | Rationale | Source |
|----|-----------|-----------|--------|
| S1 | MultipleBins_FromDiverseInput | Diverse contigs create multiple bins | Wikipedia |
| S2 | NumBins_AffectsBinCount | k-means with k clusters produces at most k bins; verified both limits | Implementation |
| S3 | LargeDataset_CompletesWithinTimeout | Performance acceptable | Practical requirement |

### Could Tests (Optional)

| ID | Test Name | Rationale |
|----|-----------|-----------|
| C1 | AllSameGc_ClustersTogether | Uniform features converge to few bins |
| C2 | VeryShortSequences_HandledGracefully | Robustness |

## Invariants

1. **Empty Preservation**: `BinContigs([]) → []`
2. **Unique IDs**: `∀ bin1, bin2 ∈ result: bin1.BinId ≠ bin2.BinId`
3. **Metric Ranges**:
   - Completeness ∈ [0, 100]
   - Contamination ∈ [0, 100]
   - GcContent ∈ [0, 1]
4. **Size Filtering**: `∀ bin ∈ result: bin.TotalLength ≥ minBinSize`
5. **Length Accuracy**: `bin.TotalLength = Σ len(contig) for contig in bin`

## Test Data Design

### High GC Contigs (for separation test)
- Sequence: 70% G/C content (e.g., `GGGGGGGGCC` repeated)
- Expected: Group together in high-GC bin

### Low GC Contigs (for separation test)
- Sequence: 30% G/C content (e.g., `AAAATTTTGC` repeated)
- Expected: Group together in low-GC bin

### Mixed Contigs
- Combination of high and low GC
- Expected: Separate into different bins

## Existing Test Status

| Location | Status | Action |
|----------|--------|--------|
| MetagenomicsAnalyzer_GenomeBinning_Tests.cs | All Must/Should/Could tests implemented | Complete |
| MetagenomicsAnalyzerTests.cs | Smoke test removed (duplicate of M3) | Cleaned |
| MetagenomicsSnapshotTests.cs | Snapshot test for deterministic output | Independent |

## Coverage Classification Summary

| Classification | Count | Details |
|---------------|-------|----------|
| ✅ Covered | 11 | M1, M2, M4, M8, M9, M10, S1, S3, C1, C2, extra:ZeroCoverage |
| ⚠ Strengthened | 7 | M3 (Not.Empty), M5 (exact formula), M6 (exact formula), M7 (exact GC), M11 (disjoint), M12 (exact average), S2 (both limits) |
| 🔁 Removed | 1 | BinContigs_BasicFunctionality_ReturnsResults (subsumed by M3) |
| ❌ Missing | 0 | None |

## Audit Notes

- Current smoke tests adequate for basic functionality
- Need comprehensive invariant testing
- Need GC-based separation validation
- Need quality metric range validation

## Consolidation Plan

All tests implemented and verified in `MetagenomicsAnalyzer_GenomeBinning_Tests.cs`.
Duplicate smoke test in `MetagenomicsAnalyzerTests.cs` has been removed.
Snapshot test in `MetagenomicsSnapshotTests.cs` remains independent.

---

**Specification Version:** 2.0
**Created:** 2026-02-04
**Updated:** 2026-03-09
**Status:** Complete — all tests implemented, strengthened, and verified
