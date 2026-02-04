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
    double minBinSize = 500000)
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
| M3 | ValidContigs_ReturnsNonNullBins | Basic functionality | Wikipedia |
| M4 | BinIds_AreUnique | Each bin must have unique ID | Invariant |
| M5 | Completeness_InValidRange | Must be [0, 100] | CheckM standard |
| M6 | Contamination_InValidRange | Must be [0, 100] | CheckM standard |
| M7 | GcContent_InValidRange | Must be [0, 1] | Mathematical invariant |
| M8 | TotalLength_EqualsContigLengthSum | Length accounting correct | Mathematical invariant |
| M9 | MinBinSize_FiltersSmallBins | Bins below threshold excluded | Implementation spec |
| M10 | HighGcVsLowGc_SeparatesIntoDifferentBins | GC-based separation works | TETRA paper |
| M11 | ContigIds_PreservedInBins | All binned contigs traceable | Data integrity |
| M12 | Coverage_PreservedInBins | Coverage values averaged correctly | Implementation spec |

### Should Tests (Recommended)

| ID | Test Name | Rationale | Source |
|----|-----------|-----------|--------|
| S1 | MultipleBins_FromDiverseInput | Diverse contigs create multiple bins | Wikipedia |
| S2 | NumBins_AffectsBinCount | Parameter respected | Implementation |
| S3 | LargeDataset_CompletesWithinTimeout | Performance acceptable | Practical requirement |

### Could Tests (Optional)

| ID | Test Name | Rationale |
|----|-----------|-----------|
| C1 | AllSameGc_ClustersTogether | Edge case validation |
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
| MetagenomicsAnalyzerTests.cs | Smoke tests only | Keep as reference, expand in dedicated file |

## Audit Notes

- Current smoke tests adequate for basic functionality
- Need comprehensive invariant testing
- Need GC-based separation validation
- Need quality metric range validation

## Consolidation Plan

1. Create new dedicated test file: `MetagenomicsAnalyzer_GenomeBinning_Tests.cs`
2. Keep existing smoke tests in `MetagenomicsAnalyzerTests.cs` as integration reference
3. Implement all Must tests with evidence-based rationale
4. Remove duplicates (none found)

---

**Specification Version:** 1.0
**Created:** 2026-02-04
**Status:** Ready for Implementation
