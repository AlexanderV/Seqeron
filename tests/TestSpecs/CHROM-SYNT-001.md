# TestSpec: CHROM-SYNT-001 - Synteny Analysis

**Test Unit ID:** CHROM-SYNT-001  
**Area:** Chromosome Analysis  
**Date:** 2026-03-08  
**Status:** Complete

---

## 1. Scope

### 1.1 Methods Under Test

| Method | Class | Type | Testing Strategy |
|--------|-------|------|------------------|
| `FindSyntenyBlocks(orthologPairs, minGenes, maxGap)` | ChromosomeAnalyzer | Canonical | Deep testing |
| `DetectRearrangements(syntenyBlocks)` | ChromosomeAnalyzer | Canonical | Deep testing |

### 1.2 Supporting Types

- `SyntenyBlock` record: Species1Chromosome, Species1Start, Species1End, Species2Chromosome, Species2Start, Species2End, Strand, GeneCount, SequenceIdentity
- `ChromosomalRearrangement` record: Type, Chromosome1, Position1, Chromosome2, Position2, Size, Description

---

## 2. Test Classification

### 2.1 MUST Tests (Required for Definition of Done)

| ID | Test Name | Rationale | Source |
|----|-----------|-----------|--------|
| M1 | FindSyntenyBlocks_CollinearForward_ReturnsBlockWithPlusStrand | Core functionality: detect forward collinearity | Wikipedia (Synteny) |
| M2 | FindSyntenyBlocks_CollinearReverse_ReturnsBlockWithMinusStrand | Core functionality: detect inverted blocks | Wikipedia (Synteny) |
| M3 | FindSyntenyBlocks_TooFewGenes_ReturnsEmpty | minGenes threshold enforcement | Definition |
| M4 | FindSyntenyBlocks_EmptyInput_ReturnsEmpty | Edge case: empty input | Implementation |
| M5 | FindSyntenyBlocks_ExactlyMinGenes_ReturnsBlock | Boundary: exactly minGenes | Definition |
| M6 | FindSyntenyBlocks_MultipleChromosomePairs_SeparateBlocks | Chromosome pair grouping | Definition |
| M9 | DetectRearrangements_Inversion_DetectsInversion | Core functionality: inversion detection | Wikipedia (Chromosomal rearrangement) |
| M10 | DetectRearrangements_Translocation_DetectsTranslocation | Core functionality: translocation detection | Wikipedia (Chromosomal rearrangement) |
| M11 | DetectRearrangements_EmptyInput_ReturnsEmpty | Edge case: empty input | Implementation |
| M12 | DetectRearrangements_SingleBlock_ReturnsEmpty | Edge case: no adjacent pairs | Definition |
| M13 | DetectRearrangements_NoRearrangements_ReturnsEmpty | Collinear genome returns empty | Definition |
| M14 | DetectRearrangements_Deletion_DetectsDeletion | Core functionality: deletion detection | Wikipedia (Chromosomal rearrangement) |
| M15 | DetectRearrangements_Duplication_DetectsDuplication | Core functionality: duplication detection | Wikipedia (Chromosomal rearrangement) |
| M16 | FindSyntenyBlocks_GapExceedsMaxGap_SplitsIntoSeparateBlocks | maxGap parameter splits collinear runs | MCScanX (Wang et al. 2012) |

### 2.2 SHOULD Tests (Recommended)

| ID | Test Name | Rationale |
|----|-----------|-----------|
| S1 | FindSyntenyBlocks_StrandIsValidChar | Invariant: strand is '+' or '-' |
| S2 | FindSyntenyBlocks_SequenceIdentityIsNaN | Invariant: NaN when not computable |
| S3 | FindSyntenyBlocks_CoordinatesValid | Invariant: Start ≤ End |
| S4 | DetectRearrangements_TypeIsRecognizedValue | Invariant: valid type string (Inversion/Translocation/Deletion/Duplication) |
| S5 | DetectRearrangements_Position1AlwaysSet | Invariant: non-null position |

### 2.3 COULD Tests (Nice to Have)

| ID | Test Name | Rationale |
|----|-----------|-----------|
| C1 | FindSyntenyBlocks_LargeDataset_CompletesReasonably | Performance baseline |

---

## 3. Test Data

### 3.1 Collinear Forward Block (M1)

```csharp
var orthologPairs = new List<(string, int, int, string, string, int, int, string)>
{
    ("chr1", 1000, 2000, "gene1", "chrA", 1000, 2000, "geneA"),
    ("chr1", 3000, 4000, "gene2", "chrA", 3000, 4000, "geneB"),
    ("chr1", 5000, 6000, "gene3", "chrA", 5000, 6000, "geneC"),
    ("chr1", 7000, 8000, "gene4", "chrA", 7000, 8000, "geneD"),
};
// Expected: 1 block, Strand='+', GeneCount=4, Species1: 1000-8000, Species2: 1000-8000
```

### 3.2 Inverted Block (M2)

```csharp
var orthologPairs = new List<(string, int, int, string, string, int, int, string)>
{
    ("chr1", 1000, 2000, "gene1", "chrA", 8000, 9000, "geneA"),
    ("chr1", 3000, 4000, "gene2", "chrA", 6000, 7000, "geneB"),
    ("chr1", 5000, 6000, "gene3", "chrA", 4000, 5000, "geneC"),
    ("chr1", 7000, 8000, "gene4", "chrA", 2000, 3000, "geneD"),
};
// Expected: Block with Strand='-', GeneCount=4, Species2: 2000-9000
```

### 3.3 Inversion Detection (M9)

```csharp
var blocks = new List<ChromosomeAnalyzer.SyntenyBlock>
{
    new("chr1", 1000, 50000, "chrA", 1000, 50000, '+', 10, 0.95),
    new("chr1", 60000, 100000, "chrA", 60000, 100000, '-', 8, 0.93)
};
// Expected: Inversion detected, Position1=50000, Position2=60000, Size=10000
```

### 3.4 Translocation Detection (M10)

```csharp
var blocks = new List<ChromosomeAnalyzer.SyntenyBlock>
{
    new("chr1", 1000, 50000, "chrA", 1000, 50000, '+', 10, 0.95),
    new("chr1", 60000, 100000, "chrB", 1000, 40000, '+', 8, 0.93)
};
// Expected: Translocation detected (chrA → chrB), Position1=50000, Chromosome2="chrB"
```

### 3.5 Deletion Detection (M14)

```csharp
var blocks = new List<ChromosomeAnalyzer.SyntenyBlock>
{
    new("chr1", 1000, 50000, "chrA", 1000, 50000, '+', 10, 0.95),
    new("chr1", 150000, 200000, "chrA", 55000, 100000, '+', 8, 0.93)
};
// gap1=100000, gap2=5000
// Asymmetric gap → Deletion, Size = 100000-5000 = 95000
```

### 3.6 Duplication Detection (M15)

```csharp
var blocks = new List<ChromosomeAnalyzer.SyntenyBlock>
{
    new("chr1", 1000, 50000, "chrA", 1000, 50000, '+', 10, 0.95),
    new("chr1", 20000, 70000, "chrA", 200000, 250000, '+', 8, 0.93)
};
// Overlapping species 1 region: 20000-50000
// Different species 2 locations → Duplication, Size = 30000
```

---

## 4. Consolidation Plan

### 4.1 Current State

All synteny tests consolidated in `ChromosomeAnalyzer_Synteny_Tests.cs` (19 tests).
`ChromosomeAnalyzerTests.cs` contains only a redirect comment.

### 4.2 Changes (2026-03-08)

| Action | Tests | Rationale |
|--------|-------|-----------|
| ⚠ Weak → Strengthened | M1, M2, M5, M6, M9, M10, M14, M15 | Replaced range/permissive assertions with exact hand-calculated values |
| 🔁 Duplicate → Removed | M7 (GeneCountMatchesInput), M8 (CoordinatesSpanAllGenes) | Subsumed by strengthened M1 which asserts exact GeneCount=4, coordinates 1000-8000 |
| ❌ Missing → Implemented | M16 (GapExceedsMaxGap_SplitsIntoSeparateBlocks) | maxGap parameter behavior was untested |

**Total:** 8 strengthened, 2 removed, 1 added

---

## 5. Invariants

### 5.1 FindSyntenyBlocks

- All blocks have `GeneCount >= minGenes`
- `Species1Start <= Species1End` and `Species2Start <= Species2End`
- `Strand` is either `'+'` or `'-'`
- `SequenceIdentity` is `NaN` (not computable from coordinate-only input per MCScanX)

### 5.2 DetectRearrangements

- `Type` is one of: "Inversion", "Translocation", "Deletion", "Duplication"
- `Position1` is always set
- For translocations, `Chromosome2` differs from the expected continuation chromosome

---

## 6. Audit Results

### 6.1 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1: Forward collinearity | ✅ Covered | Exact: count=1, GeneCount=4, strand='+', coordinates 1000-8000 both species |
| M2: Reverse collinearity | ✅ Covered | Exact: count=1, GeneCount=4, strand='-', Species2: 2000-9000 |
| M3: Below minGenes | ✅ Covered | Asserts empty |
| M4: Empty input (Find) | ✅ Covered | Asserts empty |
| M5: Exactly minGenes | ✅ Covered | Exact: count=1, GeneCount=3, strand='+', coordinates 1000-6000 |
| M6: Multi-chromosome blocks | ✅ Covered | Exact: 2 blocks, each with correct chr pair, GeneCount=3 |
| M7: GeneCount accuracy | 🔁 Removed | Subsumed by strengthened M1 (GeneCount=4 asserted) |
| M8: Coordinates span | 🔁 Removed | Subsumed by strengthened M1 (coordinates 1000-8000 asserted) |
| M9: Inversion detection | ✅ Covered | Exact: count=1, Position1=50000, Position2=60000, Size=10000 |
| M10: Translocation detection | ✅ Covered | Exact: count=1, Chromosome2="chrB", Position1=50000, Position2=1000, Size=null |
| M11: Empty input (Detect) | ✅ Covered | Asserts empty |
| M12: Single block | ✅ Covered | Asserts empty |
| M13: No rearrangements | ✅ Covered | 3 collinear blocks → empty |
| M14: Deletion detection | ✅ Covered | Exact: count=1, Size=95000 (gap1-gap2=100000-5000) |
| M15: Duplication detection | ✅ Covered | Exact: count=1, Size=30000 (overlap 20000-50000) |
| M16: maxGap split | ✅ Covered | New: 6 genes, 3MB gap, maxGap=2 → 2 blocks of 3 genes |
| S1: Strand invariant | ✅ Covered | All blocks strand ∈ {'+', '-'} |
| S2: Identity NaN | ✅ Covered | SequenceIdentity = NaN (coordinate-only input) |
| S3: Coordinates valid | ✅ Covered | Start ≤ End for both species |
| S4: Type invariant | ✅ Covered | All types ∈ {Inversion, Translocation, Deletion, Duplication} |
| S5: Position1 set | ✅ Covered | Position1 > 0 for all rearrangements |

### 6.2 Final State

| File | Role | Test Count |
|------|------|------------|
| `ChromosomeAnalyzer_Synteny_Tests.cs` | Canonical | 19 |
| Invariant tests | ✅ Covered | Added |

### 6.2 Test Quality

| Criterion | Status |
|-----------|--------|
| Naming convention | ✅ Method_Scenario_ExpectedResult |
| Assert.Multiple usage | ✅ Used where appropriate |
| Edge case coverage | ✅ Complete |
| Invariant verification | ✅ Covered |
