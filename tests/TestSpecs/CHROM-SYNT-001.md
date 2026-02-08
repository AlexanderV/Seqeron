# TestSpec: CHROM-SYNT-001 - Synteny Analysis

**Test Unit ID:** CHROM-SYNT-001  
**Area:** Chromosome Analysis  
**Date:** 2026-02-01  
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
| M7 | FindSyntenyBlocks_GeneCountMatchesInput | GeneCount accuracy | Definition |
| M8 | FindSyntenyBlocks_CoordinatesSpanAllGenes | Block boundaries correctness | Definition |
| M9 | DetectRearrangements_Inversion_DetectsInversion | Core functionality: inversion detection | Wikipedia (Chromosomal rearrangement) |
| M10 | DetectRearrangements_Translocation_DetectsTranslocation | Core functionality: translocation detection | Wikipedia (Chromosomal rearrangement) |
| M11 | DetectRearrangements_EmptyInput_ReturnsEmpty | Edge case: empty input | Implementation |
| M12 | DetectRearrangements_SingleBlock_ReturnsEmpty | Edge case: no adjacent pairs | Definition |
| M13 | DetectRearrangements_NoRearrangements_ReturnsEmpty | Collinear genome returns empty | Definition |

### 2.2 SHOULD Tests (Recommended)

| ID | Test Name | Rationale |
|----|-----------|-----------|
| S1 | FindSyntenyBlocks_StrandIsValidChar | Invariant: strand is '+' or '-' |
| S2 | FindSyntenyBlocks_SequenceIdentityInRange | Invariant: 0 ≤ identity ≤ 1 |
| S3 | FindSyntenyBlocks_CoordinatesValid | Invariant: Start ≤ End |
| S4 | DetectRearrangements_TypeIsRecognizedValue | Invariant: valid type string |
| S5 | DetectRearrangements_Position1AlwaysSet | Invariant: non-null position |

### 2.3 COULD Tests (Nice to Have)

| ID | Test Name | Rationale |
|----|-----------|-----------|
| C1 | FindSyntenyBlocks_LargeDataset_CompletesReasonably | Performance baseline |
| C2 | FindSyntenyBlocks_CustomMaxGap_AffectsBlocks | Parameter behavior |

---

## 3. Test Data

### 3.1 Collinear Forward Block (M1, M7, M8)

```csharp
var orthologPairs = new List<(string, int, int, string, string, int, int, string)>
{
    ("chr1", 1000, 2000, "gene1", "chrA", 1000, 2000, "geneA"),
    ("chr1", 3000, 4000, "gene2", "chrA", 3000, 4000, "geneB"),
    ("chr1", 5000, 6000, "gene3", "chrA", 5000, 6000, "geneC"),
    ("chr1", 7000, 8000, "gene4", "chrA", 7000, 8000, "geneD"),
};
// Expected: 1 block, Strand='+', GeneCount=4
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
// Expected: Block with Strand='-'
```

### 3.3 Inversion Detection (M9)

```csharp
var blocks = new List<ChromosomeAnalyzer.SyntenyBlock>
{
    new("chr1", 1000, 50000, "chrA", 1000, 50000, '+', 10, 0.95),
    new("chr1", 60000, 100000, "chrA", 60000, 100000, '-', 8, 0.93)
};
// Expected: Inversion detected
```

### 3.4 Translocation Detection (M10)

```csharp
var blocks = new List<ChromosomeAnalyzer.SyntenyBlock>
{
    new("chr1", 1000, 50000, "chrA", 1000, 50000, '+', 10, 0.95),
    new("chr1", 60000, 100000, "chrB", 1000, 40000, '+', 8, 0.93)
};
// Expected: Translocation detected (chrA → chrB)
```

---

## 4. Consolidation Plan

### 4.1 Current State

Existing tests in `ChromosomeAnalyzerTests.cs`:
- `FindSyntenyBlocks_WithCollinearGenes_FindsBlocks` → Keep, enhance
- `FindSyntenyBlocks_WithInvertedGenes_DetectsInversion` → Keep, enhance
- `FindSyntenyBlocks_TooFewGenes_ReturnsEmpty` → Keep
- `DetectRearrangements_WithTranslocation_DetectsTranslocation` → Keep
- `DetectRearrangements_WithInversion_DetectsInversion` → Keep

### 4.2 Action Plan

1. **Create** new dedicated test file: `ChromosomeAnalyzer_Synteny_Tests.cs`
2. **Move** existing synteny tests from `ChromosomeAnalyzerTests.cs` to new file
3. ~~**Add** missing MUST tests (M4, M5, M6, M7, M8, M11, M12, M13)~~ ✅ All added
4. **Add** SHOULD tests (S1-S5)
5. **Remove** synteny region from general `ChromosomeAnalyzerTests.cs`
6. **Rename** tests to follow `Method_Scenario_ExpectedResult` convention

---

## 5. Invariants

### 5.1 FindSyntenyBlocks

- All blocks have `GeneCount >= minGenes`
- `Species1Start <= Species1End` and `Species2Start <= Species2End`
- `Strand` is either `'+'` or `'-'`
- `SequenceIdentity` is in range `[0, 1]`

### 5.2 DetectRearrangements

- `Type` is one of: "Inversion", "Translocation"
- `Position1` is always set
- For translocations, `Chromosome2` differs from the expected continuation chromosome

---

## 6. Audit Results

### 6.1 Coverage Assessment

| Test Area | Status | Notes |
|-----------|--------|-------|
| Forward collinearity | Covered | Existing test |
| Reverse collinearity | Covered | Existing test |
| minGenes threshold | Covered | Existing test |
| Empty input (FindSynteny) | ✅ Covered | Added |
| Boundary (exactly minGenes) | ✅ Covered | Added |
| Multi-chromosome pairs | ✅ Covered | Added |
| Inversion detection | ✅ Covered | Existing test |
| Translocation detection | ✅ Covered | Existing test |
| Empty input (DetectRearr) | ✅ Covered | Added |
| Single block | ✅ Covered | Added |
| No rearrangements | ✅ Covered | Added |
| Invariant tests | ✅ Covered | Added |

### 6.2 Test Quality

| Criterion | Status |
|-----------|--------|
| Naming convention | Needs update |
| Assert.Multiple usage | Needs addition |
| Edge case coverage | Incomplete |
| Invariant verification | ✅ Covered |

---

## 7. Open Questions / Decisions

1. **Q:** Should we test large datasets for performance?
   **A:** COULD test only; not required for DoD.

2. **Q:** Should deleted/duplicated region detection be tested?
   **A:** No; current implementation only detects inversions and translocations.

3. **Q:** Is the 0.9 placeholder identity acceptable?
   **A:** Yes; documented as implementation limitation.
