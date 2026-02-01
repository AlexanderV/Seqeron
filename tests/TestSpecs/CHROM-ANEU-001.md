# TestSpec: CHROM-ANEU-001 - Aneuploidy Detection

**Test Unit ID:** CHROM-ANEU-001  
**Area:** Chromosome Analysis  
**Date:** 2026-02-01  
**Status:** Complete

---

## 1. Scope

### 1.1 Methods Under Test

| Method | Class | Type | Testing Strategy |
|--------|-------|------|------------------|
| `DetectAneuploidy(depthData, medianDepth, binSize)` | ChromosomeAnalyzer | Canonical | Deep testing |
| `IdentifyWholeChromosomeAneuploidy(cnStates, minFraction)` | ChromosomeAnalyzer | Classification | Deep testing |

### 1.2 Supporting Types

- `CopyNumberState` record: Chromosome, Start, End, CopyNumber, LogRatio, Confidence

---

## 2. Test Classification

### 2.1 MUST Tests (Required for Definition of Done)

| ID | Test Name | Rationale | Source |
|----|-----------|-----------|--------|
| M1 | DetectAneuploidy_NormalDepth_ReturnsCopyNumber2 | Verify diploid detection | Wikipedia |
| M2 | DetectAneuploidy_Trisomy_ReturnsCopyNumber3 | Verify trisomy detection (1.5× depth) | Wikipedia |
| M3 | DetectAneuploidy_Monosomy_ReturnsCopyNumber1 | Verify monosomy detection (0.5× depth) | Wikipedia |
| M4 | DetectAneuploidy_EmptyInput_ReturnsEmpty | Edge case: empty input | Implementation |
| M5 | DetectAneuploidy_ZeroMedianDepth_ReturnsEmpty | Edge case: prevents division by zero | Implementation |
| M6 | DetectAneuploidy_Nullisomy_ReturnsCopyNumber0 | Verify CN=0 detection (0× depth) | Wikipedia |
| M7 | DetectAneuploidy_Tetrasomy_ReturnsCopyNumber4 | Verify CN=4 detection (2× depth) | Wikipedia |
| M8 | DetectAneuploidy_MultipleChromosomes_GroupsCorrectly | Verify chromosome grouping | Implementation |
| M9 | DetectAneuploidy_BinAggregation_AveragesDepth | Verify binning logic | Implementation |
| M10 | DetectAneuploidy_CopyNumberClamped_MaximumIs10 | Verify upper bound | Implementation |
| M11 | IdentifyWholeChromosomeAneuploidy_Trisomy_IdentifiesTrisomy | Classification correctness | Wikipedia |
| M12 | IdentifyWholeChromosomeAneuploidy_Monosomy_IdentifiesMonosomy | Classification correctness | Wikipedia |
| M13 | IdentifyWholeChromosomeAneuploidy_Normal_ReturnsEmpty | Normal chromosome not flagged | Definition |
| M14 | IdentifyWholeChromosomeAneuploidy_Nullisomy_IdentifiesNullisomy | Classification correctness | Wikipedia |
| M15 | IdentifyWholeChromosomeAneuploidy_Tetrasomy_IdentifiesTetrasomy | Classification correctness | Wikipedia |
| M16 | IdentifyWholeChromosomeAneuploidy_MixedCN_RequiresMinFraction | Mosaicism threshold | Implementation |

### 2.2 SHOULD Tests (Recommended)

| ID | Test Name | Rationale |
|----|-----------|-----------|
| S1 | DetectAneuploidy_Confidence_InRangeZeroToOne | Output invariant |
| S2 | DetectAneuploidy_LogRatio_MatchesFormula | Mathematical correctness |
| S3 | DetectAneuploidy_OutputOrdered_ByPosition | Ordering invariant |
| S4 | IdentifyWholeChromosomeAneuploidy_HighCopyNumber_FormatsCorrectly | CN > 4 formatting |

### 2.3 COULD Tests (Nice to Have)

| ID | Test Name | Rationale |
|----|-----------|-----------|
| C1 | DetectAneuploidy_LargeBinSize_ReducesOutput | Parameter behavior |
| C2 | IdentifyWholeChromosomeAneuploidy_CustomMinFraction_Works | Parameter behavior |

---

## 3. Test Data

### 3.1 Standard Depth Values

```
medianDepth = 30.0 (arbitrary baseline)

Normal (CN=2):    depth = 30.0    (ratio = 1.0)
Monosomy (CN=1):  depth = 15.0    (ratio = 0.5)
Trisomy (CN=3):   depth = 45.0    (ratio = 1.5)
Tetrasomy (CN=4): depth = 60.0    (ratio = 2.0)
Nullisomy (CN=0): depth = 0.0     (ratio = 0.0)
```

### 3.2 Chromosome Names

- Use standard human chromosome naming: chr1, chr21, chrX, etc.
- Down syndrome example: chr21 with trisomy
- Turner syndrome example: chrX with monosomy

---

## 4. Edge Cases

| Case | Input | Expected Behavior | Source |
|------|-------|-------------------|--------|
| Empty input | No depth data | Empty result | Implementation |
| Zero median depth | medianDepth=0 | Empty result (avoid ÷0) | Implementation |
| Negative median depth | medianDepth<0 | Empty result | Implementation |
| Very high depth | ratio >> 5 | CN clamped to 10 | Implementation |
| Very low depth | ratio near 0 | CN = 0 | Implementation |
| Mixed CN per chromosome | Variable depth | Uses dominant CN | Implementation |
| Single data point | 1 bin | Works correctly | Implementation |

---

## 5. Invariants

### 5.1 DetectAneuploidy Output Invariants

```csharp
Assert.Multiple(() => {
    Assert.That(result.CopyNumber, Is.InRange(0, 10));
    Assert.That(result.Confidence, Is.InRange(0.0, 1.0));
    Assert.That(result.Start, Is.LessThan(result.End));
});
```

### 5.2 IdentifyWholeChromosomeAneuploidy Output Invariants

```csharp
Assert.Multiple(() => {
    Assert.That(result.CopyNumber, Is.Not.EqualTo(2));  // Only aneuploid chromosomes
    Assert.That(result.Type, Is.Not.Null.And.Not.Empty);
    Assert.That(result.Chromosome, Is.Not.Null.And.Not.Empty);
});
```

---

## 6. Audit of Existing Tests

### 6.1 Location
- **File:** `ChromosomeAnalyzerTests.cs`
- **Region:** `#region Aneuploidy Detection Tests`

### 6.2 Current Tests

| Test | Status | Notes |
|------|--------|-------|
| DetectAneuploidy_WithNormalDepth_ReturnsCopyNumber2 | ✓ Present | Rename for consistency |
| DetectAneuploidy_WithTrisomy_ReturnsCopyNumber3 | ✓ Present | Rename for consistency |
| DetectAneuploidy_WithMonosomy_ReturnsCopyNumber1 | ✓ Present | Rename for consistency |
| DetectAneuploidy_EmptyInput_ReturnsEmpty | ✓ Present | Keep |
| IdentifyWholeChromosomeAneuploidy_WithTrisomy_IdentifiesTrisomy | ✓ Present | Rename for consistency |
| IdentifyWholeChromosomeAneuploidy_WithNormalChromosome_ReturnsEmpty | ✓ Present | Rename for consistency |

### 6.3 Consolidation Plan

1. **Extract** aneuploidy tests from `ChromosomeAnalyzerTests.cs` to dedicated `ChromosomeAnalyzer_Aneuploidy_Tests.cs`
2. **Add** missing MUST tests (M5, M6, M7, M8, M9, M10, M12, M14, M15, M16)
3. **Add** SHOULD tests for invariants
4. **Rename** existing tests to follow `Method_Scenario_ExpectedResult` pattern without `With` prefix
5. **Keep** note in original file pointing to dedicated test file

---

## 7. Decisions

| Decision | Rationale |
|----------|-----------|
| Extract to dedicated file | Consistent with other CHROM-* test units |
| 80% default minFraction | Implementation default; test with explicit values |
| No sex chromosome special case testing | Implementation does not differentiate; document as limitation |

---

## 8. Open Questions

None remaining after evidence analysis.
