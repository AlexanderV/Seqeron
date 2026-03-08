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
| M17 | IdentifyWholeChromosomeAneuploidy_Pentasomy_IdentifiesPentasomy | Classification correctness | Wikipedia |

### 2.2 SHOULD Tests (Recommended)

| ID | Test Name | Rationale |
|----|-----------|-----------|
| S1 | DetectAneuploidy_Confidence_InRangeZeroToOne | Output invariant |
| S2 | DetectAneuploidy_LogRatio_MatchesFormula | Mathematical correctness |
| S3 | DetectAneuploidy_OutputOrdered_ByPosition | Ordering invariant |
| S4 | IdentifyWholeChromosomeAneuploidy_HighCopyNumber_FormatsCorrectly | CN > 5 formatting |

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

## 6. Coverage Classification

### 6.1 Location
- **File:** `ChromosomeAnalyzer_Aneuploidy_Tests.cs`
- **Note in:** `ChromosomeAnalyzerTests.cs` (redirect to dedicated file)

### 6.2 DetectAneuploidy — 16 test methods

| # | Test Method | Spec ID | Status |
|---|-------------|---------|--------|
| 1 | `DetectAneuploidy_NormalDepth_ReturnsCopyNumber2` | M1 | ✅ |
| 2 | `DetectAneuploidy_Trisomy_ReturnsCopyNumber3` | M2 | ✅ |
| 3 | `DetectAneuploidy_Monosomy_ReturnsCopyNumber1` | M3 | ✅ |
| 4 | `DetectAneuploidy_EmptyInput_ReturnsEmpty` | M4 | ✅ |
| 5 | `DetectAneuploidy_ZeroMedianDepth_ReturnsEmpty` | M5 | ✅ |
| 6 | `DetectAneuploidy_NegativeMedianDepth_ReturnsEmpty` | M5+ | ✅ |
| 7 | `DetectAneuploidy_Nullisomy_ReturnsCopyNumber0` | M6 | ✅ |
| 8 | `DetectAneuploidy_Tetrasomy_ReturnsCopyNumber4` | M7 | ✅ |
| 9 | `DetectAneuploidy_VeryHighDepth_CopyNumberClampedTo10` | M10 | ✅ |
| 10 | `DetectAneuploidy_VeryLowDepth_CopyNumberMinimumIsZero` | M10+ | ✅ |
| 11 | `DetectAneuploidy_MultipleChromosomes_GroupsCorrectly` | M8 | ✅ |
| 12 | `DetectAneuploidy_BinAggregation_AveragesDepth` | M9 | ✅ |
| 13 | `DetectAneuploidy_SingleDataPoint_WorksCorrectly` | Edge | ✅ |
| 14 | `DetectAneuploidy_LargeBinSize_ReducesOutput` | C1 | ✅ |
| 15 | `DetectAneuploidy_Confidence_ExactForBoundaryRatios` | S1 | ✅ |
| 16 | `DetectAneuploidy_OutputPositions_StartLessThanEnd` | Inv | ✅ |
| 17 | `DetectAneuploidy_LogRatio_MatchesFormula` | S2 | ✅ |
| 18 | `DetectAneuploidy_OutputOrdered_ByPosition` | S3 | ✅ |

### 6.3 IdentifyWholeChromosomeAneuploidy — 11 test methods

| # | Test Method | Spec ID | Status |
|---|-------------|---------|--------|
| 19 | `IdentifyWholeChromosomeAneuploidy_Trisomy_IdentifiesTrisomy` | M11 | ✅ |
| 20 | `IdentifyWholeChromosomeAneuploidy_Monosomy_IdentifiesMonosomy` | M12 | ✅ |
| 21 | `IdentifyWholeChromosomeAneuploidy_Nullisomy_IdentifiesNullisomy` | M14 | ✅ |
| 22 | `IdentifyWholeChromosomeAneuploidy_Tetrasomy_IdentifiesTetrasomy` | M15 | ✅ |
| 23 | `IdentifyWholeChromosomeAneuploidy_Normal_ReturnsEmpty` | M13 | ✅ |
| 24 | `IdentifyWholeChromosomeAneuploidy_Pentasomy_IdentifiesPentasomy` | M17 | ✅ |
| 25 | `IdentifyWholeChromosomeAneuploidy_HighCopyNumber_FormatsCorrectly` | S4 | ✅ |
| 26 | `IdentifyWholeChromosomeAneuploidy_MixedCN_BelowThreshold_ReturnsEmpty` | M16a | ✅ |
| 27 | `IdentifyWholeChromosomeAneuploidy_MixedCN_AtThreshold_IdentifiesAneuploidy` | M16b | ✅ |
| 28 | `IdentifyWholeChromosomeAneuploidy_CustomMinFraction_Works` | C2 | ✅ |
| 29 | `IdentifyWholeChromosomeAneuploidy_MultipleChromosomes_IndependentClassification` | Edge | ✅ |
| 30 | `IdentifyWholeChromosomeAneuploidy_EmptyInput_ReturnsEmpty` | Edge | ✅ |
| 31 | `IdentifyWholeChromosomeAneuploidy_SingleBin_ClassifiedCorrectly` | Edge | ✅ |

### 6.4 Classification Summary

- ✅ Covered: 31 tests (1 duplicate removed, 1 new test added — net 31)
- ❌ Missing: 0
- ⚠ Weak: 0 (M1–M3, M6–M7: LogRatio + Confidence assertions added; S1: exact values; VeryLowDepth: exact CN)
- 🔁 Duplicate: 0 (`CopyNumber_AlwaysInValidRange` removed — subsumed by M10 + VeryLowDepth)

---

## 7. Decisions

| Decision | Rationale |
|----------|-----------|
| Extract to dedicated file | Consistent with other CHROM-* test units |
| 80% default minFraction | Implementation default; test with explicit values |
| No sex chromosome special case testing | Implementation does not differentiate; documented as limitation |
| Pentasomy explicitly mapped | Wikipedia defines Pentasomy as 5 copies; CN>5 uses fallback format |

---

## 8. Documented Limitations

| Limitation | Source | Impact |
|------------|--------|--------|
| Sex chromosomes treated same as autosomes | Wikipedia (males: 1 copy X/Y is normal) | Male monosomic X/Y flagged as aneuploidy |
| Partial aneuploidy | Wikipedia (translocations) | Detected at bin level, not classified as partial monosomy/trisomy |

---

## 9. Validation Checklist

- [x] All MUST tests (M1–M17) implemented with evidence source
- [x] All SHOULD tests (S1–S4) implemented
- [x] COULD tests (C1, C2) implemented
- [x] Edge cases documented and tested (empty, zero/negative median, single point, large bin)
- [x] Invariants verified (Start < End, output ordered, confidence = 1.0 at boundaries)
- [x] CN detection formula derived from theory: CN = round(depth/median × 2)
- [x] LogRatio values verified against mathematical log₂ identities
- [x] Confidence formula verified: 1 − min(1, |CN/2 − ratio|)
- [x] Classification names match Wikipedia (Nullisomy, Monosomy, Trisomy, Tetrasomy, Pentasomy)
- [x] Threshold boundary (≥ 0.8) tested at exact boundary and below
- [x] No assumptions — all behaviors sourced from Wikipedia or mathematical definitions
- [x] No duplicates — each test serves a distinct purpose
- [x] Coverage classification complete: 0 missing, 0 weak, 0 duplicate
