# Test Specification: CHROM-CENT-001

## Test Unit Information
- **ID:** CHROM-CENT-001
- **Area:** Chromosome
- **Canonical Method:** `ChromosomeAnalyzer.AnalyzeCentromere(chromosomeName, sequence, windowSize, minAlphaSatelliteContent)`
- **Complexity:** O(n)

## Methods Under Test

| Method | Class | Type | Deep Test |
|--------|-------|------|-----------|
| `AnalyzeCentromere(chrName, seq, windowSize)` | ChromosomeAnalyzer | Canonical | Yes |

## Test Classification

### Must Tests (Evidence-backed)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| M1 | Empty sequence returns Unknown type with null boundaries | Standard edge case | Wikipedia |
| M2 | Sequence shorter than window size returns Unknown | Cannot perform analysis without sufficient data | Implementation logic |
| M3 | Chromosome name is preserved in result | Basic correctness | API contract |
| M4 | Repetitive region is detected as centromeric | Centromeres are characterized by repetitive DNA | Wikipedia |
| M5 | CentromereResult invariants: Start <= End when found | Structural invariant | Logic |
| M6 | Length equals (End - Start) when centromere found | Structural invariant | Logic |
| M7 | Type classification based on position matches biology | Metacentric = middle, Acrocentric = near end | Levan (1964) |
| M8 | IsAcrocentric flag matches type | Consistency check | API contract |

### Should Tests (Recommended)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| S1 | Different window sizes affect detection sensitivity | Parameter behavior |
| S2 | Different alpha-satellite thresholds affect detection | Parameter behavior |
| S3 | Case insensitivity in sequence | Robustness |

### Could Tests (Optional)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| C1 | Performance with large sequences | Non-functional |

## Edge Cases

| Case | Input | Expected Output |
|------|-------|-----------------|
| Empty sequence | `""` | Type = "Unknown", Start = null, End = null |
| Null sequence | `null` | Type = "Unknown", Start = null, End = null |
| Very short sequence | `"ATCG"` | Type = "Unknown" |
| All same base | `"AAAA..."` (long) | May detect repetitive pattern |
| No repetitive regions | Random non-repetitive | Type = "Unknown" or low score |

## Test Invariants

1. **Result chromosome name equals input name**
2. **When Start/End are both null, Length = 0**
3. **When Start/End are both non-null, Length = End - Start**
4. **Type is one of: Metacentric, Submetacentric, Acrocentric, Unknown**
5. **IsAcrocentric = true iff Type = "Acrocentric"**
6. **AlphaSatelliteContent >= 0**

## Audit of Existing Tests

### Existing Test Location
- File: `ChromosomeAnalyzerTests.cs`
- Region: `#region Centromere Analysis Tests`

### Coverage Assessment

| Test | Status | Notes |
|------|--------|-------|
| `AnalyzeCentromere_WithRepetitiveRegion_FindsCentromere` | Weak | Only checks chromosome name, not actual detection |
| `AnalyzeCentromere_EmptySequence_ReturnsUnknown` | Covered | Correct |
| `AnalyzeCentromere_ShortSequence_HandlesGracefully` | Covered | Correct |

### Consolidation Plan

1. Move centromere tests to dedicated file: `ChromosomeAnalyzer_Centromere_Tests.cs`
2. Remove weak test `AnalyzeCentromere_WithRepetitiveRegion_FindsCentromere` (only checks name, not detection)
3. Add comprehensive tests for all Must cases
4. Add invariant-based assertions using Assert.Multiple

## Open Questions / Decisions

1. **Decision**: The implementation uses simplified heuristics for centromere detection. Real centromere detection requires specialized databases (alpha-satellite consensus). Tests will verify the heuristic logic works correctly.

2. **Decision**: Position-based classification thresholds are based on the implementation's `DetermineCentromereType` method, which aligns with biological nomenclature but uses simplified position ratios.
