# TestSpec: CODON-USAGE-001 - Codon Usage Analysis

## Test Unit ID
CODON-USAGE-001

## Overview
Tests for codon usage calculation and comparison methods in `CodonOptimizer`.

## Methods Under Test

| Method | Class | Type | Complexity |
|--------|-------|------|------------|
| `CalculateCodonUsage(string)` | CodonOptimizer | Canonical | O(n) |
| `CompareCodonUsage(string, string)` | CodonOptimizer | Canonical | O(n+m) |

## Test Categories

### Must Tests (Required)

| ID | Test Name | Description | Evidence |
|----|-----------|-------------|----------|
| M1 | CalculateCodonUsage_SimpleCodingSequence_CountsCorrectly | Verify basic codon counting with known input | Wikipedia, Kazusa |
| M2 | CalculateCodonUsage_EmptySequence_ReturnsEmptyDictionary | Empty input returns empty result | Standard edge case |
| M3 | CalculateCodonUsage_IncompleteCodon_IgnoresTrailing | Trailing 1-2 nucleotides ignored | Kazusa format |
| M4 | CalculateCodonUsage_RepeatedCodons_CountsAccurately | Repeated codons counted correctly | Mathematical invariant |
| M5 | CalculateCodonUsage_AllSixtyFourCodons_HandlesAll | All valid codons recognized | Genetic code completeness |
| M6 | CompareCodonUsage_IdenticalSequences_ReturnsOne | Same sequence → similarity = 1.0 | Mathematical identity |
| M7 | CompareCodonUsage_DifferentCodons_ReturnsLessThanOne | Different codon usage → < 1.0 | Distance metric property |
| M8 | CompareCodonUsage_EmptySequences_ReturnsZero | Empty inputs → 0.0 | Edge case handling |
| M9 | CompareCodonUsage_Symmetric_SameResult | Sim(a,b) = Sim(b,a) | Metric symmetry property |
| M10 | CompareCodonUsage_ResultRange_ZeroToOne | Result always in [0, 1] | Mathematical range |

### Should Tests (Important)

| ID | Test Name | Description | Evidence |
|----|-----------|-------------|----------|
| S1 | CalculateCodonUsage_DnaInput_ConvertsTToU | DNA T converted to RNA U | Implementation behavior |
| S2 | CalculateCodonUsage_MixedCase_HandlesCorrectly | Case-insensitive processing | Robustness |
| S3 | CalculateCodonUsage_Invariant_SumEqualsTotalCodons | Sum of counts = n/3 | Mathematical invariant |
| S4 | CompareCodonUsage_OneEmptySequence_ReturnsZero | One empty, one non-empty → 0 | Edge case |
| S5 | CompareCodonUsage_NoOverlappingCodons_LowSimilarity | Different codon sets → low similarity | Manhattan distance |

### Could Tests (Nice to have)

| ID | Test Name | Description | Evidence |
|----|-----------|-------------|----------|
| C1 | CalculateCodonUsage_LongSequence_PerformsWell | Performance with large sequences | Non-functional |
| C2 | CompareCodonUsage_PartialOverlap_IntermediateSimilarity | Partial overlap → intermediate result | ASSUMPTION |

## Test Data

### Simple Test Cases

```
Sequence: "AUGGCUGCU" (M-A-A)
Expected: {"AUG": 1, "GCU": 2}

Sequence: "" (empty)
Expected: {}

Sequence: "AUGGCUG" (AUG + incomplete)
Expected: {"AUG": 1, "GCU": 1}
```

### Comparison Test Cases

```
seq1 = seq2 = "AUGGCUGCACUG"
Expected similarity: 1.0

seq1 = "CUGCUGCUGCUG" (all CUG)
seq2 = "CUACUACUACUA" (all CUA)
Expected similarity: < 1.0 (no shared codons)

seq1 = "", seq2 = ""
Expected similarity: 0.0
```

## Invariants to Verify

1. **Count sum invariant**: `sum(counts.Values) == sequence.Length / 3`
2. **Identity**: `CompareCodonUsage(s, s) == 1.0` for non-empty s
3. **Symmetry**: `CompareCodonUsage(a, b) == CompareCodonUsage(b, a)`
4. **Range**: `0 <= CompareCodonUsage(a, b) <= 1`

## Audit of Existing Tests

### Current Tests in CodonOptimizerTests.cs

| Test | Status | Action |
|------|--------|--------|
| `CalculateCodonUsage_CountsCodons` | Weak | Move to dedicated file, enhance |
| `CalculateCodonUsage_EmptySequence_ReturnsEmptyDictionary` | Adequate | Move to dedicated file |
| `CompareCodonUsage_IdenticalSequences_HighSimilarity` | Adequate | Move to dedicated file |
| `CompareCodonUsage_DifferentSequences_LowerSimilarity` | Adequate | Move to dedicated file |
| `CompareCodonUsage_EmptySequences_ReturnsZero` | Adequate | Move to dedicated file |

### Consolidation Plan

1. ~~Create new file: `CodonOptimizer_CodonUsage_Tests.cs`~~ ✅ Done
2. ~~Move existing 5 tests from `CodonOptimizerTests.cs`~~ ✅ Done
3. ~~Add missing Must/Should tests~~ ✅ Done
4. ~~Remove duplicates from original file~~ ✅ Done
5. ~~Keep other test categories in `CodonOptimizerTests.cs`~~ ✅ Done

## Open Questions

None - algorithm behavior is well-documented.

## Decisions

1. **D1**: Tests focus on `CodonOptimizer` methods, not `CodonUsageAnalyzer` (separate test unit)
2. **D2**: Similarity range [0, 1] is implementation-specific (Manhattan-based metric)

## Test File Location
`tests/Seqeron/Seqeron.Genomics.Tests/CodonOptimizer_CodonUsage_Tests.cs`

## Last Updated
2026-02-04
