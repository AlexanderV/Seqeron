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
| M5 | CalculateCodonUsage_AllSixtyFourCodons_HandlesAll | All 64 standard codons counted (each ×1) | Genetic code completeness |
| M6 | CompareCodonUsage_IdenticalSequences_ReturnsOne | Same sequence → similarity = 1.0 | Mathematical identity |
| M7 | CompareCodonUsage_DifferentCodons_ReturnsLessThanOne | Partial overlap → 0.5 (exact TVD value) | TVD formula derivation |
| M8 | CompareCodonUsage_EmptySequences_ReturnsZero | Empty inputs → 0.0 | Edge case handling |
| M9 | CompareCodonUsage_Symmetric_SameResult | Sim(a,b) = Sim(b,a) = 0.75 for non-trivial overlap | TVD symmetry + exact derivation |
| M10 | CompareCodonUsage_ResultRange_ZeroToOne | Exact TVD values: 0.0, 0.25, 0.75 | TVD formula derivation |

### Should Tests (Important)

| ID | Test Name | Description | Evidence |
|----|-----------|-------------|----------|
| S1 | CalculateCodonUsage_DnaInput_ConvertsTToU | DNA T→U: full dict {AUG, GCU, UAA}, no T keys | Biological equivalence |
| S2 | CalculateCodonUsage_MixedCase_HandlesCorrectly | Case-insensitive processing | Robustness |
| S3 | CalculateCodonUsage_Invariant_SumEqualsTotalCodons | Sum of counts = n/3 | Mathematical invariant |
| S4 | CompareCodonUsage_OneEmptySequence_ReturnsZero | One empty, one non-empty → 0 | Edge case |
| S5 | CompareCodonUsage_NoOverlappingCodons_ZeroSimilarity | Disjoint codon distributions → 0.0 | TVD formula: disjoint → Σ = 2 → sim = 0 |
| S6 | CompareCodonUsage_PartialOverlap_IntermediateSimilarity | 2/3 codons shared → similarity = 2/3 | TVD formula derivation |

### Could Tests (Nice to have)

| ID | Test Name | Description | Evidence |
|----|-----------|-------------|----------|
| C1 | CalculateCodonUsage_LongSequence_PerformsWell | Performance with large sequences | Non-functional |

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
# M6: Identity
seq1 = seq2 = "AUGGCUGCACUG"
Expected similarity: 1.0

# M7: Partial overlap (exact TVD derivation)
seq1 = "CUGCUGCUGCUA"  → f(CUG)=3/4, f(CUA)=1/4
seq2 = "CUACUACUACUG"  → f(CUA)=3/4, f(CUG)=1/4
Σ|f₁-f₂| = |3/4-1/4| + |1/4-3/4| = 1
Expected similarity: 1 - 1/2 = 0.5

# M9: Symmetry with exact value
seq1 = "AUGAUGCCCUUU"  → f(AUG)=1/2, f(CCC)=1/4, f(UUU)=1/4
seq2 = "AUGUUUUUUCCC"  → f(AUG)=1/4, f(UUU)=1/2, f(CCC)=1/4
Σ|f₁-f₂| = 1/4 + 0 + 1/4 = 1/2
Expected: sim(a,b) = sim(b,a) = 3/4

# M10 case 1: Disjoint
AUG vs CCC → sim = 0.0

# M10 case 2: High difference
seq1 = "AUGAUGAUGAUG"  (AUG×4)  
seq2 = "AUGCCCCCCCCC"  (AUG×1, CCC×3)
Σ = |1-1/4| + |0-3/4| = 3/4+3/4 = 3/2
Expected similarity: 1 - 3/4 = 0.25

# M10 case 3: Low difference
seq1 = "AUGAUGAUGCCC"  (AUG×3, CCC×1)  
seq2 = "AUGAUGCCCCCC"  (AUG×2, CCC×2)
Σ = |3/4-1/2| + |1/4-1/2| = 1/4+1/4 = 1/2
Expected similarity: 1 - 1/4 = 0.75

# S5: Disjoint codons → zero
seq1 = "UUUUUUUUU" (all UUU)
seq2 = "GGGGGGGGG" (all GGG)
Σ|f₁-f₂| = |1-0| + |0-1| = 2
Expected similarity: 1 - 2/2 = 0.0

# S6: Partial overlap (exact TVD derivation)
seq1 = "AUGGCUAUG"  → f(AUG)=2/3, f(GCU)=1/3
seq2 = "AUGUUUAUG"  → f(AUG)=2/3, f(UUU)=1/3
Σ|f₁-f₂| = 0 + 1/3 + 1/3 = 2/3
Expected similarity: 1 - (2/3)/2 = 2/3

# M8: Empty
seq1 = "", seq2 = ""
Expected similarity: 0.0
```

## Invariants to Verify

1. **Count sum invariant**: `sum(counts.Values) == sequence.Length / 3`
2. **Identity**: `CompareCodonUsage(s, s) == 1.0` for non-empty s
3. **Symmetry**: `CompareCodonUsage(a, b) == CompareCodonUsage(b, a)`
4. **Range**: `0 <= CompareCodonUsage(a, b) <= 1`

## Deviations and Assumptions

**None.** All tests and implementation are grounded in external authoritative sources.

- Codon usage tables (E. coli K12, S. cerevisiae, H. sapiens) verified against Kazusa Codon Usage Database (March 2026).
- Comparison metric: Total Variation Distance similarity `1 - Σ|f₁(c)-f₂(c)|/2` — standard metric from probability theory for comparing discrete distributions.
- All expected test values derived analytically from the TVD formula; zero internal assumptions.

## Open Questions

None — algorithm behavior is well-documented.

## Decisions

1. **D1**: Tests focus on `CodonOptimizer` methods, not `CodonUsageAnalyzer` (separate test unit)
2. **D2**: Similarity uses Total Variation Distance (TVD): `1 - Σ|f₁-f₂|/2`, range [0, 1]

## Test File Location
`tests/Seqeron/Seqeron.Genomics.Tests/CodonOptimizer_CodonUsage_Tests.cs`

## Coverage Classification

### Canonical (`CodonOptimizer_CodonUsage_Tests.cs`) — 22 test instances

| # | Test Method | Spec ID | Status |
|---|-------------|---------|--------|
| 1 | `CalculateCodonUsage_SimpleCodingSequence_CountsCorrectly` | M1 | ✅ |
| 2 | `CalculateCodonUsage_EmptySequence_ReturnsEmptyDictionary` | M2 | ✅ |
| 3 | `CalculateCodonUsage_IncompleteCodon_IgnoresTrailing` | M3 | ✅ |
| 4 | `CalculateCodonUsage_RepeatedCodons_CountsAccurately` | M4 | ✅ |
| 5 | `CalculateCodonUsage_AllSixtyFourCodons_HandlesAll` | M5 | ✅ |
| 6 | `CompareCodonUsage_IdenticalSequences_ReturnsOne` | M6 | ✅ |
| 7 | `CompareCodonUsage_DifferentCodons_ReturnsLessThanOne` | M7 | ✅ |
| 8 | `CompareCodonUsage_EmptySequences_ReturnsZero` | M8 | ✅ |
| 9 | `CompareCodonUsage_Symmetric_SameResultBothDirections` | M9 | ✅ |
| 10 | `CompareCodonUsage_ResultRange_ZeroToOne` (×3 cases) | M10 | ✅ |
| 11 | `CalculateCodonUsage_DnaInput_ConvertsTToU` | S1 | ✅ |
| 12 | `CalculateCodonUsage_MixedCase_HandlesCorrectly` | S2 | ✅ |
| 13 | `CalculateCodonUsage_Invariant_SumEqualsTotalCodons` | S3 | ✅ |
| 14 | `CompareCodonUsage_OneEmptySequence_ReturnsZero` | S4 | ✅ |
| 15 | `CompareCodonUsage_NoOverlappingCodons_ZeroSimilarity` | S5 | ✅ |
| 16 | `CompareCodonUsage_PartialOverlap_IntermediateSimilarity` | S6 | ✅ |
| 17 | `CalculateCodonUsage_NullSequence_ReturnsEmpty` | edge | ✅ |
| 18 | `CalculateCodonUsage_TooShort_ReturnsEmpty` (×2 cases) | edge | ✅ |
| 19 | `CompareCodonUsage_SingleCodon_IdenticalReturnsOne` | edge | ✅ |

### Coverage Summary

| Category | Total | ✅ Covered | ⚠ Weak | ❌ Missing | 🔁 Duplicate |
|----------|-------|-----------|--------|-----------|-------------|
| Must | 10 | 10 | 0 | 0 | 0 |
| Should | 6 | 6 | 0 | 0 | 0 |
| Edge | 3 | 3 | 0 | 0 | 0 |
| **Total** | **19** | **19** | **0** | **0** | **0** |

### TVD Exact Values Verified

All `CompareCodonUsage` tests use exact values derived analytically from the TVD formula:

| Similarity | Test(s) | Derivation |
|------------|---------|------------|
| 0.0 | M8, M10①, S5 | Disjoint: Σ=2 → 1−1=0 |
| 0.25 | M10② | f(AUG)=1 vs f(AUG)=1/4,f(CCC)=3/4: Σ=3/2 → 1−3/4 |
| 0.5 | M7 | Mirrored CUG/CUA: Σ=1 → 1−1/2 |
| 2/3 | S6 | Shared AUG, disjoint GCU/UUU: Σ=2/3 → 1−1/3 |
| 0.75 | M9, M10③ | Various partial overlaps: Σ=1/2 → 1−1/4 |
| 1.0 | M6, edge | Identity: Σ=0 → 1−0=1 |

## Last Updated
2026-03-11
