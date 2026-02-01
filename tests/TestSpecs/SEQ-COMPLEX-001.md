# TestSpec: SEQ-COMPLEX-001 - Linguistic Complexity

**Test Unit ID:** SEQ-COMPLEX-001
**Area:** Sequence Composition
**Created:** 2026-01-22
**Status:** Complete

## Scope

This TestSpec covers the Linguistic Complexity (LC) calculation for DNA sequences, including:
- `SequenceComplexity.CalculateLinguisticComplexity(DnaSequence, int)` — Canonical
- `SequenceComplexity.CalculateLinguisticComplexity(string, int)` — String overload
- `SequenceComplexity.FindLowComplexityRegions(...)` — Region detection
- `SequenceComplexity.MaskLowComplexity(...)` — Low-complexity masking

## Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| Wikipedia "Linguistic sequence complexity" | Encyclopedia | Formula definitions, range [0,1] |
| Troyanskaya et al. (2002) Bioinformatics | Peer-reviewed | Summation formula: LC = Σ(observed) / Σ(possible) |
| Orlov & Potapov (2004) NAR | Peer-reviewed | V_max = min(4^i, N-i+1), complexity profiles |
| Trifonov (1990) | Original paper | Concept introduction, vocabulary richness |

## Test Categories

### MUST Tests (Evidence-backed, Required)

| ID | Test Case | Invariant/Property | Source |
|----|-----------|-------------------|--------|
| M1 | Empty sequence returns 0 | No vocabulary possible | Implementation definition |
| M2 | Result in range [0, 1] | Mathematical property | Troyanskaya et al. (2002) |
| M3 | Homopolymer has low complexity | Single N-mer type per length | Orlov & Potapov (2004) |
| M4 | Random-like sequence has high complexity | Maximum vocabulary diversity | Orlov & Potapov (2004) |
| M5 | Null DnaSequence throws ArgumentNullException | Guard clause | Implementation contract |
| M6 | Invalid maxWordLength throws ArgumentOutOfRangeException | Guard clause | Implementation contract |
| M7 | String overload produces same result as DnaSequence | API consistency | Implementation contract |
| M8 | Single nucleotide returns positive value | Has vocabulary | Definition |

### SHOULD Tests (Recommended)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| S1 | Repetitive dinucleotide pattern has lower complexity than random | Reduced vocabulary |
| S2 | Longer sequences maintain valid range | Scalability |
| S3 | Case insensitivity (lowercase DNA valid) | Robustness |
| S4 | maxWordLength parameter affects result | Parameter functionality |

### COULD Tests (Nice-to-have)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| C1 | Performance benchmark for large sequences | Verify O(n×k) complexity |
| C2 | Comparison with known reference values | Implementation verification |

## FindLowComplexityRegions Tests

### MUST Tests

| ID | Test Case | Invariant/Property | Source |
|----|-----------|-------------------|--------|
| LCR-M1 | Null sequence throws ArgumentNullException | Guard clause | Contract |
| LCR-M2 | Invalid windowSize throws ArgumentOutOfRangeException | Guard clause | Contract |
| LCR-M3 | High-complexity uniform sequence returns empty | No regions below threshold | Definition |
| LCR-M4 | Sequence with poly-A tract returns region(s) | Detects low-complexity | Orlov & Potapov (2004) |
| LCR-M5 | Returned regions contain correct sequence substring | Data integrity | Contract |

### SHOULD Tests

| ID | Test Case | Rationale |
|----|-----------|-----------|
| LCR-S1 | Threshold parameter affects detection sensitivity | Parameter functionality |
| LCR-S2 | Window size affects region boundaries | Parameter functionality |

## MaskLowComplexity Tests

### MUST Tests

| ID | Test Case | Invariant/Property | Source |
|----|-----------|-------------------|--------|
| MASK-M1 | Null sequence throws ArgumentNullException | Guard clause | Contract |
| MASK-M2 | Low-complexity regions are masked with specified character | Core functionality | Definition |
| MASK-M3 | High-complexity regions are preserved | Core functionality | Definition |
| MASK-M4 | Custom mask character is applied | Parameter functionality | Contract |
| MASK-M5 | Result length equals input length | Invariant | Definition |

## Audit of Existing Tests

### Current Test File
`Seqeron.Genomics.Tests/SequenceComplexityTests.cs`

### Coverage Analysis

| Test | Maps To | Status | Notes |
|------|---------|--------|-------|
| `CalculateLinguisticComplexity_HighComplexity_ReturnsHigh` | M4 | ✓ Covered | Good |
| `CalculateLinguisticComplexity_LowComplexity_ReturnsLow` | M3 | ✓ Covered | Good |
| `CalculateLinguisticComplexity_EmptySequence_ReturnsZero` | M1 | ✓ Covered | Good |
| `CalculateLinguisticComplexity_RangeIsZeroToOne` | M2 | Weak | Only checks one sequence |
| `CalculateLinguisticComplexity_StringOverload_Works` | M7 | Weak | Only checks > 0 |
| `CalculateLinguisticComplexity_NullSequence_ThrowsException` | M5 | ✓ Covered | Good |
| `CalculateLinguisticComplexity_ZeroWordLength_ThrowsException` | M6 | ✓ Covered | Good |

### Missing Coverage

| ID | Gap | Priority |
|----|-----|----------|
| GAP-1 | M8: Single nucleotide test | MUST |
| GAP-2 | M2: Range invariant with multiple sequences (Assert.Multiple) | MUST |
| GAP-3 | M7: String overload equality test | MUST |
| GAP-4 | S1: Dinucleotide repeat vs random comparison | SHOULD |
| GAP-5 | S4: maxWordLength parameter effect | SHOULD |

### Weak Tests to Strengthen

| Test | Issue | Fix |
|------|-------|-----|
| `CalculateLinguisticComplexity_RangeIsZeroToOne` | Single sequence | Add multiple sequences with Assert.Multiple |
| `CalculateLinguisticComplexity_StringOverload_Works` | Only checks > 0 | Verify equality with DnaSequence version |

## Consolidation Plan

1. **Canonical file:** `SequenceComplexityTests.cs` (existing)
2. **No wrapper/delegate separation needed** — all methods are in SequenceComplexity class
3. **Actions:**
   - Strengthen existing weak tests
   - Add missing MUST tests
   - Add relevant SHOULD tests
   - Organize tests into clear regions

## Open Questions / Decisions

None — algorithm is well-defined by authoritative sources.

## ASSUMPTIONS

None — all test rationale is backed by evidence from sources listed above.
