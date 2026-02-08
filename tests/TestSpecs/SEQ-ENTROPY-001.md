# SEQ-ENTROPY-001: Shannon Entropy Test Specification

**Test Unit ID:** SEQ-ENTROPY-001
**Area:** Composition
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-01-22

---

## Scope

### Canonical Methods
| Method | Class | Type |
|--------|-------|------|
| `CalculateShannonEntropy(DnaSequence)` | SequenceComplexity | Canonical |
| `CalculateShannonEntropy(string)` | SequenceComplexity | Overload |
| `CalculateKmerEntropy(DnaSequence, k)` | SequenceComplexity | K-mer variant |

### Related Methods (Delegate/Wrapper - smoke only)
| Method | Class | Notes |
|--------|-------|-------|
| `CalculateShannonEntropy(string)` | SequenceStatistics | General-purpose, different behavior |

---

## Evidence

### Sources
1. **Wikipedia - Entropy (information theory)**
   - Formula: $H(X) = -\sum p(x) \log_2 p(x)$
   - For n symbols with uniform distribution: $H_{max} = \log_2(n)$
   - For DNA (n=4): $H_{max} = 2$ bits
   - For empty distribution or single symbol: H = 0

2. **Wikipedia - Sequence logo**
   - Information content at position: $R_i = \log_2(4) - H_i$
   - Confirms DNA max entropy = 2 bits

3. **Wikipedia - K-mer**
   - K-mers are overlapping substrings of length k
   - Used in bioinformatics for alignment-free sequence analysis

4. **Shannon (1948) - A Mathematical Theory of Communication**
   - Original definition of information entropy
   - Convention: 0 × log(0) = 0 (limit)

### Test Datasets from Sources
- Uniform distribution (ATGC) → H = 2.0 bits
- Single base (AAAA) → H = 0.0 bits
- Two bases 50/50 (ATAT) → H = 1.0 bit
- Three bases 1/3 each → H = log₂(3) ≈ 1.585 bits

---

## Invariants

### INV-ENT-001: Entropy Range for DNA
For any DNA sequence: $0 \leq H \leq 2$ bits

### INV-ENT-002: Maximum Entropy
Uniform distribution of 4 bases yields exactly 2.0 bits

### INV-ENT-003: Zero Entropy for Homopolymers
Single-base sequences (homopolymers) yield exactly 0.0 bits

### INV-ENT-004: K-mer Entropy Range
For any DNA sequence with k-mers: $0 \leq H_k \leq 2k$ bits
(Max when all 4^k possible k-mers appear equally)

### INV-ENT-005: Empty Sequence
Empty sequence → H = 0.0 (no information)

---

## Test Cases

### MUST Tests (Evidence-Based)

| ID | Test Case | Input | Expected | Source |
|----|-----------|-------|----------|--------|
| M01 | Uniform 4 bases | "ATGC" | 2.0 bits | Wikipedia - max entropy |
| M02 | Homopolymer | "AAAA" (or any single base) | 0.0 bits | Wikipedia - min entropy |
| M03 | Two bases 50/50 | "ATAT" or "ATATATAT" | 1.0 bit | Binary entropy formula |
| M04 | Empty sequence | "" | 0.0 | Convention: no data = no entropy |
| M05 | Range invariant | Any DNA sequence | 0 ≤ H ≤ 2 | Mathematical bound |
| M06 | K-mer empty result | seq.Length < k | 0.0 | No k-mers extractable |
| M07 | K-mer single pattern | "AAAAAAAAAA", k=2 | 0.0 | Only "AA" k-mer |
| M08 | DnaSequence overload | DnaSequence("ATGC") | 2.0 | API consistency |
| M09 | Null DnaSequence | null | ArgumentNullException | Defensive programming |
| M10 | Invalid k | k < 1 | ArgumentOutOfRangeException | Parameter validation |

### SHOULD Tests (Good Practice)

| ID | Test Case | Input | Expected | Rationale |
|----|-----------|-------|----------|-----------|
| S01 | Case insensitivity | "atgc" vs "ATGC" | Equal results | Robustness |
| S02 | String overload matches DnaSequence | Same sequence | Equal results | API consistency |
| S03 | Three bases equal | "ATGATGATG" (A,T,G=33% each) | ~1.585 bits | log₂(3) |
| S04 | K-mer varied | "ATGCATGCATGC", k=2 | > 1.5 bits | High dinucleotide diversity |
| S05 | K-mer moderate variety | "ATGCTAGCATGC", k=2 | > 0, reasonable | Sanity check |

### COULD Tests (Extended Coverage)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| C01 | Very long sequence | Performance consistency |
| C02 | K-mer k=1 equivalence | Should approximate base entropy |
| C03 | Non-DNA characters ignored | Verify only ATGC counted |

---

## Audit Results

### Existing Test Coverage (Pre-Consolidation)

**SequenceComplexityTests.cs** (Canonical Location):
- `CalculateShannonEntropy_EqualBases_ReturnsTwo` → M01 ✓
- `CalculateShannonEntropy_SingleBase_ReturnsZero` → M02 ✓
- `CalculateShannonEntropy_TwoBases_ReturnsOne` → M03 ✓
- `CalculateShannonEntropy_EmptySequence_ReturnsZero` → M04 ✓
- `CalculateShannonEntropy_StringOverload_Works` → S02 partial ✓
- `CalculateShannonEntropy_NullSequence_ThrowsException` → M09 ✓
- `CalculateKmerEntropy_VariedDinucleotides_ReturnsHigh` → S04 ✓
- `CalculateKmerEntropy_RepeatedDinucleotides_ReturnsLow` → M07 ✓
- `CalculateKmerEntropy_SequenceShorterThanK_ReturnsZero` → M06 ✓

**SequenceStatisticsTests.cs** (Delegate/Wrapper):
- `CalculateShannonEntropy_UniformDistribution_ReturnsHighEntropy` → Smoke
- `CalculateShannonEntropy_HomopolymerRun_ReturnsZero` → Smoke
- `CalculateShannonEntropy_EmptyString_ReturnsZero` → Smoke
- Tests null handling in general null-safety test → Smoke

### Gap Analysis
| ID | Coverage Status | Action |
|----|-----------------|--------|
| M05 | ✅ Covered | Range invariant test added |
| M08 | ✅ Covered | Strengthened to verify exact equality |
| M10 | ✅ Covered | k parameter validation test added |
| S01 | ✅ Covered | Case insensitivity test added |
| S03 | ✅ Covered | Three-base distribution test added |

### Consolidation Plan
1. **Canonical file:** `SequenceComplexityTests.cs` - deep evidence-based tests
2. **Wrapper smoke:** `SequenceStatisticsTests.cs` - keep existing 3 tests as delegation verification
3. **No duplicates** - different implementations (SequenceComplexity counts only ATGC, SequenceStatistics counts all letters)

---

## Open Questions / Decisions

1. **Q:** Should non-ATGC characters affect the entropy calculation?
   **A:** Current implementation ignores them (counts only ATGC). This is documented behavior. ✓

2. **Q:** Is SequenceStatistics.CalculateShannonEntropy a delegate?
   **A:** No - it's a separate general-purpose implementation that counts all letters. Tests remain separate.

---

## Validation Checklist

- [x] All MUST tests have evidence source
- [x] Invariants are mathematically verifiable
- [x] Edge cases documented
- [x] API consistency verified
- [x] Parameter validation tested
