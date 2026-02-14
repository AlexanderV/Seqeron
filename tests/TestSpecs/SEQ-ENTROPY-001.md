# SEQ-ENTROPY-001: Shannon Entropy Test Specification

**Test Unit ID:** SEQ-ENTROPY-001
**Area:** Composition
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-02-14

---

## Scope

### Canonical Methods
| Method | Class | Type |
|--------|-------|------|
| `CalculateShannonEntropy(DnaSequence)` | SequenceComplexity | Canonical |
| `CalculateShannonEntropy(string)` | SequenceComplexity | Overload |
| `CalculateKmerEntropy(DnaSequence, k)` | SequenceComplexity | K-mer variant |

### Related Methods (Delegate/Wrapper — smoke only)
| Method | Class | Notes |
|--------|-------|-------|
| `CalculateShannonEntropy(string)` | SequenceStatistics | General-purpose, counts all letters (separate implementation) |

---

## Sources

1. **Shannon (1948) — A Mathematical Theory of Communication**
   - Formula: $H(X) = -\sum p(x) \log_2 p(x)$
   - Convention: $0 \times \log(0) = 0$ (limit as $x \to 0$)

2. **Wikipedia — Entropy (information theory)**
   - For $n$ symbols with uniform distribution: $H_{\max} = \log_2(n)$
   - For DNA ($n = 4$): $H_{\max} = 2$ bits
   - "In the case of $p(x) = 0$, the value of the corresponding summand
     $0 \log_b(0)$ is taken to be 0, which is consistent with the limit
     $\lim_{p \to 0^{+}} p \log(p) = 0$"
   - "The maximal entropy of an event with $n$ different outcomes is $\log_b(n)$:
     it is attained by the uniform probability distribution"
   - "Adding or removing an event with probability zero does not contribute
     to the entropy: $H_{n+1}(p_1, \dots, p_n, 0) = H_n(p_1, \dots, p_n)$"

3. **Wikipedia — Sequence logo**
   - Information content: $R_i = \log_2(4) - H_i$
   - Confirms DNA max entropy = 2 bits

4. **Wikipedia — K-mer**
   - K-mers are overlapping substrings of length $k$

---

## Invariants

### INV-ENT-001: Entropy Range for DNA
For any DNA sequence: $0 \leq H \leq 2$ bits
(Source: $H_{\max} = \log_2(4) = 2$ — Wikipedia)

### INV-ENT-002: Maximum Entropy
Uniform distribution of 4 bases yields exactly $H = 2.0$ bits
(Source: $H_{\max} = \log_2(n)$ — Wikipedia)

### INV-ENT-003: Zero Entropy for Homopolymers
Single-base sequences yield exactly $H = 0.0$ bits
(Source: $p = 1 \Rightarrow -1 \cdot \log_2(1) = 0$ — Shannon 1948)

### INV-ENT-004: K-mer Entropy Range
For any DNA sequence with k-mers: $0 \leq H_k \leq \log_2(4^k) = 2k$ bits
(Source: maximum when all $4^k$ possible k-mers appear equally — Shannon 1948)

### INV-ENT-005: Empty / Degenerate Input
Empty sequence → $H = 0.0$ (empty sum = 0)
No extractable k-mers (length < k) → $H = 0.0$ (empty sum = 0)
(Source: empty sum convention, consistent with Wikipedia expansibility axiom)

---

## Test Cases

### MUST Tests (Evidence-Based)

| ID | Test Case | Input | Expected | Source |
|----|-----------|-------|----------|--------|
| M01 | Uniform 4 bases (DnaSequence) | `DnaSequence("ATGCATGCATGCATGC")` | $2.0$ bits exact | Shannon; $p = 0.25$, $\log_2(0.25) = -2$ |
| M02 | Homopolymer | `DnaSequence("AAAAAAA")` | $0.0$ exact | Shannon; $p = 1$, $\log_2(1) = 0$ |
| M03 | Two bases 50/50 | `DnaSequence("ATATATAT")` | $1.0$ exact | Shannon; $p = 0.5$, $\log_2(0.5) = -1$ |
| M04 | Empty sequence | `""` | $0.0$ | Empty sum = 0 (Wikipedia expansibility) |
| M05 | Range invariant | Any DNA sequence | $0 \leq H \leq 2$ | $H_{\max} = \log_2(4) = 2$ (Wikipedia) |
| M06 | K-mer: length < k | `DnaSequence("AT")`, k=5 | $0.0$ | No k-mers extractable → empty sum |
| M07 | K-mer: homopolymer | `DnaSequence("AAAAAAAAAA")`, k=2 | $0.0$ exact | Only "AA" → $p = 1$, $H = 0$ |
| M08 | Non-ATGC chars ignored | `"ATGCNN"`, `"NNNNNN"` | $2.0$, $0.0$ | Alphabet = {A,T,G,C}; non-DNA excluded from numerator and denominator |
| M09 | Null DnaSequence | `null` | `ArgumentNullException` | Guard clause (.NET convention) |
| M10 | Invalid k | k < 1 | `ArgumentOutOfRangeException` | Guard clause (.NET convention) |

### SHOULD Tests (Good Practice)

| ID | Test Case | Input | Expected | Source |
|----|-----------|-------|----------|--------|
| S01 | Case insensitivity | `"atgc"` vs `"ATGC"` | Bitwise equal | `ToUpperInvariant()` produces same core input |
| S02 | String = DnaSequence overload | Same uppercase seq | Bitwise equal | Same core function, same input |
| S03 | Three bases equal | `DnaSequence("ATGATGATG")` | $\log_2(3) \approx 1.58496$ | Shannon; $n = 3$ uniform |
| S04 | K-mer varied | `DnaSequence("ATGCATGCATGCATGC")`, k=2 | hand-calculated exact | Shannon formula on k-mer distribution |
| S05 | K-mer uniform | `DnaSequence("ATCG")`, k=2 | $\log_2(3) \approx 1.58496$ | 3 unique dinucleotides, uniform |

---

## Hand-Calculated Cross-Verification Table

### Base Shannon Entropy

| Input | Bases | Distribution | Formula | H (bits) | Float-exact |
|-------|-------|-------------|---------|-----------|-------------|
| `"ATGCATGCATGCATGC"` | A=T=G=C=4 | $p = 1/4$ | $-4 \times (0.25 \times \log_2 0.25) = -4 \times (-0.5)$ | $2.0$ | ✓ yes |
| `"AAAAAAA"` | A=7 | $p = 1$ | $-(1 \times \log_2 1) = 0$ | $0.0$ | ✓ yes |
| `"ATATATAT"` | A=T=4 | $p = 1/2$ | $-2 \times (0.5 \times \log_2 0.5) = -2 \times (-0.5)$ | $1.0$ | ✓ yes |
| `"ATGATGATG"` | A=T=G=3 | $p = 1/3$ | $-3 \times (\frac{1}{3} \log_2 \frac{1}{3})$ | $\log_2(3) \approx 1.58496$ | ✗ tolerance |
| `""` | — | — | empty sum | $0.0$ | ✓ yes |
| `"ATGCNN"` | A=T=G=C=1 | $p = 1/4$ (N ignored) | $-4 \times (0.25 \times \log_2 0.25)$ | $2.0$ | ✓ yes |
| `"NNNNNN"` | total=0 | — | total=0 → 0 | $0.0$ | ✓ yes |

### K-mer Entropy

| Input | k | K-mers (count) | Total | Formula | H |
|-------|---|----------------|-------|---------|---|
| `"AAAAAAAAAA"` | 2 | AA(9) | 9 | $-(1 \cdot \log_2 1)$ | $0.0$ exact |
| `"ATCG"` | 2 | AT(1), TC(1), CG(1) | 3 | $-3 \times (\frac{1}{3} \log_2 \frac{1}{3})$ | $\log_2(3)$ |
| `"ATGCATGCATGCATGC"` | 2 | AT(4), TG(4), GC(4), CA(3) | 15 | $-3 \cdot \frac{4}{15} \log_2 \frac{4}{15} - \frac{3}{15} \log_2 \frac{3}{15}$ | $\approx 1.98082$ |
| `"AT"` | 5 | — | 0 | empty (length < k) | $0.0$ |

---

## Coverage Classification

### Canonical (`SequenceComplexityTests.cs`) — 17 test methods

| # | Test Method | Spec ID | Status |
|---|-------------|---------|--------|
| 1 | `CalculateShannonEntropy_EqualBases_ReturnsTwo` | M01 | ✅ |
| 2 | `CalculateShannonEntropy_SingleBase_ReturnsZero` | M02 | ✅ |
| 3 | `CalculateShannonEntropy_TwoBases_ReturnsOne` | M03 | ✅ |
| 4 | `CalculateShannonEntropy_EmptySequence_ReturnsZero` | M04 | ✅ |
| 5 | `CalculateShannonEntropy_RangeIsZeroToTwo_ForDnaSequences` | M05 | ✅ |
| 6 | `CalculateShannonEntropy_ThreeBases_ReturnsLog2Of3` | S03 | ✅ |
| 7 | `CalculateShannonEntropy_StringOverload_MatchesDnaSequenceOverload` | S02 | ✅ |
| 8 | `CalculateShannonEntropy_LowercaseInput_HandledCorrectly` | S01 | ✅ |
| 9 | `CalculateShannonEntropy_NonDnaCharacters_Ignored` | M08 | ✅ |
| 10 | `CalculateShannonEntropy_NullSequence_ThrowsException` | M09 | ✅ |
| 11 | `CalculateKmerEntropy_VariedDinucleotides_ReturnsExact` | S04 | ✅ |
| 12 | `CalculateKmerEntropy_RepeatedDinucleotides_ReturnsZero` | M07 | ✅ |
| 13 | `CalculateKmerEntropy_SequenceShorterThanK_ReturnsZero` | M06 | ✅ |
| 14 | `CalculateKmerEntropy_InvalidK_ThrowsException` | M10 | ✅ |
| 15 | `CalculateKmerEntropy_NullSequence_ThrowsException` | M09 | ✅ |
| 16 | `CalculateKmerEntropy_RangeIsNonNegativeAndBounded_ForDnaSequences` | INV-004 | ✅ |
| 17 | `CalculateKmerEntropy_UniformDinucleotides_ReturnsLog2Of3` | S05 | ✅ |

### Wrapper smoke (`SequenceStatisticsTests.cs`) — 3 tests
Delegation verification only (uniform → 2.0, homopolymer → 0.0, empty → 0.0).

### Classification Summary
- ✅ Covered: 17 canonical + 3 smoke = 20 total
- ❌ Missing: 0
- ⚠ Weak: 0
- 🔁 Duplicate: 0

---

## Validation Checklist

- [x] All MUST tests have evidence source
- [x] Invariants are mathematically verifiable from Shannon (1948) and Wikipedia
- [x] Edge cases documented (empty, null, invalid k, homopolymer, non-DNA chars)
- [x] API consistency verified (string vs DnaSequence overload)
- [x] Parameter validation tested (null, k < 1)
- [x] Cross-verified against hand calculations using Shannon formula
- [x] Float-exact values use exact assertions; non-exact use $10^{-10}$ tolerance
- [x] Bitwise-identical computations use exact equality (no tolerance)
- [x] No assumptions — all behaviors sourced from Shannon (1948) or Wikipedia
- [x] No duplicates — each test serves a distinct purpose
- [x] Coverage classification complete: 0 missing, 0 weak, 0 duplicate
