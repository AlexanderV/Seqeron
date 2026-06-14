# Test Specification: SEQ-COMPLEX-KMER-001

**Test Unit ID:** SEQ-COMPLEX-KMER-001
**Area:** Complexity
**Algorithm:** K-mer Entropy (Shannon entropy of the overlapping k-mer frequency distribution)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Li, H. (2025). Finding low-complexity DNA sequences with longdust. arXiv:2509.07357 | 1 | https://arxiv.org/pdf/2509.07357 | 2026-06-14 |
| 2 | Çakır et al. (2025). Entropy–Rank Ratio. arXiv:2511.05300 | 1 | https://arxiv.org/html/2511.05300 | 2026-06-14 |
| 3 | Shannon, C.E. (1948) A Mathematical Theory of Communication (via citing secondaries) | 4 | https://en.wikipedia.org/wiki/Entropy_(information_theory) ; https://tcosmo.github.io/2019/04/21/shannon-entropy.html | 2026-06-14 |

### 1.2 Key Evidence Points

1. H = −Σ p_i log₂(p_i) where p_i is the frequency of the i-th k-mer — Li 2025.
2. K-mers are overlapping (sliding window, step 1); a length-L sequence has N = L−k+1 k-mers; p_i = n_i/N — Li 2025.
3. Logarithm base 2 → entropy in bits; single-nucleotide max = log₂(4) = 2 bits — Çakır 2025.
4. Bounds: 0 ≤ H ≤ log_b(n); H = 0 for a deterministic distribution; H = log_b(n) for uniform over n symbols — Shannon 1948 (via secondaries).

### 1.3 Documented Corner Cases

- Single repeated k-mer (homopolymer / tandem repeat) → p=1 → H = 0 (Li 2025).
- All k-mers distinct → p_i = 1/N → H = log₂(N) (Li 2025; Shannon uniform bound).

### 1.4 Known Failure Modes / Pitfalls

1. Confusing non-overlapping tuples (Çakır 2025, M=⌊L/n⌋) with overlapping k-mers (Li 2025, N=L−k+1). This unit uses the **overlapping** convention. — Li 2025 vs Çakır 2025.
2. Wrong log base (must be 2 → bits). — Çakır 2025.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateKmerEntropy(DnaSequence, int k = 2)` | SequenceComplexity | **Canonical** | Deep evidence-based testing |
| `CalculateKmerEntropy(string, int k = 2)` | SequenceComplexity | **Delegate** | String overload; upper-cases then delegates to same core; smoke + agreement |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | For any valid input, 0 ≤ H ≤ log₂(N) where N = L−k+1 (≥ 0; equals log₂ of the number of distinct k-mers' upper bound). | Yes | Shannon bounds (src 3); Li 2025 |
| INV-2 | A single distinct k-mer (deterministic distribution) ⇒ H = 0. | Yes | Shannon (src 3); Li 2025 |
| INV-3 | All-distinct k-mers (uniform) ⇒ H = log₂(N). | Yes | Shannon uniform bound (src 3); Li 2025 |
| INV-4 | Result is independent of letter case (upper/lower agree), since DnaSequence normalizes to upper-case. | Yes | DnaSequence contract; src 1 formula symbol-agnostic |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | `CalculateKmerEntropy_UniformMonomers_ReturnsLog2Of4` | `ACGT`, k=1 — 4 distinct monomers, uniform | 2.0 | Çakır 2025 (max=log₂4); Shannon uniform bound |
| M2 | `CalculateKmerEntropy_AllDistinctDimers_ReturnsLog2OfN` | `ACGT`, k=2 — 3 distinct dimers, uniform | log₂(3) = 1.5849625007211562 | Li 2025 all-distinct; Shannon uniform |
| M3 | `CalculateKmerEntropy_NonUniformDimers_ReturnsExact` | `ATATAT`, k=2 — AT=3,TA=2 (binary entropy of 0.6) | 0.9709505944546686 | Li 2025 formula H=−Σ p log₂ p |
| M4 | `CalculateKmerEntropy_SingleRepeatedDimer_ReturnsZero` | `AAAA`, k=2 — only AA (deterministic) | 0.0 | Shannon H=0 for certainty; Li 2025 |
| M5 | `CalculateKmerEntropy_MixedCounts_ReturnsExact` | `AAACGT`, k=2 — AA=2,AC=1,CG=1,GT=1 | 1.9219280948873623 (= log₂5 − 0.4) | Li 2025 formula |
| M6 | `CalculateKmerEntropy_SequenceShorterThanK_ReturnsZero` | `AC`, k=5 — no k-mers | 0.0 | No-k-mers boundary (ASSUMPTION; sibling contract) |
| M7 | `CalculateKmerEntropy_InvalidK_Throws` | k=0 | `ArgumentOutOfRangeException` | Contract (ASSUMPTION; sibling guards) |
| M8 | `CalculateKmerEntropy_NullDnaSequence_Throws` | null DnaSequence | `ArgumentNullException` | Contract (ASSUMPTION; sibling guards) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | `CalculateKmerEntropy_StringOverload_AgreesWithDnaSequence` | `ATATAT` string vs DnaSequence, k=2 | both 0.9709505944546686 | Delegate smoke + INV-4 |
| S2 | `CalculateKmerEntropy_LowercaseString_EqualsUppercase` | `atatat` vs `ATATAT`, k=2 | equal | INV-4 (case normalization) |
| S3 | `CalculateKmerEntropy_NullOrEmptyString_ReturnsZero` | null and "" string | 0.0 | String overload contract |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | `CalculateKmerEntropy_BoundsInvariant_WithinRange` | several sequences/k | 0 ≤ H ≤ log₂(L−k+1) | INV-1 property test |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/SequenceComplexityTests.cs` contains pre-existing `CalculateKmerEntropy_*` tests (under prior units SEQ-COMPLEX-001 / SEQ-ENTROPY-001), but no canonical `SequenceComplexity_CalculateKmerEntropy_Tests.cs` file existed for this unit.
- No prior tests cover the new `CalculateKmerEntropy(string, k)` overload (added by this unit).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | New canonical file |
| M2 | ❌ Missing | |
| M3 | ❌ Missing | |
| M4 | ❌ Missing | |
| M5 | ❌ Missing | |
| M6 | ❌ Missing | |
| M7 | ❌ Missing | |
| M8 | ❌ Missing | |
| S1 | ❌ Missing | string overload new |
| S2 | ❌ Missing | |
| S3 | ❌ Missing | |
| C1 | ❌ Missing | invariant property |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceComplexity_CalculateKmerEntropy_Tests.cs` — all M/S/C cases for this unit.
- **Remove:** nothing. Pre-existing `CalculateKmerEntropy_*` tests in `SequenceComplexityTests.cs` belong to earlier units (SEQ-COMPLEX-001 / SEQ-ENTROPY-001) and are left untouched (out of scope for this unit).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceComplexity_CalculateKmerEntropy_Tests.cs` | Canonical for SEQ-COMPLEX-KMER-001 | 12 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented exact-value test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented exact-value test | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented exact-value test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented exact-value test | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented exact-value test | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented boundary test | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented throw test | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented throw test | ✅ Done |
| 9 | S1 | ❌ Missing | Implemented agreement test | ✅ Done |
| 10 | S2 | ❌ Missing | Implemented case test | ✅ Done |
| 11 | S3 | ❌ Missing | Implemented null/empty test | ✅ Done |
| 12 | C1 | ❌ Missing | Implemented bounds invariant property test | ✅ Done |

**Total items:** 12
**✅ Done:** 12 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | Exact 2.0 |
| M2 | ✅ Covered | Exact log₂3 |
| M3 | ✅ Covered | Exact 0.9709505944546686 |
| M4 | ✅ Covered | Exact 0.0 |
| M5 | ✅ Covered | Exact 1.9219280948873623 |
| M6 | ✅ Covered | 0.0 (L<k) |
| M7 | ✅ Covered | Throws ArgumentOutOfRangeException |
| M8 | ✅ Covered | Throws ArgumentNullException |
| S1 | ✅ Covered | string == DnaSequence |
| S2 | ✅ Covered | lowercase == uppercase |
| S3 | ✅ Covered | null/empty → 0 |
| C1 | ✅ Covered | 0 ≤ H ≤ log₂N |

**Total in-scope cases:** 12 | **✅:** 12

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | L < k returns 0 (no k-mers; entropy of empty multiset = 0). Resolved by sibling contract. | M6 |
| 2 | Invalid k (<1) → ArgumentOutOfRangeException; null DnaSequence → ArgumentNullException; null/empty string → 0. API contract, matches sibling guards. | M7, M8, S3 |

---

## 7. Open Questions / Decisions

1. **Overlapping vs non-overlapping k-mers** — Decision: overlapping (Li 2025, N=L−k+1), which the existing implementation already uses; Çakır 2025's non-overlapping tuples are an alternative convention not adopted here. Documented in §1.4.
