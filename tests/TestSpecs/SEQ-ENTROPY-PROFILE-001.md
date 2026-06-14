# Test Specification: SEQ-ENTROPY-PROFILE-001

**Test Unit ID:** SEQ-ENTROPY-PROFILE-001
**Area:** Statistics
**Algorithm:** Shannon Entropy Profile (sliding-window Shannon entropy)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Shannon C. E. (1948). A Mathematical Theory of Communication. Bell Syst. Tech. J. 27(3):379–423. | 1 | https://doi.org/10.1002/j.1538-7305.1948.tb01338.x | 2026-06-14 |
| 2 | Wikipedia — Entropy (information theory) (citing Shannon 1948) | 4 | https://en.wikipedia.org/wiki/Entropy_(information_theory) | 2026-06-14 |
| 3 | Entropy-Based Biological Sequence Study, IntechOpen | 3–4 | https://www.intechopen.com/chapters/75997 | 2026-06-14 |

### 1.2 Key Evidence Points

1. Shannon entropy H(X) = −Σᵢ pᵢ log_b pᵢ; with b = 2 the unit is bits (shannons) — Source 2 (citing Source 1).
2. Maximum entropy of n equally-likely outcomes is log₂(n); attained by the uniform distribution — Source 2.
3. For the 4-letter DNA alphabet maximum entropy is 2 bits (log₂4); pᵢ are per-symbol frequencies, computed per sliding window of width W — Source 3 (Eq. 3 yᵢ = −Σⱼ pᵢⱼ log pᵢⱼ).
4. Zero-probability terms contribute 0 (0·log0 ≡ 0); a homopolymer window has H = 0 — Source 1/2.

### 1.3 Documented Corner Cases

- Homopolymer / single-symbol window → H = 0 (Source 1/2, zero-probability convention).
- Uniform window over k symbols → H = log₂ k (Source 2 maximum-entropy property).
- Window width W must not exceed sequence length for a full window to exist (Source 3 sliding-window method).

### 1.4 Known Failure Modes / Pitfalls

1. Using log base 10 or e instead of base 2 changes units (Source 2 — bits require b=2).
2. Counting over a wrong alphabet (block/k-mer vs mono-symbol) changes pᵢ and the maximum (Source 3 — 4-letter alphabet, 2-bit max).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateEntropyProfile(string, int, int)` | `SequenceStatistics` | Canonical | Sliding-window profile; yields one H per window |
| `CalculateShannonEntropy(string)` | `SequenceStatistics` | Internal | Per-window kernel exercised via the profile and directly for the H formula |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Each profile value H ≥ 0 | Yes | Source 2 (entropy is non-negative) |
| INV-2 | Each profile value H ≤ log₂ k, where k = number of distinct symbols in the window (≤ 2 for DNA) | Yes | Source 2 (max = log₂ n); Source 3 (2-bit DNA max) |
| INV-3 | A homopolymer window yields H = 0 | Yes | Source 1/2 (zero-probability convention) |
| INV-4 | A window with all symbols equally frequent yields H = log₂ k | Yes | Source 2 (uniform maximizes) |
| INV-5 | Number of windows = ⌊(n − windowSize)/stepSize⌋ + 1 when windowSize ≤ n, else 0 | Yes | Source 3 (sliding window of width W) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Uniform window | `CalculateShannonEntropy("ATGC")` | 2.0 (= log₂4) | Source 2/3 |
| M2 | Two-symbol equal | `CalculateShannonEntropy("AATT")` | 1.0 | Source 2 |
| M3 | Skewed 3:1 | `CalculateShannonEntropy("AAAT")` | 0.8112781244591328 | Source 1 formula |
| M4 | Three-symbol 3:2:1 | `CalculateShannonEntropy("AAATTC")` | 1.4591479170272448 | Source 1 formula |
| M5 | Homopolymer | `CalculateShannonEntropy("AAAA")` | 0.0 | Source 1/2 (INV-3) |
| M6 | Profile step 1 | `CalculateEntropyProfile("AAATGC",4,1)` | [0.8112781244591328, 1.5, 2.0] | Source 3 sliding window + Source 1 formula |
| M7 | Profile step 2 | `CalculateEntropyProfile("AAATGCAA",4,2)` | [0.8112781244591328, 2.0, 1.5] | Source 3 sliding window + Source 1 formula |
| M8 | Window count | profile length for n=8,w=4,step=2 | 3 windows (INV-5) | Source 3 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | windowSize > length | `CalculateEntropyProfile("AAA",4,1)` | empty | INV-5 (no full window) |
| S2 | windowSize == length | `CalculateEntropyProfile("AATT",4,1)` | [1.0] | single window |
| S3 | null / empty input | `CalculateEntropyProfile(null/"" ,4,1)` | empty | guarded input |
| S4 | Max-entropy invariant | every profile value ≤ 2.0 for DNA | holds | INV-2 |
| S5 | Non-negativity invariant | every profile value ≥ 0 | holds | INV-1 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Case-insensitivity | lowercase `"aaatgc"` profile equals uppercase | equal profiles | implementation case-folds |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing test file targets `CalculateEntropyProfile`. `CalculateShannonEntropy` is exercised indirectly elsewhere but no canonical evidence-based fixture exists for this unit.
- Canonical file created: `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateEntropyProfile_Tests.cs`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new unit |
| M2 | ❌ Missing | new unit |
| M3 | ❌ Missing | new unit |
| M4 | ❌ Missing | new unit |
| M5 | ❌ Missing | new unit |
| M6 | ❌ Missing | new unit |
| M7 | ❌ Missing | new unit |
| M8 | ❌ Missing | new unit |
| S1 | ❌ Missing | new unit |
| S2 | ❌ Missing | new unit |
| S3 | ❌ Missing | new unit |
| S4 | ❌ Missing | new unit |
| S5 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateEntropyProfile_Tests.cs` — all MUST/SHOULD/COULD cases.
- **Remove:** none (no pre-existing tests for this unit).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceStatistics_CalculateEntropyProfile_Tests.cs` | Canonical fixture for SEQ-ENTROPY-PROFILE-001 | 14 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented | ✅ Done |
| 9 | S1 | ❌ Missing | Implemented | ✅ Done |
| 10 | S2 | ❌ Missing | Implemented | ✅ Done |
| 11 | S3 | ❌ Missing | Implemented | ✅ Done |
| 12 | S4 | ❌ Missing | Implemented | ✅ Done |
| 13 | S5 | ❌ Missing | Implemented | ✅ Done |
| 14 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 14
**✅ Done:** 14 | **⛔ Blocked:** 0 | **Remaining:** must be 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | exact value asserted |
| M2 | ✅ Covered | exact value asserted |
| M3 | ✅ Covered | exact value asserted |
| M4 | ✅ Covered | exact value asserted |
| M5 | ✅ Covered | exact value asserted |
| M6 | ✅ Covered | exact profile asserted |
| M7 | ✅ Covered | exact profile asserted |
| M8 | ✅ Covered | window count asserted |
| S1 | ✅ Covered | empty profile |
| S2 | ✅ Covered | single window |
| S3 | ✅ Covered | null+empty empty |
| S4 | ✅ Covered | ≤ 2.0 invariant |
| S5 | ✅ Covered | ≥ 0 invariant |
| C1 | ✅ Covered | case-insensitivity |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Per-symbol (mono-nucleotide) alphabet for pᵢ; letters only, case-folded | M1–M8, INV-2 |

---

## 7. Open Questions / Decisions

1. None. The formula, base (2/bits), maximum (log₂k), and zero-probability convention are all source-confirmed; the only modelling choice (mono-symbol alphabet) is recorded as an assumption and matches the 4-letter / 2-bit DNA application in Source 3.
