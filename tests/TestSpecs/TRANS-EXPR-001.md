# Test Specification: TRANS-EXPR-001

**Test Unit ID:** TRANS-EXPR-001
**Area:** Transcriptome
**Algorithm:** Expression Quantification (TPM, FPKM/RPKM, Quantile Normalization)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Wagner, Kin & Lynch (2012), Theory in Biosciences — TPM introduction | 1 | https://doi.org/10.1007/s12064-012-0162-3 | 2026-06-13 |
| 2 | Zhao, Ye & Stanton (2020), RNA — TPM/RPKM formulas | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC7373998/ | 2026-06-13 |
| 3 | Pimentel (2014) — What the FPKM? (corroborating review) | 3 | https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/ | 2026-06-13 |
| 4 | Wikipedia "Quantile normalization" (citing Bolstad 2003) | 4 | https://en.wikipedia.org/wiki/Quantile_normalization | 2026-06-13 |

### 1.2 Key Evidence Points

1. `TPM_i = (X_i/l_i) / Σ_j(X_j/l_j) · 10^6` — Source 2 (verbatim), Source 3.
2. `FPKM_i = X_i · 10^9 / (l_i · N)` (N = total reads) — Source 2 (`RPKM = 10^9·reads/(total·length)`), Source 3.
3. TPM within a sample sums to exactly 10^6 (average TPM = 10^6 / #transcripts) — Sources 1, 2.
4. Quantile normalization: sort each column, set each rank to the mean across columns, re-place at original positions — Source 4.
5. Tie rule: tied values receive the mean of the rank means they would otherwise span — Source 4 (verbatim).

### 1.3 Documented Corner Cases

- TPM all-zero counts → denominator 0/0, undefined; convention: emit 0 (ASSUMPTION-01).
- Non-positive length / total reads → FPKM undefined; convention: 0.
- Quantile normalization with tied values → averaged rank means (Source 4).
- Empty matrix → no rank means; empty output.

### 1.4 Known Failure Modes / Pitfalls

1. TPM/RPKM are within-sample relative measures; misused across samples/protocols — Sources 1, 2.
2. Quantile normalization that assigns distinct ranks to ties (ignoring the tie-average rule) produces incorrect values — Source 4.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateTPM(geneCounts)` | TranscriptomeAnalyzer | **Canonical** | TPM per Source 2 formula |
| `CalculateFPKM(rawCount, length, totalReads)` | TranscriptomeAnalyzer | **Canonical** | FPKM per Source 2; promoted to public per Registry signature |
| `QuantileNormalize(samples)` | TranscriptomeAnalyzer | **Canonical** | Bolstad 2003 with tie-average rule |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | TPM values of a non-degenerate sample sum to 10^6 | Yes | Sources 1, 2 |
| INV-2 | TPM_i ≥ 0 for all i; equal RPK ⇒ equal TPM | Yes | Source 2 formula |
| INV-3 | FPKM_i ≥ 0; FPKM = 0 when length ≤ 0 or N ≤ 0 | Yes | Source 2; degenerate convention |
| INV-4 | Quantile normalization preserves rank order within each column (monotone) | Yes | Source 4 (rank → rank-mean is non-decreasing) |
| INV-5 | After QN, every column is a permutation of the same multiset of rank means (untied case) | Yes | Source 4 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | TPM three-gene | A(10,2000),B(20,4000),C(30,1000) | TPM = 125000,125000,750000 | Source 2/3 formula derivation |
| M2 | TPM sums to million | M1 result summed | 1,000,000 (±1e-6) | Sources 1,2 (INV-1) |
| M3 | TPM equal RPK ⇒ equal TPM | A,B both RPK 0.005 | TPM_A = TPM_B = 125000 | Source 2 (INV-2) |
| M4 | FPKM single gene | X=1000,l=2000,N=10^6 | FPKM = 500 | Source 2/3 formula |
| M5 | FPKM scales linearly with count | X=2000 (else as M4) | FPKM = 1000 | Source 2 formula |
| M6 | Quantile norm worked example | C1=(5,2,3,4),C2=(4,1,4,2),C3=(3,4,6,8) | Output matrix per Evidence (incl. tie) | Source 4 |
| M7 | Quantile norm tie handling | C2=(4,1,4,2) tied 4s | C2 = (31/6, 2, 31/6, 3); both tied 4s = 5.17 | Source 4 tie rule |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | TPM empty input | no genes | empty sequence | degenerate |
| S2 | FPKM non-positive length | l=0 | 0 | degenerate convention |
| S3 | FPKM non-positive total | N=0 | 0 | degenerate convention |
| S4 | QN empty input | no samples | empty sequence | degenerate |
| S5 | QN rank-order preserved | random column | output non-decreasing where input increasing | INV-4 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | TPM all-zero counts | every count 0 | all TPM = 0 | ASSUMPTION-01 |
| C2 | QN identical columns | all columns equal | output equals input | mean of equal values |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Existing fixture: `tests/Seqeron/Seqeron.Genomics.Tests/TranscriptomeAnalyzerTests.cs` — checked for TPM/FPKM/QuantileNormalize coverage.
- No `{Class}_{Method}_Tests.cs` file for this unit existed prior to this work.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new canonical file |
| M2 | ❌ Missing | |
| M3 | ❌ Missing | |
| M4 | ❌ Missing | |
| M5 | ❌ Missing | |
| M6 | ❌ Missing | |
| M7 | ❌ Missing | exposes tie-handling defect in existing impl |
| S1 | ❌ Missing | |
| S2 | ❌ Missing | |
| S3 | ❌ Missing | |
| S4 | ❌ Missing | |
| S5 | ❌ Missing | |
| C1 | ❌ Missing | |
| C2 | ❌ Missing | |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/TranscriptomeAnalyzer_ExpressionQuantification_Tests.cs` — all M/S/C cases for this unit.
- **Remove:** none. Pre-existing `TranscriptomeAnalyzerTests.cs` covers other methods and is out of scope.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `TranscriptomeAnalyzer_ExpressionQuantification_Tests.cs` | Canonical for TRANS-EXPR-001 | 14 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented (after impl fix) | ✅ Done |
| 8 | S1 | ❌ Missing | Implemented | ✅ Done |
| 9 | S2 | ❌ Missing | Implemented | ✅ Done |
| 10 | S3 | ❌ Missing | Implemented | ✅ Done |
| 11 | S4 | ❌ Missing | Implemented | ✅ Done |
| 12 | S5 | ❌ Missing | Implemented | ✅ Done |
| 13 | C1 | ❌ Missing | Implemented | ✅ Done |
| 14 | C2 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 14
**✅ Done:** 14 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | exact values |
| M2 | ✅ Covered | INV-1 |
| M3 | ✅ Covered | INV-2 |
| M4 | ✅ Covered | exact value |
| M5 | ✅ Covered | linearity |
| M6 | ✅ Covered | exact matrix |
| M7 | ✅ Covered | tie rule (impl corrected) |
| S1 | ✅ Covered | empty |
| S2 | ✅ Covered | length ≤ 0 |
| S3 | ✅ Covered | total ≤ 0 |
| S4 | ✅ Covered | empty |
| S5 | ✅ Covered | INV-4 |
| C1 | ✅ Covered | all-zero |
| C2 | ✅ Covered | identical columns |

**In-scope cases:** 14 | **✅:** 14

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| ASSUMPTION-01 | All-zero TPM denominator → TPM = 0 for every gene (0/0 undefined; degenerate convention) | C1 |
| ASSUMPTION-02 | Effective length = annotated length (no fragment-length correction) | TPM/FPKM formulas |

---

## 7. Open Questions / Decisions

1. **Decision:** `CalculateFPKM` was private in the existing class but the Registry lists it as a canonical method `CalculateFPKM(count, length, total)`. It is promoted to a public static method (in-scope API correction) so it can be tested directly.
2. **Decision (defect fixed):** the existing `QuantileNormalize` assigned distinct ranks to tied values, violating the Bolstad tie-average rule (Source 4). The implementation is corrected to average rank means over tied positions.
