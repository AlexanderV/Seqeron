# Test Specification: SEQ-GC-PROFILE-001

**Test Unit ID:** SEQ-GC-PROFILE-001
**Area:** Statistics
**Algorithm:** GC Content Profile (sliding-window GC content)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Wikipedia — GC-content | 4 | https://en.wikipedia.org/wiki/GC-content | 2026-06-14 |
| 2 | Biopython `Bio.SeqUtils.gc_fraction` (source) | 3 | https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py | 2026-06-14 |

### 1.2 Key Evidence Points

1. GC content = `(G + C) / (A + T + G + C) × 100%` — Wikipedia GC-content (source 1).
2. Denominator is the standard bases A+T+G+C; N is excluded — Wikipedia (source 1); Biopython `gc_fraction(..., "remove")` default, `"ACTGN" → 0.50` (source 2).
3. `gc_fraction` returns a fraction in [0, 1]; ×100 gives the percentage. `GC123` returns "percentages between 0 and 100", confirming the percentage form (source 2).
4. RNA U behaves as a non-GC base like T: `gc_fraction("GGAUCUUCGGAUCU") → 0.50` (source 2).

### 1.3 Documented Corner Cases

- Ambiguous N: `remove` excludes from denominator (`ACTGN → 0.50`), `ignore` includes it (`→ 0.40`) — Biopython (source 2). This unit uses the `remove` convention.
- A+T+G+C = 0 (window with no standard base) → GC content undefined (division by zero); repository convention returns 0 (Assumption A1).

### 1.4 Known Failure Modes / Pitfalls

1. Counting N in the denominator deflates GC% (the `ignore` value, not the `remove` value) — Biopython (source 2).
2. Returning a fraction (0–1) instead of a percentage (0–100) — divergence from Wikipedia and the repository's locked SEQ-GC-ANALYSIS-001 convention.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateGcContentProfile(string sequence, int windowSize = 100, int stepSize = 1)` | SequenceStatistics | **Canonical** | Sliding-window GC%; deep evidence-based testing |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-01 | Each window value = `(G+C)/(A+T+G+C)×100`, in [0, 100] | Yes | Wikipedia formula (source 1) |
| INV-02 | Denominator counts only A/T/U/G/C (N and other symbols excluded) | Yes | Biopython `remove` default `ACTGN → 0.50` (source 2) |
| INV-03 | Window count = ⌊(n − w)/step⌋ + 1 for w ≤ n; 0 otherwise; offsets 0, step, 2·step, … | Yes | Sliding-window definition (source 1 windowed application) |
| INV-04 | U is treated as a non-GC base equivalent to T | Yes | Biopython RNA `GGAUCUUCGGAUCU → 0.50` (source 2) |
| INV-05 | A window with no A/T/U/G/C base yields 0 | Yes | **ASSUMPTION** A1 (repository convention) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | All-GC window | `GGGGGGGGGG` w=10 → 100% | profile = [100.0] | Wikipedia: 10/10×100 (source 1) |
| M2 | Half-GC window | `ATGCATGCAT...` ATGC repeats, w=4 step=4 → 50% each | every value = 50.0 | Wikipedia: 2/4×100 (source 1) |
| M3 | Exact mixed profile | `GGGAAATGCC` w=4 step=3 → windows GGGA, AAAT, TGCC | [75.0, 0.0, 75.0] | Wikipedia formula per window (source 1) |
| M4 | N excluded from denominator | `GGAN` w=4 → 2/3×100 | [66.66666666666666] | Biopython `ACTGN` remove → 0.50 (source 2) |
| M5 | RNA U is non-GC | `GGAU` w=4 → 2/3×... no: 2 GC over 4 bases = 50% | [50.0] | Biopython RNA → 0.50 (source 2) |
| M6 | Window count & offsets | n=10,w=4,step=3 → 3 windows; step=1 → 7 | counts 3 and 7 | Sliding-window (source 1) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | windowSize > length | `ATGC` w=100 | empty profile | INV-03 |
| S2 | windowSize == length | `GGCC` w=4 | [100.0] (single window) | INV-03 |
| S3 | null / empty | null, "" | empty profile | guarded input |
| S4 | all-N window | `NNNN` w=4 | [0.0] | INV-05 / Assumption A1 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Case-insensitivity | `ggggcccc` vs uppercase | equal profiles | case-folded counting |
| C2 | Bound [0,100] | mixed sequence | all values in [0,100] | INV-01 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatisticsTests.cs` (lines 408–444): three pre-template tests for `CalculateGcContentProfile` (`_ReturnsCorrectCount`, `_AllGc_Returns100Percent`, `_WindowTooLarge_ReturnsEmpty`). They use permissive assertions (`GreaterThan(0)`, `>= 0.99`) and assume the fraction (0–1) output.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 (all-GC %) | ⚠ Weak | `_AllGc_Returns100Percent` uses `>= 0.99` (fraction); no exact 100.0 |
| M2 (half-GC) | ❌ Missing | not covered |
| M3 (mixed profile) | ❌ Missing | not covered |
| M4 (N excluded) | ❌ Missing | not covered |
| M5 (RNA U) | ❌ Missing | not covered |
| M6 (window count) | ⚠ Weak | `_ReturnsCorrectCount` only asserts `> 0` |
| S1 (window > length) | ✅ Covered | `_WindowTooLarge_ReturnsEmpty` (kept logic, re-implemented in canonical file) |
| S2 (window == length) | ❌ Missing | not covered |
| S3 (null/empty) | ❌ Missing | not covered |
| S4 (all-N) | ❌ Missing | not covered |
| C1 (case) | ❌ Missing | not covered |
| C2 (bound) | ❌ Missing | not covered |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateGcContentProfile_Tests.cs` — all M/S/C cases with exact evidence-based values.
- **Remove:** the three legacy `CalculateGcContentProfile_*` tests from `SequenceStatisticsTests.cs` (weak/duplicate; superseded by the canonical file).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceStatistics_CalculateGcContentProfile_Tests.cs` | Canonical | 12 |
| `SequenceStatisticsTests.cs` | Legacy profile tests removed | 0 (for this unit) |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ⚠ Weak | rewrite exact 100.0 | ✅ Done |
| 2 | M2 | ❌ Missing | implement | ✅ Done |
| 3 | M3 | ❌ Missing | implement | ✅ Done |
| 4 | M4 | ❌ Missing | implement | ✅ Done |
| 5 | M5 | ❌ Missing | implement | ✅ Done |
| 6 | M6 | ⚠ Weak | rewrite exact counts | ✅ Done |
| 7 | S1 | ✅ Covered | re-implement in canonical file | ✅ Done |
| 8 | S2 | ❌ Missing | implement | ✅ Done |
| 9 | S3 | ❌ Missing | implement | ✅ Done |
| 10 | S4 | ❌ Missing | implement | ✅ Done |
| 11 | C1 | ❌ Missing | implement | ✅ Done |
| 12 | C2 | ❌ Missing | implement | ✅ Done |
| 13 | legacy tests | 🔁 Duplicate | removed from SequenceStatisticsTests.cs | ✅ Done |

**Total items:** 13
**✅ Done:** 13 | **⛔ Blocked:** 0 | **Remaining:** must be 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | exact 100.0 |
| M2 | ✅ | exact 50.0 |
| M3 | ✅ | [75.0, 0.0, 75.0] |
| M4 | ✅ | 66.66666666666666 |
| M5 | ✅ | 50.0 (U non-GC) |
| M6 | ✅ | counts 3 and 7 |
| S1 | ✅ | empty |
| S2 | ✅ | [100.0] |
| S3 | ✅ | empty |
| S4 | ✅ | [0.0] |
| C1 | ✅ | equal profiles |
| C2 | ✅ | all in [0,100] |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| A1 | A window with no A/T/U/G/C base returns 0 (repository convention, matches SEQ-GC-ANALYSIS-001) | INV-05, S4 |

---

## 7. Open Questions / Decisions

1. **Decision (resolved):** the implementation previously returned a fraction (0–1). It is corrected to return a percentage (×100) to conform to the Wikipedia GC-content definition and the repository's locked SEQ-GC-ANALYSIS-001 convention; checklist behavioral note ("(G+C)/window×100") agrees.
