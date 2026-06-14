# Test Specification: ONCO-CNA-001

**Test Unit ID:** ONCO-CNA-001
**Area:** Oncology
**Algorithm:** Copy-Number Alteration Classification (log2 copy ratio → absolute copy number → CNA state)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Mermel et al. (2011), GISTIC2.0, Genome Biology 12(4):R41 | 1 | https://doi.org/10.1186/gb-2011-12-4-r41 | 2026-06-14 |
| 2 | CNVkit `cnvlib/call.py` (`absolute_threshold`, `_log2_ratio_to_absolute_pure`) | 3 | https://raw.githubusercontent.com/etal/cnvkit/master/cnvlib/call.py | 2026-06-14 |
| 3 | CNVkit docs — `call` threshold method | 3 | https://cnvkit.readthedocs.io/en/stable/pipeline.html | 2026-06-14 |
| 4 | GISTIC2 docs — `-ta`/`-td` thresholds | 2 | https://broadinstitute.github.io/gistic2/ | 2026-06-14 |

### 1.2 Key Evidence Points

1. Absolute copy number from a log2 ratio (pure, diploid): `n = ploidy · 2^log2 = 2 · 2^log2` — CNVkit `_log2_ratio_to_absolute_pure`.
2. Hard-threshold integer calling: CN = index of the first threshold the log2 value is `<=`; above the last threshold CN = `ceil(2 · 2^log2)` — CNVkit `absolute_threshold`.
3. Default thresholds `(-1.1, -0.25, 0.2, 0.7)` → states `[0, 1, 2, 3, 4+]` — CNVkit `do_call`.
4. Verbatim cutoffs: `DEL(0) < -1.1`, `LOSS(1) < -0.25`, `GAIN(3) >= +0.2`, `AMP(4) >= +0.7` — CNVkit `absolute_threshold` docstring.
5. GISTIC2 ±0.1 noise band and high-amplitude amp/del thresholds (0.848 / −0.737) corroborate a neutral band bounded by amplification (high positive) and deletion (high negative) — Mermel et al. (2011).

### 1.3 Documented Corner Cases

- Boundary inclusivity: `log2 <= thresh`, so a value exactly on a threshold gets the LOWER CN state of that bin (CNVkit binning loop).
- NaN log2 ratio: no-call → neutral reference copy number (CN 2) (CNVkit `absolute_threshold`).
- Above last threshold: CN grows as `ceil(2·2^log2)`, not a fixed value (CNVkit `absolute_threshold`).

### 1.4 Known Failure Modes / Pitfalls

1. Using simple rounding (`round(2·2^log2)`) instead of the threshold bins would mis-call values near bin edges — CNVkit uses hard thresholds, not rounding, for state calling.
2. Treating thresholds as exclusive (`<`) would mis-classify exact-boundary log2 ratios — the source uses `<=`.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `Log2RatioToCopyNumber(log2Ratio, ploidy)` | OncologyAnalyzer | Canonical | `n = ploidy·2^log2`; continuous absolute copy number |
| `CallCopyNumber(log2Ratio, thresholds, ploidy)` | OncologyAnalyzer | Canonical | hard-threshold integer CN (CNVkit `absolute_threshold`) |
| `ClassifyCopyNumber(log2Ratio, thresholds, ploidy)` | OncologyAnalyzer | Canonical | 5-state CNA classification + integer CN |
| `ClassifyCopyNumbers(log2Ratios, ...)` | OncologyAnalyzer | Delegate | per-element map over `ClassifyCopyNumber`; order/length preserving |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | `Log2RatioToCopyNumber(0, 2) = 2.0` (neutral diploid) | Yes | CNVkit `_log2_ratio_to_absolute_pure` |
| INV-2 | Called integer CN is monotonically non-decreasing in log2 ratio | Yes | CNVkit `absolute_threshold` (thresholds ascending; ceil monotone) |
| INV-3 | Integer CN ≥ 0 for all finite log2 ratios | Yes | CNVkit (CN 0 is the lowest state) |
| INV-4 | State ↔ CN mapping: CN 0→DeepDeletion, 1→Loss, 2→Neutral, 3→Gain, ≥4→Amplification | Yes | CNVkit `absolute_threshold` docstring |
| INV-5 | `ClassifyCopyNumbers` output length = input length, element i ↔ input i | Yes | per-element deterministic map |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Deep deletion | log2 = −2.0 | CN 0, DeepDeletion | CNVkit `≤ −1.1 → 0` |
| M2 | Loss | log2 = −1.0 | CN 1, Loss | CNVkit `(−1.1, −0.25] → 1` |
| M3 | Neutral | log2 = 0.0 | CN 2, Neutral | CNVkit `(−0.25, 0.2] → 2` |
| M4 | Gain | log2 = log2(3/2)=0.5849625 | CN 3, Gain | CNVkit `(0.2, 0.7] → 3` |
| M5 | Amplification | log2 = 1.0 | CN 4, Amplification | CNVkit else `ceil(2·2^1)=4` |
| M6 | Amplification high | log2 = 2.0 | CN 8, Amplification | CNVkit `ceil(2·2^2)=8` |
| M7 | Boundary −1.1 | log2 = −1.1 | CN 0, DeepDeletion | `log2 <= thresh` (inclusive) |
| M8 | Boundary −0.25 | log2 = −0.25 | CN 1, Loss | inclusive |
| M9 | Boundary 0.2 | log2 = 0.2 | CN 2, Neutral | inclusive |
| M10 | Boundary 0.7 | log2 = 0.7 | CN 3, Gain | inclusive |
| M11 | Absolute CN formula | log2 0,1,−1,log2(3/2) | n = 2.0, 4.0, 1.0, 3.0 | `_log2_ratio_to_absolute_pure` |
| M12 | NaN log2 | log2 = NaN | CN 2, Neutral (no-call) | CNVkit NaN→neutral |
| M13 | Batch classify | {−2,0,1} | [DeepDeletion, Neutral, Amplification], length 3, order preserved | per-element map |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Custom thresholds | thresholds (−0.4,−0.1,0.1,0.4), log2 = −0.3 | CN 1, Loss | CNVkit germline-tuned alt |
| S2 | Monotonicity (INV-2) | ascending log2 sequence | CN non-decreasing | property check |
| S3 | Amplification ceil | log2 = 0.8 | CN ceil(2·2^0.8)=4, Amplification | `ceil` not `round` |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Empty batch | ClassifyCopyNumbers([]) | empty list | trivial |

### 4.4 Error / Validation Tests

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| E1 | Null thresholds | thresholds = null | ArgumentNullException | input validation |
| E2 | Wrong threshold count | thresholds length ≠ 4 | ArgumentException | 4 cutoffs define 5 states |
| E3 | Non-ascending thresholds | thresholds not strictly ascending | ArgumentException | bins must be ordered |
| E4 | Non-positive ploidy | ploidy ≤ 0 | ArgumentOutOfRangeException | n = ploidy·2^log2 needs ploidy > 0 |
| E5 | Null batch | ClassifyCopyNumbers(null) | ArgumentNullException | input validation |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No prior tests for CNA classification in `OncologyAnalyzer`. SV-CNV-001 tests cover `StructuralVariantAnalyzer.DetectCNV`/`SegmentCopyNumber` (integer CN via rounding), which is a separate unit — no overlap with the 5-state classification.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M13 | ❌ Missing | new unit; no existing tests |
| S1–S3 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |
| E1–E5 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_CopyNumberClassification_Tests.cs` — all ONCO-CNA-001 cases.
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| OncologyAnalyzer_CopyNumberClassification_Tests.cs | Canonical, all cases | 22 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | implemented | ✅ Done |
| 2 | M2 | ❌ Missing | implemented | ✅ Done |
| 3 | M3 | ❌ Missing | implemented | ✅ Done |
| 4 | M4 | ❌ Missing | implemented | ✅ Done |
| 5 | M5 | ❌ Missing | implemented | ✅ Done |
| 6 | M6 | ❌ Missing | implemented | ✅ Done |
| 7 | M7 | ❌ Missing | implemented | ✅ Done |
| 8 | M8 | ❌ Missing | implemented | ✅ Done |
| 9 | M9 | ❌ Missing | implemented | ✅ Done |
| 10 | M10 | ❌ Missing | implemented | ✅ Done |
| 11 | M11 | ❌ Missing | implemented | ✅ Done |
| 12 | M12 | ❌ Missing | implemented | ✅ Done |
| 13 | M13 | ❌ Missing | implemented | ✅ Done |
| 14 | S1 | ❌ Missing | implemented | ✅ Done |
| 15 | S2 | ❌ Missing | implemented | ✅ Done |
| 16 | S3 | ❌ Missing | implemented | ✅ Done |
| 17 | C1 | ❌ Missing | implemented | ✅ Done |
| 18 | E1 | ❌ Missing | implemented | ✅ Done |
| 19 | E2 | ❌ Missing | implemented | ✅ Done |
| 20 | E3 | ❌ Missing | implemented | ✅ Done |
| 21 | E4 | ❌ Missing | implemented | ✅ Done |
| 22 | E5 | ❌ Missing | implemented | ✅ Done |

**Total items:** 22
**✅ Done:** 22 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1–M13 | ✅ Covered | exact evidence-based asserts |
| S1–S3 | ✅ Covered | edge cases + monotonicity property |
| C1 | ✅ Covered | empty batch |
| E1–E5 | ✅ Covered | validation throws |

Total in-scope cases: 22. ✅: 22.

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Diploid (autosomal) reference ploidy = 2 (CNVkit default; sex/haploid out of scope) | Log2RatioToCopyNumber, CallCopyNumber, ClassifyCopyNumber |

---

## 7. Open Questions / Decisions

1. Default thresholds use the CNVkit **source-code** tumor defaults `(-1.1, -0.25, 0.2, 0.7)` (rank-3 reference implementation), not the docs-page germline-tuned variant `(-1.1, -0.4, 0.3, 0.7)`; callers can override via the thresholds parameter. Decision recorded; no open issue.
