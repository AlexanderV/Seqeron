# Test Specification: QUALITY-PHRED-001

**Test Unit ID:** QUALITY-PHRED-001
**Area:** Quality
**Algorithm:** Phred Score Handling (FASTQ quality string parse / encode / Phred+33 ↔ Phred+64 conversion)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Cock et al. (2010), *Nucleic Acids Research* 38(6):1767–1771 — Sanger FASTQ format & Solexa/Illumina variants | 1 | https://doi.org/10.1093/nar/gkp1137 (via https://pmc.ncbi.nlm.nih.gov/articles/PMC2847217/) | 2026-06-13 |

### 1.2 Key Evidence Points

1. Sanger / Phred+33: ASCII offset 33; ASCII 33–126 encodes Phred Q 0–93 — Cock et al. (2010).
2. Illumina 1.3+ / Phred+64: ASCII offset 64; ASCII 64–126 encodes Phred Q 0–62 — Cock et al. (2010).
3. Decode: Q = ord(char) − offset; Encode: char = chr(Q + offset) — Cock et al. (2010).
4. Phred score is invariant across the two variants; conversion is a pure byte re-offset (shift by ±31) — Cock et al. (2010).

### 1.3 Documented Corner Cases

- Char below the offset decodes to negative Phred Q → malformed for that variant (Phred Q ≥ 0) — Cock et al. (2010).
- Phred+64 (Q 0–62) → Phred+33 (Q 0–93): always representable — Cock et al. (2010).
- Phred+33 Q ∈ (62, 93] → Phred+64: not representable (overflow) — Cock et al. (2010).

### 1.4 Known Failure Modes / Pitfalls

1. Mis-decoding by applying the wrong offset (off by 31) silently produces wrong scores — Cock et al. (2010).
2. Phred+33→Phred+64 of high-quality scores (Q>62) overflows the Phred+64 range — Cock et al. (2010).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `ParseQualityString(qualStr, encoding)` | QualityScoreAnalyzer | Canonical | Decode ASCII chars → Phred scores with range validation |
| `ToQualityString(scores, encoding)` | QualityScoreAnalyzer | Canonical | Encode Phred scores → ASCII chars with range validation |
| `ConvertEncoding(qualStr, from, to)` | QualityScoreAnalyzer | Canonical | Re-offset chars preserving Phred score |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | For Phred+33, Q = ord(char) − 33, with Q ∈ [0, 93] for ASCII ∈ [33, 126] | Yes | Cock et al. (2010) |
| INV-2 | For Phred+64, Q = ord(char) − 64, with Q ∈ [0, 62] for ASCII ∈ [64, 126] | Yes | Cock et al. (2010) |
| INV-3 | `ToQualityString(ParseQualityString(s, e), e) == s` for any valid quality string `s` under encoding `e` (round-trip identity) | Yes | Cock et al. (2010) points 1–3 |
| INV-4 | `ConvertEncoding` preserves the Phred score: ParseQualityString(out, to) == ParseQualityString(in, from) | Yes | Cock et al. (2010) point 4 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Parse Phred+33 boundaries | `"!~"` decoded as Phred+33 | `[0, 93]` | Cock et al. (2010) pt 1 |
| M2 | Parse Phred+33 interior | `"5?I"` decoded as Phred+33 | `[20, 30, 40]` | Cock et al. (2010) pt 1 |
| M3 | Parse Phred+64 boundaries/interior | `"@h~"` decoded as Phred+64 | `[0, 40, 62]` | Cock et al. (2010) pt 2 |
| M4 | Encode Phred+33 | `[0,20,30,40,93]` → string (Phred+33) | `"!5?I~"` | Cock et al. (2010) pt 3 |
| M5 | Encode Phred+64 | `[0,40,62]` → string (Phred+64) | `"@h~"` | Cock et al. (2010) pt 3 |
| M6 | Convert Phred+64→Phred+33 | `"@h~"` from Phred64 to Phred33 | `"!I_"` (Q 0,40,62 preserved) | Cock et al. (2010) pt 4 |
| M7 | Convert Phred+33→Phred+64 | `"!I"` from Phred33 to Phred64 | `"@h"` (Q 0,40 preserved) | Cock et al. (2010) pt 4 |
| M8 | Parse out-of-range char throws | char below offset (e.g. `" "`=32 as Phred+33 → Q=−1) | `ArgumentOutOfRangeException` | Cock et al. (2010) corner case 1 |
| M9 | Convert Phred+33→Phred+64 overflow throws | `"~"` (Q=93) Phred33→Phred64 | `ArgumentOutOfRangeException` | Cock et al. (2010) corner case 3 |
| M10 | Encode out-of-range score throws | `[94]` as Phred+33 (>93) | `ArgumentOutOfRangeException` | Cock et al. (2010) pt 1 (max Q 93) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Parse empty string | `""` → empty array | `[]` | Trivial boundary |
| S2 | Encode empty scores | `[]` → `""` | `""` | Trivial boundary |
| S3 | Parse null throws | `null` quality string | `ArgumentNullException` | Public-API failure mode |
| S4 | Encode null throws | `null` scores | `ArgumentNullException` | Public-API failure mode |
| S5 | Convert null throws | `null` quality string | `ArgumentNullException` | Public-API failure mode |
| S6 | Convert same encoding identity | Phred33→Phred33 of `"5?I"` | `"5?I"` | No-op shift |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Round-trip property | parse∘encode identity over a fixed-seed valid Phred+33 string | input == output | INV-3 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/QualityScoreAnalyzerTests.cs` — legacy fixture covering `QualityStringToPhred`/`PhredToQualityString` and other analyzer methods; it does NOT reference the canonical `ParseQualityString`, `ToQualityString`, or `ConvertEncoding` methods named for this unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new canonical method, no test |
| M2 | ❌ Missing | new canonical method, no test |
| M3 | ❌ Missing | new canonical method, no test |
| M4 | ❌ Missing | new canonical method, no test |
| M5 | ❌ Missing | new canonical method, no test |
| M6 | ❌ Missing | new canonical method, no test |
| M7 | ❌ Missing | new canonical method, no test |
| M8 | ❌ Missing | new canonical method, no test |
| M9 | ❌ Missing | new canonical method, no test |
| M10 | ❌ Missing | new canonical method, no test |
| S1 | ❌ Missing | new canonical method, no test |
| S2 | ❌ Missing | new canonical method, no test |
| S3 | ❌ Missing | new canonical method, no test |
| S4 | ❌ Missing | new canonical method, no test |
| S5 | ❌ Missing | new canonical method, no test |
| S6 | ❌ Missing | new canonical method, no test |
| C1 | ❌ Missing | new canonical method, no test |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/QualityScoreAnalyzer_ParseQualityString_Tests.cs` — all QUALITY-PHRED-001 cases for the three canonical methods.
- **Remove:** nothing. The legacy `QualityScoreAnalyzerTests.cs` covers other (out-of-scope) methods and is left untouched.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `QualityScoreAnalyzer_ParseQualityString_Tests.cs` | Canonical for QUALITY-PHRED-001 | 17 |
| `QualityScoreAnalyzerTests.cs` | Legacy (other methods, out of scope) | unchanged |

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
| 9 | M9 | ❌ Missing | Implemented | ✅ Done |
| 10 | M10 | ❌ Missing | Implemented | ✅ Done |
| 11 | S1 | ❌ Missing | Implemented | ✅ Done |
| 12 | S2 | ❌ Missing | Implemented | ✅ Done |
| 13 | S3 | ❌ Missing | Implemented | ✅ Done |
| 14 | S4 | ❌ Missing | Implemented | ✅ Done |
| 15 | S5 | ❌ Missing | Implemented | ✅ Done |
| 16 | S6 | ❌ Missing | Implemented | ✅ Done |
| 17 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 17
**✅ Done:** 17 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | exact values asserted |
| M2 | ✅ Covered | exact values asserted |
| M3 | ✅ Covered | exact values asserted |
| M4 | ✅ Covered | exact string asserted |
| M5 | ✅ Covered | exact string asserted |
| M6 | ✅ Covered | exact string asserted + score-preservation |
| M7 | ✅ Covered | exact string asserted |
| M8 | ✅ Covered | exception asserted |
| M9 | ✅ Covered | exception asserted |
| M10 | ✅ Covered | exception asserted |
| S1 | ✅ Covered | empty array asserted |
| S2 | ✅ Covered | empty string asserted |
| S3 | ✅ Covered | exception asserted |
| S4 | ✅ Covered | exception asserted |
| S5 | ✅ Covered | exception asserted |
| S6 | ✅ Covered | identity asserted |
| C1 | ✅ Covered | round-trip property asserted |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Exception type for malformed/out-of-range decode is `ArgumentOutOfRangeException` (range bounds are source-backed; type is API-shape) | M8, M10 |
| 2 | Exception type for Phred+33→Phred+64 overflow is `ArgumentOutOfRangeException` (non-representability is source-backed; type is API-shape) | M9 |

---

## 7. Open Questions / Decisions

1. `Auto` encoding detection is not exercised by these canonical methods (callers pass an explicit `Phred33`/`Phred64`); auto-detection is an out-of-scope existing method.
