# Test Specification: ANNOT-CODING-001

**Test Unit ID:** ANNOT-CODING-001
**Area:** Annotation
**Algorithm:** Coding Potential Calculation (CPAT hexamer usage-bias score)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Wang L et al. (2013) CPAT, Nucleic Acids Res 41(6):e74 | 1 | https://doi.org/10.1093/nar/gkt006 (https://pmc.ncbi.nlm.nih.gov/articles/PMC3616698/) | 2026-06-13 |
| 2 | CPAT/lncScore `cpmodule/FrameKmer.py` (`kmer_ratio`) | 3 | https://raw.githubusercontent.com/WGLab/lncScore/master/tools/cpmodule/FrameKmer.py | 2026-06-13 |
| 3 | Fickett & Tung (1992) Assessment of protein coding measures | 1 | https://doi.org/10.1093/nar/20.24.6441 | 2026-06-13 |
| 4 | Fickett (1982) / EMBOSS tcode (TESTCODE, alternative) | 1/3 | https://doi.org/10.1093/nar/10.17.5303 ; https://www.bioinformatics.nl/cgi-bin/emboss/help/tcode | 2026-06-13 |

### 1.2 Key Evidence Points

1. Per-hexamer score = `ln(coding[k] / noncoding[k])` when both frequencies > 0 — Source 2, `kmer_ratio`.
2. Coding-only hexamer (noncoding==0) contributes +1; noncoding-only (coding==0) contributes −1 — Source 2.
3. Sequence score = sum of per-hexamer contributions divided by the number of scored hexamers — Source 2 (`sum_of_log_ratio_0/frame0_count`).
4. Hexamers are extracted in-frame: window starts at 0, step = 3, word size = 6; only full-length words are scored — Source 2, `word_generator`.
5. A hexamer absent from either table is skipped and not counted — Source 2 (`has_key` guard).
6. `len(seq) < word_size` → score 0 — Source 2.
7. Logarithm is natural (base e) — Source 2 (`math.log`).
8. Positive score = coding, negative = non-coding — Source 1 (interpretation).

### 1.3 Documented Corner Cases

- Sequence shorter than word size → 0 (Source 2).
- Hexamer in only one table → skipped (Source 2).
- coding>0 & noncoding==0 → +1; coding==0 & noncoding>0 → −1 (Source 2).
- No scorable hexamer → reference returns −1 (caught division); C# port returns 0 (see §6 ASSUMPTION-1).

### 1.4 Known Failure Modes / Pitfalls

1. Using base-10 instead of natural log changes every score by a constant factor — Source 2 (`math.log` = ln).
2. Counting skipped hexamers in the denominator inflates/deflates the mean — Source 2 (`frame0_count` only incremented for in-both hexamers).
3. Mixing table units (counts in one, proportions in the other) shifts the score by `ln(Σc/Σn)` — Source 2 (`kmer_freq_file`), §6 ASSUMPTION-2.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateCodingPotential(string, IReadOnlyDictionary<string,double>, IReadOnlyDictionary<string,double>, int, int)` | GenomeAnnotator | Canonical | CPAT `kmer_ratio` (frame 0); deep evidence-based testing |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Score = (Σ per-hexamer contributions) / (number of scored in-frame hexamers) | Yes | Source 2 |
| INV-2 | Per-hexamer contribution = `ln(coding[k]/noncoding[k])` when both > 0 | Yes | Source 2 |
| INV-3 | Coding-only hexamer ⇒ +1; noncoding-only ⇒ −1 | Yes | Source 2 |
| INV-4 | Coding-biased tables ⇒ score > 0; noncoding-biased ⇒ score < 0 | Yes | Source 1 |
| INV-5 | Only in-frame full-length hexamers (start 0, step 3, len = wordSize) are scored | Yes | Source 2 |
| INV-6 | `sequence.Length < wordSize` ⇒ score = 0 | Yes | Source 2 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Mean log-ratio, two hexamers | `ATGAAA`/`AAACCC`, coding {8,2}, noncoding {2,4} | 0.34657359027997264 | Source 2 + worked derivation |
| M2 | Single in-both hexamer | one hexamer coding=4, noncoding=1 → ln 4 / 1 | 1.3862943611198906 | Source 2 |
| M3 | Coding-only pseudo-score | one hexamer coding=5, noncoding=0 | +1.0 | Source 2 |
| M4 | Noncoding-only pseudo-score | one hexamer coding=0, noncoding=5 | −1.0 | Source 2 |
| M5 | Sign invariant — coding | coding table dominates | score > 0 (exact 0.34657...) | Source 1, 2 |
| M6 | Sign invariant — noncoding | noncoding table dominates | score < 0 (exact −0.34657...) | Source 1, 2 |
| M7 | Too-short sequence | length < wordSize | 0.0 | Source 2 |
| M8 | Skip hexamer missing from a table | one hexamer absent in noncoding, one in both | mean over the one scored hexamer only | Source 2 |
| M9 | In-frame stepping | hexamers only at offsets 0,3,6,... | out-of-frame hexamers never scored | Source 2 |
| M10 | Null sequence | sequence = null | `ArgumentNullException` | validation contract |
| M11 | Null coding table | coding = null | `ArgumentNullException` | validation contract |
| M12 | Null noncoding table | noncoding = null | `ArgumentNullException` | validation contract |
| M13 | Non-positive wordSize | wordSize = 0 | `ArgumentOutOfRangeException` | validation contract |
| M14 | Non-positive stepSize | stepSize = 0 | `ArgumentOutOfRangeException` | validation contract |
| M15 | No scorable hexamer | all hexamers missing from tables | 0.0 (see ASSUMPTION-1) | Source 2 + §6 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Case-insensitivity | lowercase sequence | same score as uppercase | `ToUpperInvariant` |
| S2 | Empty sequence | "" | 0.0 | length 0 < wordSize |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Coding==0 & noncoding==0 entry | hexamer present in both tables with both values 0 | skipped (`continue`), NOT counted | matches reference branch `elif coding[k]==0 and noncoding[k]==0: continue` |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No prior test file for `GenomeAnnotator.CalculateCodingPotential`. Searched `tests/Seqeron/Seqeron.Genomics.Tests/` (no `*CodingPotential*` file). The prior production code used an invented heuristic (coefficients 0.7/0.6) and had no tests.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M15, S1–S2, C1 | ❌ Missing | New unit; no prior tests exist |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/GenomeAnnotator_CalculateCodingPotential_Tests.cs` — all cases for the CPAT hexamer score.
- **Remove:** nothing (no prior tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| GenomeAnnotator_CalculateCodingPotential_Tests.cs | Canonical | 18 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented exact mean-log-ratio test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented single-hexamer ln test | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented coding-only +1 test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented noncoding-only −1 test | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented coding-sign test | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented noncoding-sign test | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented too-short test | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented skip-missing test | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented in-frame-stepping test | ✅ Done |
| 10 | M10 | ❌ Missing | Implemented null-sequence test | ✅ Done |
| 11 | M11 | ❌ Missing | Implemented null-coding test | ✅ Done |
| 12 | M12 | ❌ Missing | Implemented null-noncoding test | ✅ Done |
| 13 | M13 | ❌ Missing | Implemented wordSize-0 test | ✅ Done |
| 14 | M14 | ❌ Missing | Implemented stepSize-0 test | ✅ Done |
| 15 | M15 | ❌ Missing | Implemented no-scorable-hexamer test | ✅ Done |
| 16 | S1 | ❌ Missing | Implemented case-insensitivity test | ✅ Done |
| 17 | S2 | ❌ Missing | Implemented empty-sequence test | ✅ Done |
| 18 | C1 | ❌ Missing | Implemented both-zero-skipped test (corrected from counted in ANNOT-CODING-001 validation) | ✅ Done |

**Total items:** 18
**✅ Done:** 18 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | CalculateCodingPotential_TwoInFrameHexamers_ReturnsMeanLogRatio |
| M2 | ✅ Covered | CalculateCodingPotential_SingleHexamerBothTables_ReturnsNaturalLogRatio |
| M3 | ✅ Covered | CalculateCodingPotential_CodingOnlyHexamer_ContributesPlusOne |
| M4 | ✅ Covered | CalculateCodingPotential_NoncodingOnlyHexamer_ContributesMinusOne |
| M5 | ✅ Covered | CalculateCodingPotential_CodingBiasedTables_ReturnsPositiveScore |
| M6 | ✅ Covered | CalculateCodingPotential_NoncodingBiasedTables_ReturnsNegativeScore |
| M7 | ✅ Covered | CalculateCodingPotential_SequenceShorterThanWord_ReturnsZero |
| M8 | ✅ Covered | CalculateCodingPotential_HexamerMissingFromTable_IsSkipped |
| M9 | ✅ Covered | CalculateCodingPotential_OutOfFrameHexamer_IsNotScored |
| M10 | ✅ Covered | CalculateCodingPotential_NullSequence_Throws |
| M11 | ✅ Covered | CalculateCodingPotential_NullCodingTable_Throws |
| M12 | ✅ Covered | CalculateCodingPotential_NullNoncodingTable_Throws |
| M13 | ✅ Covered | CalculateCodingPotential_NonPositiveWordSize_Throws |
| M14 | ✅ Covered | CalculateCodingPotential_NonPositiveStepSize_Throws |
| M15 | ✅ Covered | CalculateCodingPotential_NoScorableHexamer_ReturnsZero |
| S1 | ✅ Covered | CalculateCodingPotential_LowercaseSequence_MatchesUppercase |
| S2 | ✅ Covered | CalculateCodingPotential_EmptySequence_ReturnsZero |
| C1 | ✅ Covered | CalculateCodingPotential_BothTablesZero_HexamerIsSkipped |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| ASSUMPTION-1 | No-scorable-hexamer returns 0 (port choice) instead of the reference −1 sentinel; only affects inputs with zero scorable hexamers | M15, §1.3 |
| ASSUMPTION-2 | Both frequency tables must use the same units (counts or proportions); the score is unit-invariant only within consistent units | Contract, Evidence dataset |

---

## 7. Open Questions / Decisions

1. **Decision:** The Registry title "hexamer frequency bias" maps to the CPAT hexamer usage-bias log-likelihood (Wang et al. 2013), not the Fickett TESTCODE composite. TESTCODE is recorded as a related, not-implemented alternative.
2. **Decision:** The prior single-argument `CalculateCodingPotential(string)` used invented weights (0.7, 0.6) with no authoritative basis; per the implementation conformance policy it was a defect and was replaced by the table-based CPAT method (correctness-affecting constants removed). The MCP `coding_potential` tool now accepts the two hexamer tables.
3. **Decision:** Hexamer tables are algorithm **inputs** (training data), so they are passed as parameters rather than hard-coded; this matches CPAT, where tables are organism-specific and supplied separately.
