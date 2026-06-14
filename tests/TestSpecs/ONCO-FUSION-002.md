# Test Specification: ONCO-FUSION-002

**Test Unit ID:** ONCO-FUSION-002
**Area:** Oncology
**Algorithm:** Known Fusion Database Lookup (HGNC gene-fusion designation + caller-supplied known-fusion match)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Bruford et al. (2021). HGNC recommendations for the designation of gene fusions. *Leukemia* 35(11):3040–3043. | 2 | https://pmc.ncbi.nlm.nih.gov/articles/PMC8550944/ | 2026-06-14 |
| 2 | Recommendations for future extensions to the HGNC gene fusion nomenclature. *Leukemia* 35(11):3044–3045. | 2 | https://pmc.ncbi.nlm.nih.gov/articles/PMC8632684/ | 2026-06-14 |

### 1.2 Key Evidence Points

1. The separator between fused gene symbols is a double colon `::`, e.g. `BCR::ABL1` — Source 1.
2. The 5′ partner is always listed first, before the `::`, irrespective of chromosomal location/orientation — Source 1.
3. Order is therefore directional: `A::B` and `B::A` are different fusions (reciprocal) — derived from Source 1, point 2.
4. Genes are designated by HGNC approved gene symbols — Source 1.
5. Read-through transcripts use a hyphen (`INS-IGF2`), not `::` — Source 1.
6. `::` is the cross-resource standard endorsed by WHO/COSMIC/Mitelman/Atlas — Sources 1 & 2.

### 1.3 Documented Corner Cases

- Directional designation: a known-fusion lookup keyed by `5′::3′` must NOT match the reciprocal `3′::5′` (Source 1, point 2).
- Hyphen vs double colon: read-throughs use `-`; this unit emits `::` for true fusions only (Source 1).

### 1.4 Known Failure Modes / Pitfalls

1. Treating the gene pair as unordered (alphabetical/set) would collapse reciprocal fusions — violates the 5′-first rule — Source 1.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `GetFusionAnnotation(gene5p, gene3p)` | OncologyAnalyzer | Canonical | Formats the HGNC `5′::3′` designation string |
| `MatchKnownFusions(fusion, knownFusions)` | OncologyAnalyzer | Canonical | Directional lookup of a `FusionCall` against a caller-supplied known-fusion map |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Designation = `gene5p` + `::` + `gene3p` (5′ first, double colon) | Yes | Source 1, points 1–2 |
| INV-2 | Designation is directional: `Annotation(A,B) ≠ Annotation(B,A)` when A≠B | Yes | Source 1, point 2 |
| INV-3 | A match requires the *directional* key `5′::3′` to be present; the reciprocal key does not match | Yes | Source 1, point 2 |
| INV-4 | Symbol comparison is case-insensitive but order-preserving | Yes | **ASSUMPTION** (Evidence Assumption 2) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Annotation BCR/ABL1 | `GetFusionAnnotation("BCR","ABL1")` | `"BCR::ABL1"` | Source 1, points 1,2,5 |
| M2 | Annotation directional | `GetFusionAnnotation("ABL1","BCR")` | `"ABL1::BCR"` (≠ `BCR::ABL1`) | Source 1, point 2 |
| M3 | Match present | `MatchKnownFusions(EML4→ALK, {"EML4::ALK":"…"})` | matched=true, annotation returned | Source 1, points 1,2 |
| M4 | Match reciprocal absent | known set has only `"ALK::EML4"`, query EML4→ALK | matched=false | Source 1, point 2 (INV-3) |
| M5 | Match absent | known set has `"BCR::ABL1"`, query EML4→ALK | matched=false | Source 1 (directional keying) |
| M6 | Null gene5p annotation | `GetFusionAnnotation(null,"ABL1")` | `ArgumentException` | input-validation contract |
| M7 | Empty gene3p annotation | `GetFusionAnnotation("BCR","")` | `ArgumentException` | input-validation contract |
| M8 | Null known-fusion map | `MatchKnownFusions(fusion, null)` | `ArgumentNullException` | input-validation contract |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Case-insensitive match | known `"EML4::ALK"`, query `eml4`→`alk` | matched=true | INV-4 (Assumption 2) |
| S2 | Annotation preserves input case | `GetFusionAnnotation("bcr","abl1")` | `"bcr::abl1"` | designation is verbatim concatenation; matching is case-insensitive |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Integration with DetectFusions | match a `FusionCall` from `DetectFusions` | matched=true, designation `EML4::ALK` | ties ONCO-FUSION-001 → 002 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Searched `tests/Seqeron/Seqeron.Genomics.Tests/` for `MatchKnownFusions` / `GetFusionAnnotation`: none found.
- ONCO-FUSION-001 tests exist in `OncologyAnalyzer_DetectFusions_Tests.cs` but cover detection only, not lookup.

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
| C1 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_MatchKnownFusions_Tests.cs` — all cases for both methods.
- **Remove:** nothing (no pre-existing tests for these methods).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `OncologyAnalyzer_MatchKnownFusions_Tests.cs` | Canonical (M1–M8, S1–S2, C1) | 11 |

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
| 9 | S1 | ❌ Missing | implemented | ✅ Done |
| 10 | S2 | ❌ Missing | implemented | ✅ Done |
| 11 | C1 | ❌ Missing | implemented | ✅ Done |

**Total items:** 11
**✅ Done:** 11 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `GetFusionAnnotation_BcrAbl1_FivePrimeFirstDoubleColon` |
| M2 | ✅ Covered | `GetFusionAnnotation_Reciprocal_IsDirectional` |
| M3 | ✅ Covered | `MatchKnownFusions_DesignationPresent_ReturnsAnnotation` |
| M4 | ✅ Covered | `MatchKnownFusions_OnlyReciprocalPresent_NoMatch` |
| M5 | ✅ Covered | `MatchKnownFusions_DesignationAbsent_NoMatch` |
| M6 | ✅ Covered | `GetFusionAnnotation_NullFivePrime_Throws` |
| M7 | ✅ Covered | `GetFusionAnnotation_EmptyThreePrime_Throws` |
| M8 | ✅ Covered | `MatchKnownFusions_NullKnownSet_Throws` |
| S1 | ✅ Covered | `MatchKnownFusions_CaseInsensitiveSymbols_Match` |
| S2 | ✅ Covered | `GetFusionAnnotation_PreservesInputCase` |
| C1 | ✅ Covered | `MatchKnownFusions_FromDetectFusions_Match` |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Symbol comparison is case-insensitive (ordinal-ignore-case), order-preserving (Evidence Assumption 2) | INV-4, S1 |

(Known-fusion membership being caller-supplied is a scope decision per the unit mandate, not a correctness-affecting assumption about the algorithm's formal behavior — the format/keying rules are fully source-backed.)

---

## 7. Open Questions / Decisions

1. None. The designation format and directional keying are fully defined by Bruford et al. (2021); the known-fusion set is intentionally caller-supplied (Framework algorithm).
