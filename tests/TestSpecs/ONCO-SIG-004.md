# Test Specification: ONCO-SIG-004

**Test Unit ID:** ONCO-SIG-004
**Area:** Oncology
**Algorithm:** Mutational Process Classification (SBS exposure → active mutational processes)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Rosenthal et al. (2016) deconstructSigs, *Genome Biology* 17:31 | 1 | https://doi.org/10.1186/s13059-016-0893-4 | 2026-06-14 |
| 2 | deconstructSigs `whichSignatures.R` (reference implementation) | 3 | https://github.com/raerose01/deconstructSigs/blob/master/R/whichSignatures.R | 2026-06-14 |
| 3 | COSMIC Mutational Signatures — SBS aetiologies | 5 | https://cancer.sanger.ac.uk/signatures/sbs/ | 2026-06-14 |
| 4 | Alexandrov et al. (2020) *Nature* 578:94–101 | 1 | https://doi.org/10.1038/s41586-020-1943-3 | 2026-06-14 |

### 1.2 Key Evidence Points

1. Signature activities are reported as **normalized relative contributions** (weights "normalized between 0 and 1"). — Rosenthal 2016.
2. A signature is declared **present/active** only when its normalized contribution ≥ **0.06**; below that it is set to zero. — Rosenthal 2016; `signature.cutoff = 0.06`, `weights[weights < signature.cutoff] <- 0` (deconstructSigs source).
3. Cutoff comparison is **strict less-than** (`weights < signature.cutoff`), so exactly 0.06 is retained. — deconstructSigs source.
4. COSMIC SBS→aetiology map (verbatim): SBS1 "Spontaneous deamination of 5-methylcytosine (clock-like)"; SBS5 "Unknown (clock-like)" → Aging/clock-like. SBS2, SBS13 "Activity of APOBEC family of cytidine deaminases" → APOBEC. SBS4 "Tobacco smoking" → Tobacco. SBS7a–d "Ultraviolet light exposure" → UV. SBS6, SBS15, SBS26 "Defective DNA mismatch repair", SBS20 "Concurrent POLD1 mutations and defective DNA mismatch repair" → MMR deficiency. — COSMIC SBS pages.

### 1.3 Documented Corner Cases

- Sub-cutoff signatures are forced to zero, so surviving contributions can sum to < 1 (remainder is "unknown"). — Rosenthal 2016.
- Multiple processes can be active simultaneously; classification reports the full active set. — Rosenthal 2016.
- Signatures with unknown/unmapped aetiology contribute to no named process. — COSMIC.

### 1.4 Known Failure Modes / Pitfalls

1. Applying the cutoff to per-process totals instead of per-signature contributions changes which signatures survive. — resolved per ASSUMPTION-2 (cutoff per signature, then aggregate).
2. All-zero exposures: Σ = 0, no normalization possible ⇒ no active processes. — degenerate input.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `ClassifyMutationalProcess(exposures, contributionCutoff)` | OncologyAnalyzer | Canonical | Normalizes exposures, applies 6% cutoff, aggregates by COSMIC process, returns active set + dominant process |
| `GetMutationalProcess(signatureLabel)` | OncologyAnalyzer | Canonical | COSMIC SBS label → process lookup (the aetiology map) |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Each surviving signature contribution ∈ [0,1]; Σ surviving contributions ≤ 1 | Yes | Rosenthal 2016 (weights normalized 0–1; sub-cutoff dropped) |
| INV-2 | A signature is active iff normalized contribution ≥ 0.06 (strict `<` excludes) | Yes | deconstructSigs `weights < signature.cutoff`, cutoff 0.06 |
| INV-3 | Per-process contribution = Σ surviving member-signature contributions | Yes | additive weights (Rosenthal 2016) + ASSUMPTION-1 |
| INV-4 | Dominant process = the active process with maximum aggregated contribution; none when no process active | Yes | derivation |
| INV-5 | Σ exposure = 0 ⇒ no active processes, no dominant process | Yes | normalization undefined for zero total |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Normalized contribution | exposures {SBS2:50, SBS13:30, SBS1:15, SBS4:5} → normalized {0.50,0.30,0.15,0.05} | contributions = raw/100 | Rosenthal 2016 (normalized 0–1) |
| M2 | Sub-cutoff excluded | SBS4 at 0.05 (< 0.06) dropped | SBS4 not in any active process | deconstructSigs 0.06 cutoff |
| M3 | Cutoff boundary retained | a signature at exactly 0.06 normalized | retained (strict `<`) | `weights < signature.cutoff` |
| M4 | Just-below boundary excluded | a signature at 0.059… (< 0.06) | excluded | `weights < signature.cutoff` |
| M5 | APOBEC aggregation | SBS2(0.50)+SBS13(0.30) | APOBEC contribution = 0.80 | COSMIC SBS2/13 APOBEC; additive |
| M6 | Active set | dataset above | active = {APOBEC, Aging}; Tobacco excluded | hand-derived dataset |
| M7 | Dominant process | dataset above | dominant = APOBEC (0.80) | INV-4 + dataset |
| M8 | MMR aggregation | SBS6,SBS15,SBS20,SBS26 all map to MMR deficiency | all four → MMR deficiency process | COSMIC aetiologies |
| M9 | Process lookup (each etiology) | GetMutationalProcess for SBS1,SBS5,SBS2,SBS13,SBS4,SBS7a,SBS6,SBS15,SBS20,SBS26 | Aging,Aging,APOBEC,APOBEC,Tobacco,UV,MMR,MMR,MMR,MMR | COSMIC aetiologies |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | All-zero exposures | every exposure 0 | no active processes, no dominant | INV-5 |
| S2 | Unmapped signature | a label not in the COSMIC map | contributes to no process; lookup returns Unknown/null | COSMIC unknown |
| S3 | Single dominant signature | one signature at 1.0 | one active process, it is dominant | trivial dominance |
| S4 | Custom cutoff override | cutoff = 0.20 drops the 0.15 signature | Aging excluded | deconstructSigs cutoff parameter |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | UV multi-subtype | SBS7a + SBS7b both present | both → UV process, summed | COSMIC SBS7a–d all UV |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing tests for `ClassifyMutationalProcess` / `GetMutationalProcess`. Sibling signature tests exist (`OncologyAnalyzer_FitSignatures_Tests.cs`, `OncologyAnalyzer_BootstrapExposures_Tests.cs`, `OncologyAnalyzer_ClassifySbsContext_Tests.cs`) but cover other units. New canonical file `OncologyAnalyzer_ClassifyMutationalProcess_Tests.cs`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M9, S1–S4, C1 | ❌ Missing | brand-new unit; no prior tests |
| Null/empty/invalid input guards | ❌ Missing | new |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ClassifyMutationalProcess_Tests.cs` — all ONCO-SIG-004 cases.
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| OncologyAnalyzer_ClassifyMutationalProcess_Tests.cs | Canonical (ONCO-SIG-004) | 18 |

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
| 9 | M9 | ❌ Missing | implemented (10 labels) | ✅ Done |
| 10 | S1 | ❌ Missing | implemented | ✅ Done |
| 11 | S2 | ❌ Missing | implemented | ✅ Done |
| 12 | S3 | ❌ Missing | implemented | ✅ Done |
| 13 | S4 | ❌ Missing | implemented | ✅ Done |
| 14 | C1 | ❌ Missing | implemented | ✅ Done |
| 15 | Guards | ❌ Missing | null/empty/negative/bad-cutoff implemented | ✅ Done |

**Total items:** 15
**✅ Done:** 15 | **⛔ Blocked:** 0 | **Remaining:** must be 0 → 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | normalized contributions exact |
| M2 | ✅ | sub-cutoff dropped |
| M3 | ✅ | boundary 0.06 retained |
| M4 | ✅ | just-below excluded |
| M5 | ✅ | APOBEC summed = 0.80 |
| M6 | ✅ | active set {APOBEC, Aging} |
| M7 | ✅ | dominant = APOBEC |
| M8 | ✅ | MMR aggregation |
| M9 | ✅ | 10-label process lookup |
| S1 | ✅ | all-zero ⇒ empty |
| S2 | ✅ | unmapped → no process |
| S3 | ✅ | single dominant |
| S4 | ✅ | custom cutoff |
| C1 | ✅ | UV multi-subtype |
| Guards | ✅ | null/empty/negative/bad-cutoff |

**In-scope cases:** 15 — **✅:** 15.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Per-process contribution = sum of surviving member-signature contributions | M5, M8, INV-3 |
| 2 | 6% cutoff applied to per-signature normalized contribution, then aggregate by process | M2, M3, M4, INV-2 |

---

## 7. Open Questions / Decisions

1. None. SBS→process map and 6% cutoff are source-backed; reference signature profiles remain caller-supplied (only labels + exposures are consumed here).
