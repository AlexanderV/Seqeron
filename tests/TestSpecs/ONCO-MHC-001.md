# Test Specification: ONCO-MHC-001

**Test Unit ID:** ONCO-MHC-001
**Area:** Oncology
**Algorithm:** MHC-Peptide Binding Classification (length filtering + affinity/%rank thresholds)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Reynisson et al. (2020) NetMHCpan-4.1, *Nucleic Acids Res* 48(W1):W449–W454 | 1 | https://doi.org/10.1093/nar/gkaa379 (PMC7319546) | 2026-06-14 |
| 2 | Sette et al. (1994), *J Immunol* 153(12):5586–92 | 1 | https://pubmed.ncbi.nlm.nih.gov/7527444/ | 2026-06-14 |
| 3 | Roomp, Antes & Lengauer (2010), *BMC Bioinformatics* 11:90 | 1 | https://doi.org/10.1186/1471-2105-11-90 (PMC2836306) | 2026-06-14 |
| 4 | IEDB threshold help article | 5 | https://help.iedb.org/hc/en-us/articles/114094152371 | 2026-06-14 |
| 5 | IEDB MHC class II tool description | 5 | https://help.iedb.org/hc/en-us/articles/114094151731 | 2026-06-14 |

### 1.2 Key Evidence Points

1. Class I %Rank: strong binder < 0.5%, weak binder < 2% (default) — Reynisson 2020 (PMC7319546).
2. Class II %Rank: strong binder < 2%, weak binder < 10% (default) — Reynisson 2020.
3. IC50 tiers: high (strong) < 50 nM, intermediate (weak) < 500 nM, low < 5000 nM — IEDB (#4); the 50/500 nM cutoffs trace to Sette 1994 "≈500 nM (preferably 50 nM or less)".
4. 500 nM is the binder/non-binder demarcation — Roomp 2010 (PMC2836306), corroborating IEDB.
5. Class I peptide length: 8–14, default 8–11 — Reynisson 2020. Class II: 13–25 — IEDB (#5).
6. The prediction of the IC50/%Rank value itself requires a trained model (NetMHCpan) and is **caller-supplied / out of scope** — only classification of a supplied value is in scope.

### 1.3 Documented Corner Cases

- Boundary semantics are strict `<`: IC50 = 50 nM is NOT strong (→ weak); IC50 = 500 nM is NOT a binder (→ non-binder). %Rank = 0.5 is NOT strong (→ weak); %Rank = 2.0 is NOT a weak binder (→ non-binder). Source: verbatim "<" inequalities in Reynisson 2020 and IEDB tiers.
- Peptide length outside the class range is not a valid binder candidate (Reynisson 2020; IEDB #5).

### 1.4 Known Failure Modes / Pitfalls

1. IC50 must be a positive concentration (invariant IC50 > 0); zero/negative/NaN/∞ are invalid — Registry invariant; concentration semantics.
2. %Rank is a percentile in [0, 100]; values < 0, > 100, or NaN are invalid — percentile definition (Reynisson 2020 §%Rank).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `ClassifyBindingAffinity(double ic50Nm)` | OncologyAnalyzer | Canonical | IC50 → Strong/Weak/NonBinder (50/500 nM) |
| `ClassifyBindingRank(double percentRank, MhcClass mhcClass)` | OncologyAnalyzer | Canonical | %Rank → Strong/Weak/NonBinder (class I 0.5/2; class II 2/10) |
| `IsValidPeptideLength(int length, MhcClass mhcClass)` | OncologyAnalyzer | Canonical | class I 8–11, class II 13–25 |
| `ClassifyMhcBinding(int peptideLength, double ic50Nm, MhcClass mhcClass)` | OncologyAnalyzer | Delegate | length gate + `ClassifyBindingAffinity` |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | IC50 input must be > 0 (positive concentration); else exception | Yes | Registry invariant "IC50 > 0"; concentration semantics |
| INV-2 | %Rank input must be in [0, 100]; else exception | Yes | Reynisson 2020 (%Rank is a percentile) |
| INV-3 | Strong ⇒ Weak ⇒ NonBinder are mutually exclusive, monotone in the score (smaller IC50/%Rank ⇒ stronger or equal class) | Yes | Reynisson 2020; IEDB tiers |
| INV-4 | Classification cutoffs are strict `<` (boundary value excluded from the stronger class) | Yes | Verbatim "<" in Reynisson 2020 and IEDB tiers |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | IC50 strong | ic50 = 10 nM | Strong | IEDB <50; Sette 1994 |
| M2 | IC50 strong boundary | ic50 = 50 nM (strict `<`) | Weak | IEDB "<50"; INV-4 |
| M3 | IC50 weak | ic50 = 200 nM | Weak | IEDB <500 |
| M4 | IC50 weak boundary | ic50 = 500 nM (strict `<`) | NonBinder | IEDB "<500"; Roomp 2010 |
| M5 | IC50 non-binder | ic50 = 1000 nM | NonBinder | IEDB tiers |
| M6 | %Rank class I strong | rank = 0.4, class I | Strong | Reynisson 2020 (<0.5) |
| M7 | %Rank class I strong boundary | rank = 0.5, class I | Weak | Reynisson 2020 ("<0.5"); INV-4 |
| M8 | %Rank class I weak | rank = 1.0, class I | Weak | Reynisson 2020 (<2) |
| M9 | %Rank class I weak boundary | rank = 2.0, class I | NonBinder | Reynisson 2020 ("<2"); INV-4 |
| M10 | %Rank class I non-binder | rank = 5.0, class I | NonBinder | Reynisson 2020 |
| M11 | %Rank class II strong | rank = 1.5, class II | Strong | Reynisson 2020 (<2) |
| M12 | %Rank class II weak boundary | rank = 10.0, class II | NonBinder | Reynisson 2020 ("<10"); INV-4 |
| M13 | %Rank class II weak | rank = 5.0, class II | Weak | Reynisson 2020 (<10) |
| M14 | Length class I valid | len = 9, class I | true | Reynisson 2020 (8–11) |
| M15 | Length class I too short | len = 7, class I | false | Reynisson 2020 |
| M16 | Length class I above range | len = 12, class I | false | Reynisson 2020 (default 8–11) |
| M17 | Length class II valid | len = 15, class II | true | IEDB #5 (13–25) |
| M18 | Length class II too short | len = 12, class II | false | IEDB #5 |
| M19 | Length class II too long | len = 26, class II | false | IEDB #5 |
| M20 | Combined gate fails on length | len = 7, ic50 = 10, class I | NonBinder (invalid length ⇒ not a candidate) | Reynisson 2020 length gate + IEDB affinity |
| M21 | Combined gate passes length, strong affinity | len = 9, ic50 = 10, class I | Strong | combined of M1 + M14 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | IC50 ≤ 0 rejected | ic50 = 0 / −1 | ArgumentOutOfRangeException | INV-1 |
| S2 | IC50 non-finite rejected | ic50 = NaN / ∞ | ArgumentOutOfRangeException | INV-1 |
| S3 | %Rank out of [0,100] rejected | rank = −0.1 / 100.1 | ArgumentOutOfRangeException | INV-2 |
| S4 | %Rank NaN rejected | rank = NaN | ArgumentOutOfRangeException | INV-2 |
| S5 | Length ≤ 0 | len = 0 / −1 | false (not valid) | length is a count |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Monotonicity %Rank class I | rank 0.1 < 0.5 ≤ 1 < 2 ≤ 3 ⇒ classes non-increasing in strength | Strong, Weak, Weak, NonBinder, NonBinder | INV-3 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- New unit. No existing `ClassifyBindingAffinity` / `ClassifyBindingRank` / `IsValidPeptideLength` methods or tests in `OncologyAnalyzer` (grep over `src/.../Seqeron.Genomics.Oncology/OncologyAnalyzer.cs` and `tests/.../OncologyAnalyzer_*` showed only ONCO-NEO-001 windowing). All planned cases start ❌ Missing.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M21 | ❌ Missing | new unit |
| S1–S5 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ClassifyMhcBinding_Tests.cs` — all cases for this unit.
- **Remove:** none (new unit).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| OncologyAnalyzer_ClassifyMhcBinding_Tests.cs | canonical | 27 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1–M5 (IC50) | ❌ Missing | implemented | ✅ Done |
| 2 | M6–M10 (%Rank I) | ❌ Missing | implemented | ✅ Done |
| 3 | M11–M13 (%Rank II) | ❌ Missing | implemented | ✅ Done |
| 4 | M14–M19 (length) | ❌ Missing | implemented | ✅ Done |
| 5 | M20–M21 (combined) | ❌ Missing | implemented | ✅ Done |
| 6 | S1–S5 (validation) | ❌ Missing | implemented | ✅ Done |
| 7 | C1 (monotonicity) | ❌ Missing | implemented | ✅ Done |

**Total items:** 7 groups (27 tests)
**✅ Done:** 7 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1–M21 | ✅ Covered | exact expected values from Evidence |
| S1–S5 | ✅ Covered | exception-type assertions |
| C1 | ✅ Covered | monotonicity property |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Class I accepted length range = 8–11 (Reynisson default; matches ONCO-NEO-001 constants) rather than the full 8–14 | M16, `IsValidPeptideLength` class I |

---

## 7. Open Questions / Decisions

1. The peptide–MHC affinity / %Rank PREDICTION (trained NetMHCpan model, PSSM weights) is **out of scope** for this unit and is caller-supplied input. This unit classifies a supplied IC50 or %Rank; it does not predict it. No model weights are fabricated.
