# Test Specification: ONCO-MSI-001

**Test Unit ID:** ONCO-MSI-001
**Area:** Oncology
**Algorithm:** Microsatellite Instability (MSI) detection — fraction-of-unstable-loci score and status classification
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Niu et al. (2014) MSIsensor, Bioinformatics 30(7):1015–1016 | 1 | https://doi.org/10.1093/bioinformatics/btt755 | 2026-06-14 |
| 2 | niu-lab/msisensor2 README (reference implementation) | 3 | https://github.com/niu-lab/msisensor2 | 2026-06-14 |
| 3 | Boland et al. (1998) NCI Workshop, Cancer Res 58:5248–5257 | 1 | https://pubmed.ncbi.nlm.nih.gov/9823339/ | 2026-06-14 |

### 1.2 Key Evidence Points

1. MSI score = number of unstable (msi) loci / all valid loci, as a fraction/percentage — MSIsensor2 README; Niu 2014 ("percentage of microsatellite sites with a somatic indel").
2. MSI-High when msi score ≥ 20% (inclusive boundary) — MSIsensor2 README ("msi high: msi score >= 20%").
3. Bethesda 5-marker categorical: 0 unstable → MSS, exactly 1 → MSI-L, ≥ 2 → MSI-H — Boland 1998.
4. Per-locus call uses chi-square tumor-vs-normal length distribution; the unit tests the fraction/classification layer given per-locus stability flags (the statistical per-locus test is upstream and out of scope) — Niu 2014.

### 1.3 Documented Corner Cases

- Zero valid loci → score undefined (division by zero) — MSIsensor2 formula.
- MSS vs MSI-L is unreliable on small panels — Boland 1998.
- Tumor-only mode keeps the same fraction definition — MSIsensor2.

### 1.4 Known Failure Modes / Pitfalls

1. Treating the 20% cutoff as exclusive — it is inclusive (≥ 20%) — MSIsensor2 README.
2. Counting unstable > valid loci — invalid; 0 ≤ unstable ≤ valid — MSIsensor2 formula.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateMSIScore(int unstableLoci, int totalLoci)` | OncologyAnalyzer | Canonical | Fraction unstable/total in [0,1]; throws on invalid counts |
| `ClassifyMSIStatus(double score)` | OncologyAnalyzer | Canonical | MSI-H if score ≥ 0.20 (MSIsensor2), else MSS |
| `ClassifyBethesdaPanel(int unstableMarkers, int totalMarkers)` | OncologyAnalyzer | Canonical | 0→MSS, 1→MSI-L, ≥2→MSI-H (Boland 1998) |
| `DetectMSI(IEnumerable<bool> locusUnstableFlags)` | OncologyAnalyzer | Canonical | End-to-end: counts unstable flags, computes score, classifies |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | 0 ≤ MSI score ≤ 1 | Yes | MSIsensor2 (fraction unstable/valid) |
| INV-2 | score == unstableLoci / totalLoci exactly | Yes | MSIsensor2 README |
| INV-3 | ClassifyMSIStatus returns High iff score ≥ 0.20 | Yes | MSIsensor2 README |
| INV-4 | Bethesda: count 0→MSS, 1→MSI-L, ≥2→MSI-H | Yes | Boland 1998 |
| INV-5 | DetectMSI score == CalculateMSIScore(#unstable, #total) | Yes | composition of evidence |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Score basic | 5 unstable / 25 valid | 0.20 | MSIsensor2 (5/25=20%) |
| M2 | Score fraction | 3 unstable / 12 valid | 0.25 | MSIsensor2 formula |
| M3 | Score zero | 0 unstable / 25 valid | 0.0 | MSIsensor2 formula |
| M4 | Score all | 25 unstable / 25 valid | 1.0 | MSIsensor2 formula (INV-1) |
| M5 | Status at boundary | score 0.20 | MSI-H | MSIsensor2 "≥20%" inclusive |
| M6 | Status below boundary | score 0.16 (4/25) | MSS (not High) | MSIsensor2 (<20%) |
| M7 | Status well above | score 0.40 | MSI-H | MSIsensor2 (≥20%) |
| M8 | Status zero | score 0.0 | MSS | MSIsensor2 (<20%) |
| M9 | Bethesda MSS | 0 of 5 unstable | MSS | Boland 1998 |
| M10 | Bethesda MSI-L | 1 of 5 unstable | MSI-L | Boland 1998 |
| M11 | Bethesda MSI-H | 2 of 5 unstable | MSI-H | Boland 1998 |
| M12 | Bethesda MSI-H all | 5 of 5 unstable | MSI-H | Boland 1998 |
| M13 | DetectMSI end-to-end | flags with 6 unstable / 20 | score 0.30, MSI-H | composition (6/20=30% ≥20%) |
| M14 | DetectMSI MSS | flags with 2 unstable / 20 | score 0.10, MSS | composition (10% <20%) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Zero valid loci | CalculateMSIScore(0,0) | throws ArgumentOutOfRangeException | division by zero |
| S2 | unstable > valid | CalculateMSIScore(6,5) | throws ArgumentOutOfRangeException | 0 ≤ unstable ≤ valid |
| S3 | Negative counts | CalculateMSIScore(-1,5) | throws ArgumentOutOfRangeException | counts ≥ 0 |
| S4 | Status invalid score | ClassifyMSIStatus(1.5) / (-0.1) / NaN | throws ArgumentOutOfRangeException | score in [0,1] |
| S5 | Bethesda invalid | ClassifyBethesdaPanel(3,2), (-1,5), (1,0) | throws ArgumentOutOfRangeException | 0 ≤ unstable ≤ total, total>0 |
| S6 | DetectMSI null | DetectMSI(null) | throws ArgumentNullException | guard |
| S7 | DetectMSI empty | DetectMSI(empty) | throws ArgumentOutOfRangeException | no valid loci |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Larger panel 40% | ClassifyBethesdaPanel(4,10) | MSI-H | ≥2 unstable → MSI-H (Boland) |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing MSI methods in `OncologyAnalyzer.cs` (grep for MSI/Microsatellite/DetectMSI returned nothing). No `OncologyAnalyzer_DetectMSI_Tests.cs`. Net-new unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M14, S1–S7, C1 | ❌ Missing | net-new unit, no prior tests |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_DetectMSI_Tests.cs` — all MSI score/status/Bethesda/DetectMSI tests.
- **Remove:** (none — net-new)

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| OncologyAnalyzer_DetectMSI_Tests.cs | canonical | 22 |

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
| 14 | M14 | ❌ Missing | implemented | ✅ Done |
| 15 | S1 | ❌ Missing | implemented | ✅ Done |
| 16 | S2 | ❌ Missing | implemented | ✅ Done |
| 17 | S3 | ❌ Missing | implemented | ✅ Done |
| 18 | S4 | ❌ Missing | implemented | ✅ Done |
| 19 | S5 | ❌ Missing | implemented | ✅ Done |
| 20 | S6 | ❌ Missing | implemented | ✅ Done |
| 21 | S7 | ❌ Missing | implemented | ✅ Done |
| 22 | C1 | ❌ Missing | implemented | ✅ Done |

**Total items:** 22
**✅ Done:** 22 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1–M14 | ✅ Covered | exact evidence values, executed and passing |
| S1–S7 | ✅ Covered | failure modes / guards |
| C1 | ✅ Covered | larger-panel Bethesda cross-check |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | No MSI-L band on the continuous fraction score (MSIsensor2 defines only binary MSI-H ≥20%); MSI-L applies only to discrete Bethesda marker counts | ClassifyMSIStatus (binary), ClassifyBethesdaPanel (three-tier) |

---

## 7. Open Questions / Decisions

1. Decision: per-locus chi-square instability calling is upstream/out of scope; this unit tests the fraction-and-classification layer given per-locus stability flags, matching the checklist canonical signatures (`CalculateMSIScore`, `ClassifyMSIStatus`, `DetectMSI`).
2. Decision: continuous status uses the MSIsensor2 20% cutoff (rank-3 reference impl, the only retrievable continuous-score cutoff). The Bethesda categorical method uses the rank-1 Boland 1998 marker-count rule.
