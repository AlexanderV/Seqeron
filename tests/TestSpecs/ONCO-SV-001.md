# Test Specification: ONCO-SV-001

**Test Unit ID:** ONCO-SV-001
**Area:** Oncology
**Algorithm:** Somatic Complex Rearrangement Classification (Chromothripsis Inference)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-15

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Korbel JO, Campbell PJ (2013). Criteria for Inference of Chromothripsis in Cancer Genomes. Cell 152:1226–1236. | 1 | https://doi.org/10.1016/j.cell.2013.02.023 (PMID 23498933) | 2026-06-15 |
| 2 | Cortés-Ciriano I et al. (2020). Comprehensive analysis of chromothripsis in 2,658 human cancers. Nat Genet 52:331–341. | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC7058534/ | 2026-06-15 |
| 3 | Magrangeas F et al. (2011). Chromothripsis in multiple myeloma. Blood 118:675–678 (≥10 oscillating-CN-change first-pass screen). | 1 | https://doi.org/10.1182/blood-2011-03-344069 | 2026-06-15 |
| 4 | Maher CA, Wilson RK (2012). "Chromothripsis and beyond" review enumerating the Korbel & Campbell criteria. | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC3861665/ | 2026-06-15 |

### 1.2 Key Evidence Points

1. Six hallmark criteria for chromothripsis: clustering of breakpoints; oscillating CN states; interspersed LOH; haplotype-specific rearrangement; random fragment order/joins; ability to walk the derivative chromosome — Source 1 (via Source 4).
2. Canonical CN profile oscillates between **two** copy-number states (alternation of loss and retention), not progressively rising amplification — Source 1/4.
3. First-pass screen threshold: **≥ 10 oscillating copy-number changes** — Magrangeas et al. 2011 (Source 3, via Source 4).
4. High-confidence chromothripsis = oscillation between two states in **≥ 7 adjacent segments**; low-confidence = **4–6 segments** — Source 2.
5. Canonical events have **> 60%** of CN segments oscillating between two states — Source 2.
6. Focal events with **< 6 SVs** are excluded — minimum clustered intrachromosomal SV burden = 6 — Source 2.
7. Breakpoint clustering null: random breakpoints give **exponentially distributed** inter-breakpoint distances (CV = 1 for the exponential) — Source 1/4.

### 1.3 Documented Corner Cases

- Clustering of breakpoints is "necessary but not sufficient": clustered breakpoints with progressively rising CN (>3 states, monotone) are NOT chromothripsis (Source 1/4).
- Fewer than the screening minimum of oscillating CN changes / adjacent oscillating segments → not called (Sources 2, 3).

### 1.4 Known Failure Modes / Pitfalls

1. Confusing progressive amplification (BFB, many ascending CN states) with two-state oscillation — Source 1/4.
2. Calling on clustering alone without checking CN oscillation — Source 1/4.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CountCopyNumberStateOscillations(states)` | OncologyAnalyzer | Canonical | Counts per-segment CN state transitions (the "oscillating CN changes"). |
| `TestBreakpointClustering(positions)` | OncologyAnalyzer | Canonical | Exponential-null clustering summary (inter-breakpoint distances, CV). |
| `ClassifyComplexRearrangement(input)` | OncologyAnalyzer | Canonical | Applies Korbel & Campbell + Cortés-Ciriano criteria to classify chromothripsis vs not, with confidence tier. |
| `ComplexRearrangementResult` / `ComplexRearrangementType` / `ChromothripsisConfidence` | OncologyAnalyzer | Internal | Output records/enums (tested via canonical methods). |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Oscillation count = number of adjacent segments differing in CN state from predecessor; for n segments, 0 ≤ count ≤ n−1. | Yes | Magrangeas 2011 / Korbel & Campbell 2013 |
| INV-2 | A profile is classified Chromothripsis only if it oscillates between ≤ 3 distinct CN states (canonically 2) AND has ≥ 10 oscillations AND SV burden ≥ 6. | Yes | Korbel & Campbell 2013; Cortés-Ciriano 2020 |
| INV-3 | Confidence is High when adjacent oscillating segments ≥ 7, Low when 4–6, None when < 4. | Yes | Cortés-Ciriano 2020 |
| INV-4 | Monotone (strictly rising/falling, > 3 distinct states) CN profile is never Chromothripsis regardless of segment count. | Yes | Korbel & Campbell 2013 (two-state hallmark) |
| INV-5 | Exponential-null clustering: regular-spacing breakpoints give CV ≈ 0 (not clustered); over-dispersed (many short + few long gaps) give CV > 1 (clustered). | Yes | Korbel & Campbell 2013 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Oscillation count on 2,1,2,1,... | 11-segment two-state profile (2,1,…,2) | count = 10 | Magrangeas 2011 (transition count) |
| M2 | Oscillation count, all-equal | [2,2,2,2] | count = 0 | INV-1 |
| M3 | Oscillation count, monotone | [2,3,4,5,6] | count = 4 (each step is a transition) | INV-1 |
| M4 | Classify two-state, 10 osc, SV≥6 | states 2,1×11; 12 SVs | Chromothripsis, 2 states | Korbel & Campbell + Magrangeas (10 screen) |
| M5 | Classify two-state, 5 osc | states 2,1,2,1,2,1 (5 osc); 12 SVs | NotComplex (5 < 10) | Magrangeas 10 screen |
| M6 | Classify monotone rising | 2,3,4,5,6,7,8,9,10,11,12 (10 osc but 11 states) | NotComplex (>3 states, not 2-state oscillation) | Korbel & Campbell two-state hallmark (INV-4) |
| M7 | Classify, SV burden < 6 | valid 10-osc two-state profile but 5 SVs | NotComplex (focal exclusion <6 SVs) | Cortés-Ciriano 2020 |
| M8 | Confidence High | two-state oscillation across ≥7 adjacent segments (≥10 osc) | Confidence = High | Cortés-Ciriano 2020 |
| M9 | Clustering CV on regular spacing | breakpoints 100,200,300,400,500 (equal gaps) | clustered = false, CV ≈ 0 | Korbel & Campbell exponential null |
| M10 | Clustering CV on over-dispersed | tight cluster + far outlier (e.g. 100,101,102,103,5000) | clustered = true, CV > 1 | Korbel & Campbell exponential null |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Confidence Low | Two-state profile with 5 oscillations → 6 oscillating segments (in [4,6]); below the ≥10 screen | Confidence = Low AND Type = NotComplex | Cortés-Ciriano 2020 (tier); Magrangeas 2011 (screen) |
| S2 | Three-state oscillation accepted | CN bounded to {1,2,3} oscillating, 10 osc, ≥6 SV | Chromothripsis (≤3 states allowed) | Korbel & Campbell "two or three" tolerance |
| S3 | Clustering CV exact value | known distances giving an exact CV | CV matches computed value within 1e-10 | exponential CV definition |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Null states | CountCopyNumberStateOscillations(null) | ArgumentNullException | API robustness |
| C2 | Empty / single segment | [] or [2] | count = 0 | INV-1 boundary |
| C3 | Null/empty breakpoints | TestBreakpointClustering(null/empty/one) | ArgumentNullException / not clustered | API robustness |
| C4 | Null classify input | ClassifyComplexRearrangement with null states | ArgumentNullException | API robustness |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing tests for somatic complex-rearrangement / chromothripsis classification in `OncologyAnalyzer`. Generic SV detection lives in `StructuralVariantAnalyzer` (Annotation) and is out of scope. New canonical file: `OncologyAnalyzer_ClassifyComplexRearrangement_Tests.cs`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M10, S1–S3, C1–C4 | ❌ Missing | Brand-new unit; all cases start Missing. |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ClassifyComplexRearrangement_Tests.cs` — all ONCO-SV-001 cases.
- **Remove:** none (new unit).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| OncologyAnalyzer_ClassifyComplexRearrangement_Tests.cs | Canonical | 17 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented oscillation-count test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented all-equal test | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented monotone-transition test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented two-state classify test | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented below-screen test | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented monotone-not-chromothripsis test | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented SV-burden exclusion test | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented high-confidence test | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented regular-spacing clustering test | ✅ Done |
| 10 | M10 | ❌ Missing | Implemented over-dispersed clustering test | ✅ Done |
| 11 | S1 | ❌ Missing | Implemented low-confidence tier test | ✅ Done |
| 12 | S2 | ❌ Missing | Implemented three-state accepted test | ✅ Done |
| 13 | S3 | ❌ Missing | Implemented exact-CV test | ✅ Done |
| 14 | C1 | ❌ Missing | Implemented null-states test | ✅ Done |
| 15 | C2 | ❌ Missing | Implemented empty/single test | ✅ Done |
| 16 | C3 | ❌ Missing | Implemented null/empty breakpoints test | ✅ Done |
| 17 | C4 | ❌ Missing | Implemented null classify input test | ✅ Done |

**Total items:** 17
**✅ Done:** 17 | **⛔ Blocked:** 0 | **Remaining:** must be 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | oscillation count = 10 verified |
| M2 | ✅ | all-equal → 0 |
| M3 | ✅ | monotone transitions counted |
| M4 | ✅ | Chromothripsis, 2 states |
| M5 | ✅ | NotComplex below screen |
| M6 | ✅ | monotone → NotComplex |
| M7 | ✅ | SV<6 → NotComplex |
| M8 | ✅ | High confidence |
| M9 | ✅ | regular spacing not clustered |
| M10 | ✅ | over-dispersed clustered |
| S1 | ✅ | Low confidence tier |
| S2 | ✅ | three-state accepted |
| S3 | ✅ | exact CV |
| C1 | ✅ | null states throws |
| C2 | ✅ | empty/single → 0 |
| C3 | ✅ | null/empty breakpoints handled |
| C4 | ✅ | null classify input throws |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Oscillation = adjacent-segment CN-state transition count (bounded to ≤3 states for chromothripsis). | INV-1, INV-2, M1–M6 |
| 2 | Breakpoint clustering summarised via coefficient of variation of inter-breakpoint distances vs exponential null (CV=1). | INV-5, M9, M10, S3 |

---

## 7. Open Questions / Decisions

1. Exact goodness-of-fit statistic for breakpoint clustering is not fixed by Korbel & Campbell (they only state the exponential null); we expose the CV summary as a transparent, source-anchored proxy rather than a clinical caller. Documented as Assumption 2.
2. Generic SV typing (DEL/DUP/INV/TRA) is deliberately out of scope (covered by `StructuralVariantAnalyzer`); ONCO-SV-001 implements only the oncology-specific complex-rearrangement (chromothripsis) classification layer.
