# Test Specification: SEQ-SECSTRUCT-001

**Test Unit ID:** SEQ-SECSTRUCT-001
**Area:** Statistics
**Algorithm:** Protein Secondary Structure Prediction — Chou-Fasman conformational propensities (sliding-window profile)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Chou & Fasman (1978) Empirical predictions of protein conformation, Annu Rev Biochem 47:251-276 | 1 | https://pubmed.ncbi.nlm.nih.gov/354496/ | 2026-06-13 |
| 2 | Chou & Fasman (1974) Prediction of protein conformation, Biochemistry 13:222-245 | 1 | (cited by #3, #8) | 2026-06-13 |
| 3 | Wikipedia — Chou–Fasman method | 4 | https://en.wikipedia.org/wiki/Chou%E2%80%93Fasman_method | 2026-06-13 |
| 4 | Kelley bioinfo — Chou-Fasman algorithm (PDF) | 4 | https://www.kelleybioinfo.org/algorithms/background/BCho.pdf | 2026-06-13 |
| 5 | CSB\|SJU (Jakubowski) — Chou-Fasman propensities | 4 | https://employees.csbsju.edu/hjakubowski/classes/ch331/protstructure/tablechoufas.htm | 2026-06-13 |
| 6 | Przytycka (NCBI) — Protein sec. struct. prediction (PDF) | 4 | https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/download/lectures/CAMS_02_Prot_Sec_Str.pdf | 2026-06-13 |
| 7 | ravihansa3000/ChouFasman reference impl (ChouFasman.py) | 3 | https://raw.githubusercontent.com/ravihansa3000/ChouFasman/master/ChouFasman.py | 2026-06-13 |
| 8 | Chen et al. (2006) Improved Chou-Fasman method, BMC Bioinformatics 7(S4):S14 | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC1780123/ | 2026-06-13 |

### 1.2 Key Evidence Points

1. Each residue has Pα (helix), Pβ (sheet), Pt (turn) propensities; propensity = observed/expected occurrence (integer convention ×100) — sources #4, #6, #7.
2. Verbatim per-residue Pa/Pb/Pt for all 20 residues match the implemented table, e.g. A 1.42/0.83/0.66, E 1.51/0.37/0.74, V 1.06/1.70/0.50, N 0.67/0.89/1.56 — sources #5, #6, #7.
3. Lysine Pα conflict: source #5 = 1.16; sources #6 and #7 = 1.14 → adopt 1.14 (two-source majority + integer-convention consistency) — see §6.
4. The published method scans a 6-residue helix nucleation window (4/6) and a 5-residue sheet window (3/5); the method under test is a generic configurable sliding-window *mean* of these propensities, stepping one residue — sources #3, #4, #8.

### 1.3 Documented Corner Cases

- Parameter table covers exactly the 20 standard residues; X/B/Z/gaps have no defined value and are excluded from window means (sources #5, #7).
- A window longer than the sequence produces no scan positions (window-vs-length, source #4).

### 1.4 Known Failure Modes / Pitfalls

1. Low Q3 accuracy (~50-60%); parameters derived from a small 29-protein sample — source #3. (Accuracy caveat only; the propensity values themselves are the formal Chou-Fasman set.)
2. Using the wrong lysine Pα (1.16 vs 1.14) shifts any window mean containing K — sources #5 vs #6/#7.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `PredictSecondaryStructure(string proteinSequence, int windowSize = 7)` | SequenceStatistics | Canonical | Per-window mean Chou-Fasman (Pa, Pb, Pt) propensity profile, step 1, N→C order. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-01 | A single-residue window (windowSize = 1) returns exactly that residue's (Pa, Pb, Pt) tuple from the Chou-Fasman table. | Yes | Sources #5, #6, #7 |
| INV-02 | Each emitted tuple equals the arithmetic mean of the member residues' propensities (per component). | Yes | Definition of profile + #4 |
| INV-03 | The number of emitted windows = max(0, n − windowSize + 1) for a sequence of length n with all-known residues. | Yes | Window-scan definition (#4) |
| INV-04 | Result is case-insensitive: lower-case input yields identical tuples to upper-case. | Yes | Implementation uppercases input |
| INV-05 | Unknown residues are excluded from a window's count and mean; an all-unknown window emits nothing. | Yes | ASSUMPTION (Evidence §Assumptions 3) |
| INV-06 | Null/empty input, windowSize > n, or windowSize < 1 → empty result (no exception). | Yes | Precondition contract (#4) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Single residue A, window 1 | windowSize=1 over "A" | one tuple (1.42, 0.83, 0.66) | #5,#6,#7 |
| M2 | Single residue E, window 1 | "E" | one tuple (1.51, 0.37, 0.74) | #6,#7 |
| M3 | Single residue V, window 1 | "V" | one tuple (1.06, 1.70, 0.50) | #5,#7 |
| M4 | Lysine K, window 1 | "K" — conflict-resolved | one tuple (1.14, 0.74, 1.01); helix NOT 1.16 | #6,#7 vs #5 |
| M5 | Two-residue mean | "AE", window 2 | one tuple ((1.465), (0.60), (0.70)) | #5,#6,#7 + INV-02 |
| M6 | Three-residue mean | "AEV", window 3 | one tuple (1.33, 2.90/3, 1.90/3) | #5,#6,#7 + INV-02 |
| M7 | Sliding step + count | "AEV", window 2 | two tuples: [ (1.465,0.60,0.70), ((1.51+1.06)/2,(0.37+1.70)/2,(0.74+0.50)/2) ] | INV-02, INV-03 |
| M8 | Case-insensitive | "ae" vs "AE", window 2 | identical tuples | INV-04 |
| M9 | Unknown residue excluded | "AXE", window 3 | one tuple = mean of A,E only = (1.465,0.60,0.70) | INV-05 |
| M10 | All-unknown window | "XBZ", window 3 | empty result | INV-05 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Null input | null, window 7 | empty | INV-06 |
| S2 | Empty input | "", window 7 | empty | INV-06 |
| S3 | Window > length | "AE", window 7 | empty | INV-06 |
| S4 | Non-positive window | "AEV", window 0 | empty | INV-06 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Helix-favouring peptide | "AEMLK" window len | mean Helix > mean Sheet | Exact means from table |
| C2 | Sheet-favouring peptide | "VIY" window len | mean Sheet > mean Helix | Exact means from table |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing test file for `PredictSecondaryStructure` (`grep` over `tests/` found no `PredictSecondaryStructure` references). Production method already existed in `SequenceStatistics.cs` (Present-but-nonconforming: Lysine Pα was 1.16).

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
| M9 | ❌ Missing | new unit |
| M10 | ❌ Missing | new unit |
| S1 | ❌ Missing | new unit |
| S2 | ❌ Missing | new unit |
| S3 | ❌ Missing | new unit |
| S4 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |
| C2 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_PredictSecondaryStructure_Tests.cs` — all cases above.
- **Remove:** none (no prior tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceStatistics_PredictSecondaryStructure_Tests.cs` | Canonical | 16 |

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
| 15 | C1 | ❌ Missing | Implemented | ✅ Done |
| 16 | C2 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 16
**✅ Done:** 16 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | exact tuple |
| M2 | ✅ Covered | exact tuple |
| M3 | ✅ Covered | exact tuple |
| M4 | ✅ Covered | conflict-resolved value asserted |
| M5 | ✅ Covered | exact mean |
| M6 | ✅ Covered | exact mean |
| M7 | ✅ Covered | step + count |
| M8 | ✅ Covered | case-insensitive |
| M9 | ✅ Covered | unknown excluded |
| M10 | ✅ Covered | all-unknown empty |
| S1 | ✅ Covered | null empty |
| S2 | ✅ Covered | empty empty |
| S3 | ✅ Covered | window>len empty |
| S4 | ✅ Covered | window<1 empty |
| C1 | ✅ Covered | helix>sheet |
| C2 | ✅ Covered | sheet>helix |

All 16 in-scope cases ✅.

---

## 6. Assumption Register

**Total assumptions:** 3

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Lysine Pα = 1.14 (not 1.16) — two-source majority (#6, #7) over #5 | M4, implementation table |
| 2 | Default window size 7 is an API convenience, not a Chou-Fasman constant; tests pass window explicitly | §2, §4 |
| 3 | Unknown residues skipped and excluded from window mean; all-unknown window emits nothing | M9, M10, INV-05 |

---

## 7. Open Questions / Decisions

1. The method is a sliding-window propensity *profile*, not the full Chou-Fasman nucleation/extension/turn state machine. The full classifier (helix/sheet/turn segment assignment, 4-of-6 nucleation, p-value turn product) is **Not implemented** here and is documented as out-of-scope in the algorithm doc §5.3.
