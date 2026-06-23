# Test Specification: PRIMER-TM-001 (Primer3 weighted penalty objective)

**Test Unit ID:** PRIMER-TM-001
**Area:** MolTools
**Algorithm:** Primer3 weighted per-primer penalty (objective function)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-23

<!-- Companion to tests/TestSpecs/PRIMER-TM-001.md (which covers the Tm calculation,
     validated under SEQ-THERMO-001). This spec covers the Primer3 penalty objective
     added to resolve LIMITATIONS.md §1 PRIMER-TM-001. -->

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Untergasser et al. (2012), Primer3—new capabilities and interfaces, NAR 40(15):e115 | 1 | https://doi.org/10.1093/nar/gks596 | 2026-06-23 |
| 2 | Koressaar & Remm (2007), Bioinformatics 23(10):1289 | 1 | https://doi.org/10.1093/bioinformatics/btm091 | 2026-06-23 |
| 3 | Primer3 source `libprimer3.cc` (`p_obj_fn`, default weights/optima) + `libprimer3.h` | 3 | https://github.com/primer3-org/primer3 | 2026-06-23 |
| 4 | Primer3 manual §19 "HOW PRIMER3 CALCULATES THE PENALTY VALUE" | 2 | https://primer3.org/manual.html | 2026-06-23 |

### 1.2 Key Evidence Points

1. Per-primer penalty = weighted sum of one-sided deviations of Tm, size, GC% from their optima, plus weighted self_any, self_end and num_ns terms; lower is better — Source 3 (`p_obj_fn`), Source 4 (§19).
2. Default weights: `WT_TM_GT = WT_TM_LT = WT_SIZE_GT = WT_SIZE_LT = 1`; `WT_GC_GT = WT_GC_LT = WT_SELF_ANY = WT_SELF_END = WT_NUM_NS = 0` — Source 3.
3. Default optima: `OPT_TM = 60.0`°C, `OPT_SIZE = 20` bases, `OPT_GC_PERCENT = 50.0`% — Source 3 (60/20), Source 4 (GC 50.0).
4. `gc_content` is a percentage 0–100 (`100.0 * num_gc/num_gcat`) — Source 3 (line 3856).
5. Each term is sign- and weight-gated; the total is always ≥ 0 — Source 3.

### 1.3 Documented Corner Cases

- Parameter exactly at optimum → that parameter contributes 0 (strict `>`/`<` gates). Source 3.
- 0-weight term never contributes. Source 3.
- Penalty ≥ 0; 0 ⇔ every term at optimum. Source 3.

### 1.4 Known Failure Modes / Pitfalls

1. Treating GC as a fraction [0,1] instead of percent [0,100] scales the GC term 100× wrong — Source 3 (line 3856).
2. Using `|deviation|·singleWeight` is wrong when `WT_GT ≠ WT_LT` (one-sided weights) — Source 3/4.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `PrimerDesigner.CalculatePrimer3Penalty(Primer3PenaltyInputs, Primer3PenaltyWeights?, Primer3Optima?)` | PrimerDesigner | Canonical | Reproduces `p_obj_fn` (left/right primer) |
| `PrimerDesigner.DefaultPrimer3Weights` | PrimerDesigner | Internal | Sourced default-weight constant struct |
| `PrimerDesigner.DefaultPrimer3Optima` | PrimerDesigner | Internal | Sourced default-optima constant struct |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Penalty ≥ 0 for all inputs | Yes | Source 3 (weight·non-negative deviation) |
| INV-2 | Penalty = 0 ⇔ Tm=opt, len=opt, GC=opt, self/numNs=0 | Yes | Source 3 (sign-gated terms) |
| INV-3 | Parameter at optimum contributes 0 to its term | Yes | Source 3 (strict gates) |
| INV-4 | Term scales linearly with its weight | Yes | Source 3/4 |
| INV-5 | Default weights TM/SIZE=1, GC/SELF/NUM_NS=0; optima 60/20/50 | Yes | Source 3, Source 4 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Default penalty at optimum | Tm=60, len=20, GC=50, N=0 | 0.0 | Evidence dataset A |
| M2 | Default Tm-only deviation | Tm=63, len=20 | 3.0 | Evidence dataset B |
| M3 | Default Tm+size deviation | Tm=57, len=18 | 5.0 | Evidence dataset C |
| M4 | Default fractional deviation | Tm=62.5, len=22 | 4.5 | Evidence dataset D |
| M5 | GC `_gt` term | GC=60, WT_GC_GT=0.5 | 5.0 | Evidence dataset E |
| M6 | GC `_lt` term | GC=40, WT_GC_LT=0.5 | 5.0 | Evidence dataset F |
| M7 | self_any term | selfAny=4, WT_SELF_ANY=0.1 | 0.4 | Evidence dataset G |
| M8 | self_end term | selfEnd=3, WT_SELF_END=0.2 | 0.6 | Evidence dataset H |
| M9 | num_ns term | N=2, WT_NUM_NS=1 | 2.0 | Evidence dataset I |
| M10 | Combined multi-term | dataset J | 8.0 | Evidence dataset J |
| M11 | Default weights/optima constants | inspect default structs | TM/SIZE=1, GC/SELF/NUM_NS=0; 60/20/50 | INV-5 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Sign gate at optimum | Tm=60, WT_TM_GT=WT_TM_LT=5 | 0.0 | INV-3 |
| S2 | Asymmetric Tm weights | Tm=62, WT_TM_GT=2, WT_TM_LT=1 | 4.0 | one-sided weights |
| S3 | Linearity in weight | selfAny=4 at WT 0.1 vs 0.2 | 0.4 vs 0.8 | INV-4 |
| S4 | Non-negativity battery | several deviating inputs | all ≥ 0 | INV-1 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Lower penalty = better | two candidates | near-optimal one is smaller | "lower is better" |
| C2 | Percent-units guard | GC=50 → 0 with non-zero GC weight | 0.0 | failure mode 1 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner*.cs` cover Tm, structure, and the legacy heuristic `Score`. No test exercises a Primer3 penalty objective — the method is new for this unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M11 | ❌ Missing | new method, no prior tests |
| S1–S4 | ❌ Missing | new method |
| C1–C2 | ❌ Missing | new method |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_Primer3Penalty_Tests.cs` — all PRIMER-TM-001 penalty tests.
- **Remove:** nothing.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| PrimerDesigner_Primer3Penalty_Tests.cs | Canonical (penalty) | 17 |

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
| 11 | M11 | ❌ Missing | Implemented | ✅ Done |
| 12 | S1 | ❌ Missing | Implemented | ✅ Done |
| 13 | S2 | ❌ Missing | Implemented | ✅ Done |
| 14 | S3 | ❌ Missing | Implemented | ✅ Done |
| 15 | S4 | ❌ Missing | Implemented | ✅ Done |
| 16 | C1 | ❌ Missing | Implemented | ✅ Done |
| 17 | C2 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 17
**✅ Done:** 17 | **⛔ Blocked:** 0 | **Remaining:** must be 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | exact 0.0 |
| M2 | ✅ Covered | exact 3.0 |
| M3 | ✅ Covered | exact 5.0 |
| M4 | ✅ Covered | exact 4.5 |
| M5 | ✅ Covered | exact 5.0 |
| M6 | ✅ Covered | exact 5.0 |
| M7 | ✅ Covered | exact 0.4 |
| M8 | ✅ Covered | exact 0.6 |
| M9 | ✅ Covered | exact 2.0 |
| M10 | ✅ Covered | exact 8.0 |
| M11 | ✅ Covered | default constants asserted |
| S1 | ✅ Covered | sign gate 0.0 |
| S2 | ✅ Covered | asymmetric 4.0 |
| S3 | ✅ Covered | linearity 0.4/0.8 |
| S4 | ✅ Covered | non-negativity battery |
| C1 | ✅ Covered | lower-is-better ordering |
| C2 | ✅ Covered | percent-units guard |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | `self_any`/`self_end` alignment scores are caller-supplied (dpal local-alignment value not reproduced); the penalty arithmetic on them is exact, and default weights make these terms 0. | M7, M8, M10 |

---

## 7. Open Questions / Decisions

1. Pair-level objective (`PRIMER_PAIR_*`, Tm-difference, product-size) is out of scope; this unit reproduces the per-primer `p_obj_fn`.
