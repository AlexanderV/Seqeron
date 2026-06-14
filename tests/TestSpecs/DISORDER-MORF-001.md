# Test Specification: DISORDER-MORF-001

**Test Unit ID:** DISORDER-MORF-001
**Area:** ProteinPred
**Algorithm:** MoRF (Molecular Recognition Feature) Prediction
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Mohan et al. (2006) Analysis of MoRFs, J Mol Biol | 1 | https://doi.org/10.1016/j.jmb.2006.07.087 (PMID 16935303) | 2026-06-14 |
| 2 | Cheng/Oldfield, Mining α-MoRFs, Biochemistry | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC2570644/ | 2026-06-14 |
| 3 | Oldfield et al. (2005) Coupled folding/binding, Biochemistry | 1 | https://pubmed.ncbi.nlm.nih.gov/16156658/ | 2026-06-14 |
| 4 | Wikipedia, Molecular recognition feature (cites Mohan 2006) | 4 | https://en.wikipedia.org/wiki/Molecular_recognition_feature | 2026-06-14 |
| 5 | Campen et al. (2008) TOP-IDP scale (per-residue scores) | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC2676888/ | 2026-06-14 |

### 1.2 Key Evidence Points

1. A MoRF is a short region of relative **order** ("dip") embedded **within a longer region of disorder** — Cheng/Oldfield (PMC2570644): "short regions of order within longer regions of disorder – or 'dips'".
2. MoRF length range is **10–70 residues** — Mohan (2006): "relatively short (10-70 residues)".
3. The order/disorder boundary in the prediction profile is the **0.5 threshold** — Cheng/Oldfield (PMC2570644): "the threshold of 0.5"; below it = predicted order (the dip).
4. A MoRF is **disordered prior to binding** and undergoes a disorder-to-order transition on binding — Wikipedia (citing Mohan 2006); Oldfield (2005).
5. Per-residue disorder scores come from the normalized TOP-IDP scale already used by `PredictDisorder` (range [0,1], higher = more disordered) — Campen (2008) Table 2.

### 1.3 Documented Corner Cases

- Fully ordered protein → no surrounding disorder → no MoRF (PMC2570644).
- Fully disordered protein → no ordered dip → no MoRF (PMC2570644).
- Region outside 10–70 residues → not a MoRF (Mohan 2006).
- Dip at sequence terminus, not flanked by disorder on both sides → not "within a region of disorder" → not a MoRF (Oldfield 2005 / Mohan 2006).

### 1.4 Known Failure Modes / Pitfalls

1. Exact dip flank/run-length numeric parameters are paywalled (Oldfield 2005 Methods) and not retrievable — the implemented criterion uses the retrievable qualitative shape (ordered run flanked by disorder, 10–70 length, 0.5 threshold). See §6.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `PredictMoRFs(sequence, minLength, maxLength)` | DisorderPredictor | **Canonical** | Dip-in-disorder detector |
| `PredictDisorder(sequence)` | DisorderPredictor | **Internal** | Provides per-residue disorder scores; tested indirectly |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every MoRF has `0 ≤ Start ≤ End < sequence.Length` | Yes | Coordinate contract |
| INV-2 | Every MoRF length is within `[minLength, maxLength]` (default 10–70) | Yes | Mohan 2006 length range |
| INV-3 | Every MoRF lies inside a predicted disordered region (flanked by disorder) | Yes | PMC2570644 "within disorder" |
| INV-4 | Reported MoRFs are non-overlapping and ordered by `Start` | Yes | Distinct regions |
| INV-5 | `0 ≤ Score ≤ 1` and Score increases with dip depth below 0.5 | Yes | 0.5 threshold (PMC2570644); bounded derivation |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Ordered dip in disorder | 15 ordered L (score 0.298) flanked by 20 disordered P each side | Exactly one MoRF covering the L run (positions 20–34) | PMC2570644 "dips"; Mohan length |
| M2 | Fully ordered sequence | 40 L (all ordered) | No MoRFs (no surrounding disorder) | PMC2570644 |
| M3 | Fully disordered sequence | 40 P (all disordered) | No MoRFs (no ordered dip) | PMC2570644 |
| M4 | Dip too short | 8 ordered L flanked by 20 P each side | No MoRFs (< 10 residues) | Mohan 2006 length |
| M5 | Dip too long | 80 ordered L flanked by 20 P each side | No MoRFs (> 70 residues) | Mohan 2006 length |
| M6 | Dip at terminus | 15 ordered L at start, then 30 P | No MoRFs (not flanked by disorder both sides) | Oldfield 2005 "within disorder" |
| M7 | Score bounded & depth-monotone | Deep dip (I, score 0.213) scores higher than shallow dip (L, score 0.298); both in [0,1] | I-MoRF Score > L-MoRF Score; both within [0,1] | 0.5 threshold derivation |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Two separate dips | Two L runs separated by disorder | Two non-overlapping MoRFs | Region independence |
| S2 | Case-insensitive | Lower-case M1 sequence | Same result as upper-case | Sibling methods upper-case |
| S3 | minLength/maxLength respected | Custom bounds shrink/grow eligibility | All MoRFs within bounds | INV-2 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Null / empty | null and "" | Empty result | Standard guard |
| C2 | Shorter than minLength | 5-residue sequence | Empty result | Cannot fit a MoRF |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_MoRF_Tests.cs` existed, encoding the **previous hydropathy-enrichment heuristic** (non-canonical; no retrieved source defines hydropathy-vs-IDR-context as the MoRF criterion, and it used invented constants 3.0 and 0.01).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| Old S1/S2 (no MoRF on ordered/uniform) | ⚠ Weak | Right intent but encoded hydropathy heuristic, not the dip definition |
| Old S3 (hydrophobic island) | ⚠ Weak | Non-source criterion; permissive `GreaterThanOrEqualTo`, `InRange`, `GreaterThan` |
| Old M1 (hydropathy contrast) | ⚠ Weak | Conditional `if`, no exact value; non-source criterion |
| Old M2/M3 (min/max length) | ⚠ Weak | Permissive; reuses hydropathy construct |
| Old M4 (case) | ⚠ Weak | Count-only |
| Old INV1–3 | ⚠ Weak | Permissive; reuse hydropathy construct |
| M1–M7, S1–S3, C1–C2 (this spec) | ❌ Missing | New dip-in-disorder cases |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_MoRF_Tests.cs` — rewritten from scratch to encode the dip-in-disorder definition with exact coordinates and exact ordered/disorder scores.
- **Remove:** all previous hydropathy-enrichment tests in that file (replaced).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `DisorderPredictor_MoRF_Tests.cs` | Canonical, evidence-based | 14 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented exact-coordinate dip test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | S1 | ❌ Missing | Implemented | ✅ Done |
| 9 | S2 | ❌ Missing | Implemented | ✅ Done |
| 10 | S3 | ❌ Missing | Implemented | ✅ Done |
| 11 | C1 | ❌ Missing | Implemented | ✅ Done |
| 12 | C2 | ❌ Missing | Implemented | ✅ Done |
| 13 | INV-1..INV-5 | ❌ Missing | Implemented property tests | ✅ Done |
| 14 | Old hydropathy tests | ⚠ Weak | Removed (file rewritten) | ✅ Done |

**Total items:** 14
**✅ Done:** 14 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | Exact coordinates 20–34 |
| M2 | ✅ | Empty |
| M3 | ✅ | Empty |
| M4 | ✅ | Empty (too short) |
| M5 | ✅ | Empty (too long) |
| M6 | ✅ | Empty (terminus) |
| M7 | ✅ | Depth-monotone, bounded |
| S1 | ✅ | Two MoRFs |
| S2 | ✅ | Case-insensitive |
| S3 | ✅ | Bounds respected |
| C1 | ✅ | Null/empty empty |
| C2 | ✅ | Short empty |
| INV-1..INV-5 | ✅ | Property tests |
| Old hydropathy tests | ✅ | Removed |

In-scope cases: 12 planned + 5 invariants. All ✅.

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Exact dip flank/run-length parameters (Oldfield 2005 Methods, paywalled) approximated by "ordered run flanked by ≥1 disordered residue both sides"; threshold 0.5 and length 10–70 are source-traceable, not assumed | Implementation flank rule |

---

## 7. Open Questions / Decisions

1. Decided: reimplement `PredictMoRFs` from the previous hydropathy-enrichment heuristic to the source-traceable dip-in-disorder definition. The previous heuristic had an untraceable criterion (hydropathy vs IDR context) and invented normalization constants (3.0, 0.01) — defects under the Implementation conformance policy. The new score normalizes by the 0.5-threshold dip depth, which is source-derived.
