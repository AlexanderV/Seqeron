# Test Specification: DISORDER-MORF-001

**Test Unit ID:** DISORDER-MORF-001
**Area:** ProteinPred
**Algorithm:** Molecular Recognition Feature (MoRF) Prediction
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-02-12

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Mohan et al. (2006) Analysis of Molecular Recognition Features | 1 | https://doi.org/10.1016/j.jmb.2006.07.087 | 2026-02-12 |
| 2 | Kyte & Doolittle (1982) Hydropathy scale | 1 | https://doi.org/10.1016/0022-2836(82)90515-0 | 2026-02-12 |
| 3 | Campen et al. (2008) TOP-IDP scale | 1 | https://doi.org/10.2174/092986608785849164 | 2026-02-12 |

### 1.2 Key Evidence Points

1. MoRFs are short segments (10–70 AA) within intrinsically disordered regions that undergo disorder-to-order transition upon binding — Mohan et al. (2006)
2. MoRFs are characterized by relatively higher hydrophobicity compared to their surrounding IDR context — Mohan et al. (2006)
3. This implementation uses a heuristic: hydrophobic island detection within IDRs using Kyte-Doolittle scale — not a machine-learning predictor
4. IDR detection via TOP-IDP (Campen et al. 2008) provides the context regions

### 1.3 Documented Corner Cases

1. Ordered protein (no IDRs) → no MoRFs
2. Uniformly disordered (no hydropathy contrast) → no MoRFs
3. IDR shorter than minLength → no MoRFs
4. Floating-point noise in uniform sequences → epsilon threshold (0.01 KD units) filters artifacts

### 1.4 Known Failure Modes / Pitfalls

1. True MoRFs with only marginal hydropathy enrichment may be missed (heuristic limitation)
2. Cannot distinguish MoRF subtypes (α-MoRF, β-MoRF, ι-MoRF, complex) — that requires structural information
3. The "hydrophobic island" heuristic is an annotation tool, not validated against crystallographic MoRF datasets

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `PredictMoRFs(sequence, minLength, maxLength)` | DisorderPredictor | Public API | MoRF prediction entry point |
| `PredictDisorder(sequence)` | DisorderPredictor | Public dependency | Provides IDR context for MoRF scanning |
| `MeanHydropathy(sequence, start, length)` | DisorderPredictor | Private helper | KD hydropathy for windows |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | All MoRFs have Start ≥ 0, End < sequence.Length, Start ≤ End, Score ∈ [0, 1] | Yes | Algorithm definition |
| INV-2 | Every MoRF is contained within a predicted IDR | Yes | By construction (scans within IDRs only) |
| INV-3 | No two MoRFs overlap | Yes | Greedy non-overlapping selection |

---

## 4. Test Cases

### 4.1 MUST Tests (Required)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| S1 | OrderedProtein_NoMoRFs | 40×L → no IDRs → no MoRFs | 0 MoRFs | No IDRs from Campen (2008) |
| S2 | UniformIDR_NoMoRFs | 40×P → uniform KD → no hydropathy contrast | 0 MoRFs | Mohan (2006): MoRFs need enrichment |
| S3 | HydrophobicIsland_Detected | G-island (KD −0.4) in E/K context (KD −3.7) | ≥1 MoRF near island, Score > 0.5 | Mohan (2006) + KD scale |
| M1 | ScoreReflectsContrast | G-island (−0.4) vs S-island (−0.8) in same context | G scores higher than S | KD values |
| M2 | MaxLengthRespected | All MoRFs ≤ maxLength | Length ≤ maxLength | API contract |
| M3 | MinLengthRespected | All MoRFs ≥ minLength | Length ≥ minLength | API contract |

### 4.2 SHOULD Tests (Edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | ShortSequence | 5 AA → too short for IDR + MoRF | 0 MoRFs | Minimum viability |
| C2 | IDRShorterThanMinLen | 8-residue IDR with minLength=10 | 0 MoRFs | Window doesn't fit |
| M4 | CaseInsensitive | Lower vs upper input → same results | Same count | API contract |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_MoRF_Tests.cs` — 14 tests for DISORDER-MORF-001.
- Previous tests in `DisorderPredictorTests.cs` removed (used old fabricated scoring API).

### 5.2 Coverage Classification

| Test Case ID | Status | Notes |
|--------------|--------|-------|
| S1: Ordered protein | ✅ Covered | Is.Empty |
| S2: Uniform IDR | ✅ Covered | Is.Empty (epsilon filter) |
| S3: Hydrophobic island | ✅ Covered | Count ≥ 1, position range, Score > 0.5 |
| C1: Short sequence | ✅ Covered | Is.Empty |
| C2: IDR < minLength | ✅ Covered | Is.Empty |
| M1: Score contrast | ✅ Covered | G-score > S-score |
| M2: Max length | ✅ Covered | All lengths ≤ maxLen |
| M3: Min length | ✅ Covered | All lengths ≥ minLen |
| M4: Case insensitive | ✅ Covered | Same count |
| INV-1: Valid boundaries | ✅ Covered | 3 diverse inputs |
| INV-2: Within IDR | ✅ Covered | SeqWithIsland |
| INV-3: Non-overlapping | ✅ Covered | SeqWithIsland |

### 5.3 Final State

| File | Role | Test Count |
|------|------|------------|
| `DisorderPredictor_MoRF_Tests.cs` | DISORDER-MORF-001 canonical | 14 |

---

## 6. Evidence Traceability

| Parameter | Value | Source | Status |
|-----------|-------|--------|--------|
| MoRF definition | 10–70 AA disorder-to-order transition | Mohan et al. (2006) | ✅ Published |
| Hydropathy scale | Kyte-Doolittle 20 AA values | Kyte & Doolittle (1982) | ✅ Published |
| IDR detection | TOP-IDP, cutoff 0.542, window 21 | Campen et al. (2008) | ✅ Published |
| Detection criterion | window KD > IDR KD + 0.01 | **Internal heuristic** | ⚠️ Design decision |
| Score normalization | min(1, diff / 3.0) | **Internal heuristic** | ⚠️ Design decision |
| Default minLength | 10 | **Internal choice** — Mohan et al. report 5–25 AA range | ⚠️ Design decision |
| Default maxLength | 25 | Mohan et al. (2006) — most MoRFs ≤ 25 AA | ✅ Approximate |
| Greedy selection | Non-overlapping, best-score-first | **Internal heuristic** | ⚠️ Design decision |

---

## 7. Deviations and Assumptions

| ID | Item | Status | Detail |
|----|------|--------|--------|
| D1 | Heuristic not ML | ⚠️ Limitation | True MoRF predictors (MoRFpred, MoRFchibi, ANCHOR) use ML/energy models. This heuristic uses only KD hydropathy enrichment as an annotation tool. |
| D2 | Score normalization | ⚠️ Internal | Score = min(1, diff / 3.0) where diff = windowKD − idrKD. The divisor 3.0 is an internal design choice (roughly the maximum plausible KD difference for a hydrophobic island in an IDR context). |
| D3 | Epsilon threshold | ⚠️ Internal | windowKD > idrKD + 0.01 KD units. The 0.01 threshold prevents floating-point noise from producing false MoRFs in uniform IDRs. |
| D4 | Default maxLength = 25 | ⚠️ Approximate | Mohan et al. (2006) report MoRFs of 5–25 AA in their crystal structure survey. The 25 upper bound is approximate — some MoRFs may be longer. |
