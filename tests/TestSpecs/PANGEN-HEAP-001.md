# Test Specification: PANGEN-HEAP-001

**Test Unit ID:** PANGEN-HEAP-001
**Area:** PanGenome
**Algorithm:** Pan-Genome Growth Model (Heaps' law fit + gene presence/absence matrix)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Tettelin et al. (2005), PNAS 102(39):13950–13955 | 1 | https://doi.org/10.1073/pnas.0506758102 (full text https://pmc.ncbi.nlm.nih.gov/articles/PMC1216834/) | 2026-06-13 |
| 2 | Tettelin et al. (2008), Curr Opin Microbiol 11(5):472–477 | 1 | https://doi.org/10.1016/j.mib.2008.09.006 (https://pubmed.ncbi.nlm.nih.gov/19086349/) | 2026-06-13 |
| 3 | micropan `heaps()` (R/powerlaw.R) | 3 | https://raw.githubusercontent.com/larssnip/micropan/master/R/powerlaw.R | 2026-06-13 |

### 1.2 Key Evidence Points

1. Heaps' law model for new gene clusters: `n(N) = K · N^(−alpha)` (`y.hat <- p[1]*x^(-p[2])`) — micropan powerlaw.R [3].
2. New cluster at genome i (i ≥ 2) ⇔ first appearance: `(cm==1)[i] & (cm==0)[i−1]` over the binary presence matrix; index starts at 2 — micropan powerlaw.R [3].
3. Objective minimized: `J = sqrt(Σ(y − K·x^(−alpha))²)/|x|` over K ∈ [0,10000], alpha ∈ [0,2]; start `p0 = (mean y at x=2, 1)` — micropan powerlaw.R [3].
4. Open/closed rule (verbatim): "If alpha>1.0 the pan-genome is closed, if alpha<1.0 it is open." — micropan docs [3], Tettelin 2008 [2].
5. Presence/absence is binary: counts > 0 collapse to 1 — micropan binarization [3].
6. Real anchor: S. agalactiae added 161 new genes at genome 2, 54 at genome 5, asymptote 33 (open) — Tettelin 2005 [1].

### 1.3 Documented Corner Cases

- Fewer than 2 genomes → empty new-gene curve, no fit defined [3].
- First genome contributes no "new" count (fit uses N = 2..G) [3].
- alpha boundary at 1.0 (open/closed); alpha constrained to [0,2] [3].
- Copy-number/duplicate presence collapses to 1 (binary) [3].

### 1.4 Known Failure Modes / Pitfalls

1. Using cumulative pan-genome SIZE with a positive growth exponent (P = K·N^+gamma) instead of the new-gene DECAY curve (n = K·N^−alpha) — would invert the open/closed semantics — micropan model [3].
2. Counting genes by raw id without binarizing presence — double counts duplicates — micropan binarization [3].

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FitHeapsLaw(IEnumerable<GenePresenceRow>, int)` | PanGenomeAnalyzer | Canonical | Core micropan heaps() port; deep evidence-based tests |
| `FitHeapsLaw(IReadOnlyDictionary<...>, double, int)` | PanGenomeAnalyzer | Delegate | Clusters → matrix → canonical overload |
| `CreatePresenceAbsenceMatrix(genomes, clusters)` | PanGenomeAnalyzer | Canonical | Matrix that feeds the fit; presence/absence semantics |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Heaps' fit minimizes `J=sqrt(Σ(y−K·x^(−α))²)/|x|`; for points exactly on a power curve in-bounds, recovered (K,α) match the analytic solution | Yes | micropan powerlaw.R [3] |
| INV-2 | IsOpen ⇔ alpha < 1.0 (closed when alpha > 1.0) | Yes | micropan docs [3]; Tettelin 2008 [2] |
| INV-3 | New-gene count at genome index i (i≥2) = clusters present at i but absent in all earlier genomes (first appearance) | Yes | micropan powerlaw.R [3] |
| INV-4 | Presence is binary: a cluster present any number of times in a genome counts once | Yes | micropan binarization [3] |
| INV-5 | Fitted alpha ∈ [0,2] and Intercept ∈ [0,10000] (box constraints) | Yes | micropan optim bounds [3] |
| INV-6 | predictor(N) = Intercept · N^(−alpha) is non-increasing in N for alpha ≥ 0 | Yes | model form [3] |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Closed power curve exact fit | Fixed-order matrix giving new-gene curve x=[2,3], y=[8,4] (1 perm) | Intercept ≈ 26.1640013949735, Alpha ≈ 1.70951129135145, IsOpen=false | [3] + derived (alpha=ln2/ln1.5, K=8·2^alpha) |
| M2 | Constant curve → open | Matrix giving constant new-gene count (y=[1,1]) | Alpha = 0, Intercept = 1, IsOpen = true | [3] open rule + derived |
| M3 | Open/closed boundary semantics | alpha<1 open, alpha>1 closed (M1 closed, M2 open) | IsOpen matches alpha<1 strictly | micropan rule [3] |
| M4 | New-gene first-appearance counting | Genome2 adds 8 novel clusters, Genome3 adds 4 (over shared core) | Curve y=[8,4] drives the M1 fit | INV-3 [3] |
| M5 | Binarization / duplicate presence | A cluster present in a genome via duplicate gene ids counts once | Present flag true once; no double count | INV-4 [3] |
| M6 | Fewer than 2 genomes | Single genome (and zero genomes) | Degenerate fit Intercept=0, predictor(N)=0, IsOpen=false; no exception | corner case [3] |
| M7 | Null / empty matrix | null and empty `IEnumerable<GenePresenceRow>` | Degenerate fit, no exception | corner case [3] |
| M8 | CreatePresenceAbsenceMatrix presence/absence | 2 genomes, cluster present in one only | Row.PresentGenes and GenePresence flags exact per genome | INV-3/INV-4 [3] |
| M9 | Dictionary overload delegates | Dictionary overload fits via matrix path | Same (Intercept,Alpha,IsOpen) as matrix overload on equivalent data | wrapper delegation |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Predictor monotonicity | predictor(5) ≥ predictor(10) for alpha>0 (M1 fit) | non-increasing | INV-6 |
| S2 | Bounds respected | Fitted alpha ∈ [0,2], Intercept ∈ [0,10000] | within box | INV-5 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Determinism | Same input → identical fit across calls | equal (Intercept,Alpha) | fixed-seed permutations |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/PanGenomeAnalyzerTests.cs` — contained a "Heaps Law Tests" region (`FitHeapsLaw_WithMultipleGenomes_ReturnsValidFit`, `FitHeapsLaw_PredictorWorks`, `FitHeapsLaw_TooFewGenomes_ReturnsEmpty`) and "Presence/Absence Matrix Tests" region, all written against the old non-conforming `FitHeapsLaw` (cumulative-size, gene-id matching) with permissive assertions (`GreaterThan(0)`, `GreaterThanOrEqualTo`, `.K`, `.Gamma`, `.PredictPanGenomeSize`).
- No prior canonical `PanGenomeAnalyzer_FitHeapsLaw_Tests.cs` file existed.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 closed exact fit | ❌ Missing | no exact-value Heaps fit existed |
| M2 constant → open | ❌ Missing | open/closed never asserted |
| M3 open/closed boundary | ❌ Missing | not covered |
| M4 new-gene counting | ❌ Missing | not covered |
| M5 binarization | ❌ Missing | not covered |
| M6 < 2 genomes | ⚠ Weak | old `FitHeapsLaw_TooFewGenomes_ReturnsEmpty` checked `.K==0` against removed API |
| M7 null/empty matrix | ❌ Missing | not covered |
| M8 matrix presence/absence | ⚠ Weak | old matrix tests check counts but not against the canonical-file conventions; tied to old API region |
| M9 dictionary delegation | ❌ Missing | not covered |
| S1 predictor monotonicity | ⚠ Weak | old `FitHeapsLaw_PredictorWorks` used `GreaterThanOrEqualTo` on cumulative-size model (wrong semantics) |
| S2 bounds | ❌ Missing | not covered |
| C1 determinism | ❌ Missing | not covered |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/PanGenomeAnalyzer_FitHeapsLaw_Tests.cs` — all PANGEN-HEAP-001 cases (FitHeapsLaw + CreatePresenceAbsenceMatrix).
- **Remove:** the entire "Heaps Law Tests" region and the "Presence/Absence Matrix Tests" region from `PanGenomeAnalyzerTests.cs` (they targeted the removed/old API and used permissive assertions). Other regions of that legacy file are out of scope and untouched.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `PanGenomeAnalyzer_FitHeapsLaw_Tests.cs` | Canonical PANGEN-HEAP-001 | 13 |
| `PanGenomeAnalyzerTests.cs` | Legacy (Heaps + matrix regions removed) | n/a (other regions only) |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented exact closed-curve fit | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented constant→open fit | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented open/closed assertions in M1/M2 | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented new-gene-curve matrix builder + fit | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented duplicate-presence binarization test | ✅ Done |
| 6 | M6 | ⚠ Weak | Rewrote: single + zero genomes degenerate fit | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented null + empty matrix | ✅ Done |
| 8 | M8 | ⚠ Weak | Rewrote matrix presence/absence with exact flags | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented dictionary-overload delegation | ✅ Done |
| 10 | S1 | ⚠ Weak | Rewrote predictor monotonicity on decay model | ✅ Done |
| 11 | S2 | ❌ Missing | Implemented bounds assertions | ✅ Done |
| 12 | C1 | ❌ Missing | Implemented determinism test | ✅ Done |
| 13 | legacy removal | 🔁/⚠ | Deleted old Heaps + matrix regions from legacy file | ✅ Done |

**Total items:** 13
**✅ Done:** 13 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `FitHeapsLaw_ClosedPowerCurve_RecoversExactParameters` |
| M2 | ✅ Covered | `FitHeapsLaw_ConstantNewGeneCurve_ReturnsAlphaZeroOpen` |
| M3 | ✅ Covered | open/closed asserted in M1 (closed) and M2 (open) |
| M4 | ✅ Covered | `FitHeapsLaw_CountsNewGenesByFirstAppearance` |
| M5 | ✅ Covered | `CreatePresenceAbsenceMatrix_DuplicatePresence_CountsOnce` |
| M6 | ✅ Covered | `FitHeapsLaw_FewerThanTwoGenomes_ReturnsDegenerateFit` |
| M7 | ✅ Covered | `FitHeapsLaw_NullOrEmptyMatrix_ReturnsDegenerateFit` |
| M8 | ✅ Covered | `CreatePresenceAbsenceMatrix_PresenceAbsence_ExactFlags` |
| M9 | ✅ Covered | `FitHeapsLaw_DictionaryOverload_DelegatesToMatrixFit` |
| S1 | ✅ Covered | `FitHeapsLaw_PredictNewGenes_IsNonIncreasing` |
| S2 | ✅ Covered | `FitHeapsLaw_RespectsParameterBounds` |
| C1 | ✅ Covered | `FitHeapsLaw_IsDeterministic` |

**Total in-scope cases:** 12 — **✅ 12** / ❌ 0 / ⚠ 0.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Bounded coordinate-descent optimizer minimizes the identical micropan objective over identical bounds/start; verified to match the analytic optimum < 1e-9 (optimization method is non-correctness-affecting) | M1, M2, S2 |
| 2 | Fixed-seed permutations with natural input order as first permutation make fixed-order/permutation-invariant fits exactly reproducible | M1, M2, C1 |

---

## 7. Open Questions / Decisions

1. The micropan model uses random `sample()` permutations; exact expected (K, alpha) are asserted only on permutation-invariant / fixed-order matrices where the curve is order-independent, avoiding RNG dependence. (Resolved.)
2. The checklist `FitHeapsLaw` originally fit cumulative pan-genome SIZE with a positive growth exponent; this conflicts with the authoritative new-gene DECAY model. External evidence wins — `FitHeapsLaw` was rewritten to the micropan decay model and the registry/checklist note updated. (Resolved.)
