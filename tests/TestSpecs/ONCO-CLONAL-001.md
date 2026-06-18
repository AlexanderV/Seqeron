# Test Specification: ONCO-CLONAL-001

**Test Unit ID:** ONCO-CLONAL-001
**Area:** Oncology
**Algorithm:** Clonal vs Subclonal Mutation Classification (cancer cell fraction posterior)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Landau DA et al. (2013). Evolution and Impact of Subclonal Mutations in CLL. *Cell* 152(4):714–726. | 1 | https://doi.org/10.1016/j.cell.2013.01.019 | 2026-06-14 |
| 2 | Satas G et al. (2021). DeCiFering the Elusive Cancer Cell Fraction. *Cell Systems* 12(10):1004–1018. | 1 | https://doi.org/10.1016/j.cels.2021.07.006 | 2026-06-14 |

### 1.2 Key Evidence Points

1. Clonal iff "the CCF harboring it was >0.95 with probability > 0.5, and subclonal otherwise" — Landau (2013).
2. Expected allele fraction `f(c) = αc/(2(1−α)+αq)`; posterior `P(c) ∝ Binom(a|N,f(c))`, uniform prior, grid of 100 c values, c∈[0.01,1] — Landau (2013).
3. Multiplicity-general CCF `c ≈ (1/ρ)·(ρ·N_tot + 2(1−ρ))/M · v̂` (Eq. 1) ⇒ `f(c) = ρMc/(2(1−ρ)+ρq)` — Satas (2021).
4. Clonal ≈ present in all cells (CCF ≈ 1); subclonal ≈ present in a subpopulation (CCF ≪ 1) — Satas (2021).

### 1.3 Documented Corner Cases

- CCF grid lower bound is 0.01 (a detected mutation is in ≥ 1 cancer cell) — Landau (2013).
- Classification is by posterior mass above 0.95, not the point estimate — shallow coverage can leave a near-1 point estimate subclonal — Landau (2013).
- Multiplicity M > 1 lowers CCF for the same VAF — Satas (2021).

### 1.4 Known Failure Modes / Pitfalls

1. Using the CCF point estimate instead of the posterior probability would misclassify uncertain variants — Landau (2013).
2. Ignoring multiplicity overestimates CCF for multi-copy mutant loci — Satas (2021).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `ClassifyClonality(variants, purity)` | OncologyAnalyzer | Canonical | Posterior-grid clonal/subclonal classification per Landau (2013). |
| `IdentifyClonalMutations(ccfValues)` | OncologyAnalyzer | Canonical | Point-estimate clonal selection: CCF > 0.95. |
| `ClonalityVariant`, `ClonalityCall`, `ClonalityResult`, `ClonalityStatus` | OncologyAnalyzer | Internal | Result/record types tested via the two methods above. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | ClonalCount + SubclonalCount = total variants | Yes | Registry invariant; partition of a Boolean classification |
| INV-2 | ClonalFraction = ClonalCount / total (0 for empty set) | Yes | Definition of clonal fraction summary |
| INV-3 | Each call's CCF mean ∈ [0.01, 1] | Yes | Landau (2013) grid c ∈ [0.01, 1] |
| INV-4 | Each call's ProbabilityClonal ∈ [0, 1] | Yes | Normalised posterior probability |
| INV-5 | A variant is Clonal iff ProbabilityClonal > 0.5 (P(CCF>0.95)) | Yes | Landau (2013) rule |
| INV-6 | `IdentifyClonalMutations` returns index i iff ccf[i] > 0.95 (strict) | Yes | Landau (2013), CCF > 0.95 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Clonal pure het deep | a=300, N=300, q=2, M=1, ρ=1.0 | Clonal; CCF≈0.999486; P=1.0 | Landau (2013) |
| M2 | Clonal het impure deep | a=400, N=1000, q=2, M=1, ρ=0.8 | Clonal; CCF≈0.972455; P≈0.864167 | Landau (2013) |
| M3 | Subclonal CCF~0.6 | a=240, N=1000, q=2, M=1, ρ=0.8 | Subclonal; CCF≈0.601297; P≈0.0 | Landau (2013) |
| M4 | Subclonal CCF~0.4 | a=200, N=1000, q=2, M=1, ρ=1.0 | Subclonal; CCF≈0.401198; P≈0 | Landau (2013) |
| M5 | Multiplicity M=2 raises CCF | a=100, N=100, q=2, M=2, ρ=1.0 | Clonal; CCF≈0.994330; P≈0.998016 | Satas (2021) Eq. 1 |
| M6 | Counts partition (INV-1) | mix of M1+M3 variants | ClonalCount=1, SubclonalCount=1, total=2 | Registry invariant |
| M7 | ClonalFraction (INV-2) | 3 clonal + 1 subclonal | ClonalFraction = 0.75 | Definition |
| M8 | IdentifyClonalMutations strict >0.95 | {0.96,0.95,1.0,0.5,0.951} | indices {0,2,4} | Landau (2013) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Null variants | ClassifyClonality(null, 0.8) | ArgumentNullException | Standard validation |
| S2 | Purity out of range | purity 0 / 1.5 / NaN | ArgumentOutOfRangeException | Domain ρ ∈ (0,1] |
| S3 | Invalid read counts | TotalReads 0; AltReads > TotalReads | ArgumentException | Domain |
| S4 | Invalid copy / multiplicity | LocalCopyNumber 0; Multiplicity 3 at q=2 | ArgumentException | Domain |
| S5 | Null CCF list | IdentifyClonalMutations(null) | ArgumentNullException | Standard validation |
| S6 | CCF out of range | {1.2} or {NaN} | ArgumentException | Domain [0,1] |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Empty variant set | ClassifyClonality([], 0.8) | Empty calls, counts 0, ClonalFraction 0 | Degenerate input |
| C2 | Empty CCF list | IdentifyClonalMutations([]) | empty index list | Degenerate input |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No prior tests for `ClassifyClonality` / `IdentifyClonalMutations`; both methods are new in this unit. Canonical file created: `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ClassifyClonality_Tests.cs`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M8 | ❌ Missing | New unit; no prior tests. |
| S1–S6 | ❌ Missing | New unit. |
| C1–C2 | ❌ Missing | New unit. |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ClassifyClonality_Tests.cs` — all cases for both methods.
- **Remove:** none (no pre-existing tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| OncologyAnalyzer_ClassifyClonality_Tests.cs | Canonical | 16 |

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
| 9 | S1 | ❌ Missing | Implemented | ✅ Done |
| 10 | S2 | ❌ Missing | Implemented | ✅ Done |
| 11 | S3 | ❌ Missing | Implemented | ✅ Done |
| 12 | S4 | ❌ Missing | Implemented | ✅ Done |
| 13 | S5 | ❌ Missing | Implemented | ✅ Done |
| 14 | S6 | ❌ Missing | Implemented | ✅ Done |
| 15 | C1 | ❌ Missing | Implemented | ✅ Done |
| 16 | C2 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 16
**✅ Done:** 16 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | Clonal pure-het deep coverage |
| M2 | ✅ Covered | Clonal impure-het deep coverage |
| M3 | ✅ Covered | Subclonal CCF~0.6 |
| M4 | ✅ Covered | Subclonal CCF~0.4 |
| M5 | ✅ Covered | Multiplicity M=2 |
| M6 | ✅ Covered | Counts partition INV-1 |
| M7 | ✅ Covered | ClonalFraction INV-2 |
| M8 | ✅ Covered | IdentifyClonalMutations strict >0.95 |
| S1 | ✅ Covered | Null variants |
| S2 | ✅ Covered | Purity out of range (3 cases) |
| S3 | ✅ Covered | Invalid read counts |
| S4 | ✅ Covered | Invalid copy/multiplicity |
| S5 | ✅ Covered | Null CCF list |
| S6 | ✅ Covered | CCF out of range |
| C1 | ✅ Covered | Empty variant set |
| C2 | ✅ Covered | Empty CCF list |

Total in-scope cases: 16; ✅ = 16.

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Registry `ploidy` param = per-variant local copy number q (API shape only; not correctness-affecting). | Method signature `ClassifyClonality(variants, purity)` |

---

## 7. Open Questions / Decisions

1. ONCO-CCF-001 (`EstimateCCF`, `ClusterCCFValues`) is a separate later unit; this unit does the clonal/subclonal classification only. `IdentifyClonalMutations(ccfValues)` consumes pre-computed CCFs (the bridge to CCF-001) but does not itself estimate CCF beyond the internal posterior used by `ClassifyClonality`.
