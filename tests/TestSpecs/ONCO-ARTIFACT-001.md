# Test Specification: ONCO-ARTIFACT-001

**Test Unit ID:** ONCO-ARTIFACT-001
**Area:** Oncology
**Algorithm:** Sequencing Artifact Detection (OxoG / FFPE deamination substitution classification + Fisher strand-bias)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Chen L. et al. (2017). DNA damage is a pervasive cause of sequencing errors. Science 355:752–756 | 1 | https://www.science.org/doi/10.1126/science.aai8690 | 2026-06-14 |
| 2 | Ettwiller Damage-estimator (reference impl; GIV score) | 3 | https://github.com/Ettwiller/Damage-estimator | 2026-06-14 |
| 3 | Nature Methods (2017). DNA variants or DNA damage? (GIV thresholds) | 1 | https://www.nature.com/articles/nmeth.4254 | 2026-06-14 |
| 4 | Do & Dobrovic (2015). FFPE deamination artifacts | 1 | https://www.sciencedirect.com/science/article/pii/S152515781630188X | 2026-06-14 |
| 5 | GATK FisherStrand / StrandBiasTest (Broad) | 3 | https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/annotator/StrandBiasTest.java | 2026-06-14 |

### 1.2 Key Evidence Points

1. OxoG artifact = G>T excess in read 1 / C>A excess in read 2 (reverse complement) — Chen 2017.
2. GIV_G_T = (count G>T in R1) / (count G>T in R2); GIV = 1 undamaged, GIV > 1.5 damaged — Chen 2017 / Nature Methods.
3. FFPE deamination artifact = C>T / G>A (collectively C:G>T:A), uracil pairs with adenine — Do & Dobrovic 2015.
4. FFPE (C>T/G>A) and oxidation (G>T/C>A) are disjoint substitution classes — Do & Dobrovic 2015.
5. FisherStrand FS = -10·log10(two-sided Fisher exact p) on the 2×2 table [ref_fwd, ref_rev, alt_fwd, alt_rev]; floor MIN_PVALUE = 1E-320 — GATK.

### 1.3 Documented Corner Cases

- No strand bias (balanced table) ⇒ Fisher p ≈ 1 ⇒ FS ≈ 0 (GATK).
- Zero G>T in R2 ⇒ GIV ratio undefined; both zero ⇒ no imbalance (treated as GIV = 1); R2 zero with R1 > 0 ⇒ maximal imbalance (damaged).
- Substitution outside {C>T,G>A,G>T,C>A} ⇒ not an artifact class.

### 1.4 Known Failure Modes / Pitfalls

1. Treating C>T/G>A (FFPE) as oxidation, or G>T/C>A (OxoG) as deamination — they are disjoint classes (Do & Dobrovic 2015).
2. Division by zero in GIV when R2 has no G>T reads (Chen 2017 / Damage-estimator).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FilterArtifacts(variants, ...)` | OncologyAnalyzer | Canonical | Removes variants flagged as sequencing artifacts |
| `DetectOxoGArtifacts(variants, ...)` | OncologyAnalyzer | Canonical | OxoG G>T/C>A detection via GIV |
| `ClassifyArtifact(observation)` | OncologyAnalyzer | Canonical | Substitution-class + strand classification of one variant |
| `CalculateGivScore(r1Count, r2Count)` | OncologyAnalyzer | Canonical | GIV ratio R1/R2 |
| `CalculateStrandBias(refFwd, refRev, altFwd, altRev)` | OncologyAnalyzer | Canonical | Phred-scaled Fisher strand FS |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | `FilterArtifacts` result ⊆ input (subset, input order) | Yes | composition (Chen 2017; GATK) |
| INV-2 | GIV ≥ 0; GIV = 1 for balanced R1=R2 | Yes | Chen 2017 / Nature Methods |
| INV-3 | FS ≥ 0; FS = 0 for a perfectly balanced strand table (p = 1) | Yes | GATK FisherStrand |
| INV-4 | FFPE class {C>T, G>A} and OxoG class {G>T, C>A} are disjoint | Yes | Do & Dobrovic 2015 |
| INV-5 | Higher strand segregation ⇒ FS non-decreasing (monotone) | Yes | GATK (Fisher p decreases) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | FFPE C>T | classify C→T substitution | FfpeDeamination | Do & Dobrovic 2015 |
| M2 | FFPE G>A | classify G→A substitution | FfpeDeamination | Do & Dobrovic 2015 |
| M3 | OxoG G>T | classify G→T substitution | OxoG | Chen 2017 |
| M4 | OxoG C>A | classify C→A substitution | OxoG | Chen 2017 |
| M5 | Non-artifact A>G | classify A→G | None | classes are specific (Do & Dobrovic; Chen) |
| M6 | GIV damaged | R1=200, R2=100 | GIV = 2.0 (> 1.5) | Chen 2017 / Damage-estimator |
| M7 | GIV undamaged | R1=R2=100 | GIV = 1.0 | Chen 2017 / Nature Methods |
| M8 | FS balanced | table [10,10,10,10] | FS = 0.000 (p = 1) | GATK FisherStrand |
| M9 | FS segregated | table [20,0,0,20] | FS = -10·log10(p), p computed exactly ⇒ FS > 0 | GATK FisherStrand |
| M10 | DetectOxoGArtifacts | variant set with one G>T at GIV 2.0 | that variant flagged OxoG | Chen 2017 |
| M11 | FilterArtifacts removes | FFPE C>T artifact + real A>G variant | only A>G kept; result ⊆ input | composition; INV-1 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | GIV zero R2 | R1=50, R2=0 | GIV = +∞-sentinel (damaged); no exception | corner case |
| S2 | GIV both zero | R1=0, R2=0 | GIV = 1.0 (no imbalance) | corner case |
| S3 | null input | FilterArtifacts(null) | ArgumentNullException | API contract |
| S4 | empty input | FilterArtifacts(empty) | empty result | API contract |
| S5 | FS negative count | CalculateStrandBias(-1,...) | ArgumentOutOfRangeException | validation |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | FS monotonicity (property) | increasing segregation | FS non-decreasing | INV-5 |
| C2 | GIV symmetry | R1=100,R2=200 | GIV = 0.5 (< 1, no OxoG on this strand) | ratio is directional |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing tests for `FilterArtifacts` / `DetectOxoGArtifacts` (methods did not exist). Sibling Oncology fixtures: `OncologyAnalyzer_CallSomaticMutations_Tests.cs`, `OncologyAnalyzer_CalculateVAF_Tests.cs`, `OncologyAnalyzer_IdentifyDriverMutations_Tests.cs`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M11 | ❌ Missing | new unit, no prior tests |
| S1–S5 | ❌ Missing | new unit |
| C1–C2 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_FilterArtifacts_Tests.cs` — all cases for this unit.
- **Remove:** nothing (no prior tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `OncologyAnalyzer_FilterArtifacts_Tests.cs` | canonical | 18 |

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
| 12 | S1 | ❌ Missing | implemented | ✅ Done |
| 13 | S2 | ❌ Missing | implemented | ✅ Done |
| 14 | S3 | ❌ Missing | implemented | ✅ Done |
| 15 | S4 | ❌ Missing | implemented | ✅ Done |
| 16 | S5 | ❌ Missing | implemented | ✅ Done |
| 17 | C1 | ❌ Missing | implemented | ✅ Done |
| 18 | C2 | ❌ Missing | implemented | ✅ Done |

**Total items:** 18
**✅ Done:** 18 | **⛔ Blocked:** 0 | **Remaining:** must be 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1–M11 | ✅ Covered | exact evidence-based assertions |
| S1–S5 | ✅ Covered | corner cases / validation |
| C1–C2 | ✅ Covered | property + directional ratio |

In-scope cases: 18. ✅ = 18.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | No BAM parser — strand/read-mate counts passed on the observation record (API-shape only) | method signatures |
| 2 | GIV neutral = 1, damaged > 1.5 (verbatim from Chen 2017 / Nature Methods) | DetectOxoGArtifacts, M6/M7 |

---

## 7. Open Questions / Decisions

1. The exact two-sided Fisher p for the segregated table M9 is computed analytically in the test from the hypergeometric distribution (so the expected FS is independent of the implementation). FS = 0 for the balanced table M8 is exact (p = 1).
