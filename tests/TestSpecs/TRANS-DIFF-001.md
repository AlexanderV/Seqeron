# Test Specification: TRANS-DIFF-001

**Test Unit ID:** TRANS-DIFF-001
**Area:** Transcriptome
**Algorithm:** Differential Expression (log2 fold change, Welch's t-test, Benjamini-Hochberg FDR)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Love, Huber & Anders (2014) — DESeq2, Genome Biology 15:550 | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC4302049/ | 2026-06-13 |
| 2 | Benjamini & Hochberg (1995), JRSS-B 57(1):289–300 | 1 (via 2/4) | https://doi.org/10.1111/j.2517-6161.1995.tb02031.x ; https://en.wikipedia.org/wiki/False_discovery_rate ; https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html | 2026-06-13 |
| 3 | Welch (1947), Biometrika 34 (t-test) | 4 (via Wikipedia) | https://en.wikipedia.org/wiki/Welch%27s_t-test | 2026-06-13 |
| 4 | Student's t-distribution CDF (regularized incomplete beta) | 4 | https://en.wikipedia.org/wiki/Student%27s_t-distribution | 2026-06-13 |
| 5 | Science Park RNA-seq lesson — Differential expression | 3 | https://scienceparkstudygroup.github.io/rna-seq-lesson/06-differential-analysis/index.html | 2026-06-13 |

### 1.2 Key Evidence Points

1. log2 fold change = `log₂(mean_treatment / mean_control)`; positive ⇒ up in treatment (numerator) — DESeq2 (Love 2014); Science Park lesson.
2. DESeq2 reports model coefficients on the log₂ scale and adjusts p-values "using the procedure of Benjamini and Hochberg" — Love et al. (2014).
3. A gene is differentially expressed only when BOTH |log2FC| ≥ threshold AND adjusted p-value < alpha — Science Park lesson; DESeq2.
4. Welch's t-statistic `t = (X̄₂−X̄₁)/√(s²₁/N₁+s²₂/N₂)` with unbiased (N−1) variances; Welch-Satterthwaite df `ν = (s²₁/N₁+s²₂/N₂)² / [s⁴₁/(N²₁(N₁−1))+s⁴₂/(N²₂(N₂−1))]` — Welch (1947) via Wikipedia.
5. Two-sided p-value `P(|T|≥t) = I_{ν/(ν+t²)}(ν/2, ½)` (regularized incomplete beta) — Student's t-distribution CDF.
6. BH adjusted p-value (R `p.adjust`): sort p descending, multiply by `n/rank` (rank m..1), take running minimum, clamp to 1, restore order: `pmin(1, cummin(n/i·p[o]))[ro]` — Benjamini & Hochberg (1995) / R stats.

### 1.3 Documented Corner Cases

- Two-criterion gate: failing either |log2FC| or adjusted-p criterion ⇒ not significant (Science Park / DESeq2).
- Group with <2 replicates: variance undefined ⇒ gene not testable ⇒ p = 1 (Welch precondition).
- se = 0: equal means ⇒ p = 1; unequal means ⇒ p = 0 (limit of the t-statistic).
- Zero mean expression: `log₂(mean2/mean1)` undefined ⇒ pseudocount added to both means.

### 1.4 Known Failure Modes / Pitfalls

1. Reporting raw p-values across many genes inflates false positives — BH adjustment required (Science Park lesson; DESeq2).
2. Using the normal approximation instead of the t-distribution mis-estimates p-values at small N — addressed by the exact `I_x` t-CDF (Student's t-distribution source).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindDifferentiallyExpressed(condition1, condition2, alpha, log2FoldChangeThreshold)` | TranscriptomeAnalyzer | **Canonical** | Welch t-test + BH FDR + two-criterion DE gate |
| `CalculateFoldChange(expression1, expression2)` | TranscriptomeAnalyzer | **Canonical** | log₂((mean2+c)/(mean1+c)), sign = up-in-2 |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-01 | `CalculateFoldChange(a,b) = −CalculateFoldChange(b,a)` (sign symmetry) | Yes | log-ratio definition, DESeq2 / Science Park |
| INV-02 | Equal group means ⇒ log2FC = 0 and p-value = 1 ⇒ not significant | Yes | log-ratio = 0; t = 0 ⇒ p = 1 |
| INV-03 | Adjusted p-value ≥ raw p-value and ≤ 1 for every gene | Yes | BH step-up clamps to ≤1 and inflates by n/rank ≥ 1 |
| INV-04 | Adjusted p-values are monotone non-decreasing in raw-p ascending order | Yes | BH cumulative-minimum (R p.adjust) |
| INV-05 | A gene is significant iff |log2FC| ≥ threshold AND adjusted p < alpha | Yes | two-criterion gate, Science Park / DESeq2 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | FoldChange UP | control {10,10,10}, treatment {40,40,40} | log2FC = +1.8981204 (log₂(41/11)) | DESeq2; Science Park |
| M2 | FoldChange sign symmetry | swap arguments | DOWN = −1.8981204 = −UP (INV-01) | DESeq2; Science Park |
| M3 | FoldChange flat | equal means | log2FC = 0 | log-ratio definition |
| M4 | Welch t-stat & df | control {1,2,3}, treatment {7,8,9} | t = 7.3484692, ν = 4, p = 0.0018262607 | Welch 1947; Student t CDF; SciPy cross-check |
| M5 | BH adjusted p-values | raw (0.001,0.4,0.5,0.9) | adjusted (0.004,0.66667,0.66667,0.9) | BH 1995 / R p.adjust |
| M6 | BH monotone & ≥ raw | same set | adjusted non-decreasing in p-order; each ≥ raw, ≤1 (INV-03,04) | BH 1995 |
| M7 | DE two-criterion gate | strong DE gene vs flat gene | strong gene IsSignificant=true (|log2FC|≥thr & adjP<alpha); flat gene false | Science Park; DESeq2 |
| M8 | DE regulation label | UP gene vs DOWN gene | "Upregulated" (log2FC>0) / "Downregulated" (log2FC<0) | DESeq2 sign convention |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Empty input | no genes | empty result | degenerate |
| S2 | Identical groups | control = treatment | log2FC=0, p=1, not significant (INV-02) | both criteria fail |
| S3 | <2 replicates | one replicate per group | p = 1, not significant | Welch precondition / failure mode |
| S4 | Threshold gate | strong p but |log2FC| below threshold | not significant | two-criterion gate (INV-05) |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | FoldChange zero mean | control mean 0, treatment positive | finite log2FC (pseudocount), positive sign | degenerate pseudocount convention |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No test file exists for TRANS-DIFF-001. Searched `tests/Seqeron/Seqeron.Genomics.Tests/` — sibling file `TranscriptomeAnalyzer_ExpressionQuantification_Tests.cs` (TRANS-EXPR-001) only; no `FindDifferentiallyExpressed`/`CalculateFoldChange` tests.
- Production: `TranscriptomeAnalyzer` had `AnalyzeDifferentialExpression` (overlapping concept) but the registry-canonical methods `FindDifferentiallyExpressed` and `CalculateFoldChange` did not exist.

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
| S1 | ❌ Missing | new unit |
| S2 | ❌ Missing | new unit |
| S3 | ❌ Missing | new unit |
| S4 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/TranscriptomeAnalyzer_DifferentialExpression_Tests.cs` — all TRANS-DIFF-001 cases.
- **Remove:** nothing (no pre-existing tests for these methods).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `TranscriptomeAnalyzer_DifferentialExpression_Tests.cs` | canonical TRANS-DIFF-001 fixture | 13 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | implemented exact log2FC | ✅ Done |
| 2 | M2 | ❌ Missing | implemented sign symmetry | ✅ Done |
| 3 | M3 | ❌ Missing | implemented flat = 0 | ✅ Done |
| 4 | M4 | ❌ Missing | implemented Welch t/df/p | ✅ Done |
| 5 | M5 | ❌ Missing | implemented BH adjusted values | ✅ Done |
| 6 | M6 | ❌ Missing | implemented BH monotone/≥raw | ✅ Done |
| 7 | M7 | ❌ Missing | implemented two-criterion gate | ✅ Done |
| 8 | M8 | ❌ Missing | implemented regulation label | ✅ Done |
| 9 | S1 | ❌ Missing | implemented empty input | ✅ Done |
| 10 | S2 | ❌ Missing | implemented identical groups | ✅ Done |
| 11 | S3 | ❌ Missing | implemented <2 replicates | ✅ Done |
| 12 | S4 | ❌ Missing | implemented threshold gate | ✅ Done |
| 13 | C1 | ❌ Missing | implemented zero-mean pseudocount | ✅ Done |

**Total items:** 13
**✅ Done:** 13 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `CalculateFoldChange_TreatmentHigher_ReturnsExactLog2Ratio` |
| M2 | ✅ Covered | `CalculateFoldChange_ArgumentsSwapped_ReturnsNegatedValue` |
| M3 | ✅ Covered | `CalculateFoldChange_EqualMeans_ReturnsZero` |
| M4 | ✅ Covered | `FindDifferentiallyExpressed_WelchExample_ReturnsExactStatisticAndPValue` |
| M5 | ✅ Covered | `FindDifferentiallyExpressed_BenjaminiHochberg_ReturnsExactAdjustedPValues` |
| M6 | ✅ Covered | `FindDifferentiallyExpressed_AdjustedPValues_AreMonotoneAndGreaterEqualRaw` |
| M7 | ✅ Covered | `FindDifferentiallyExpressed_TwoCriterionGate_FlagsOnlyStrongGene` |
| M8 | ✅ Covered | `FindDifferentiallyExpressed_RegulationLabel_MatchesFoldChangeSign` |
| S1 | ✅ Covered | `FindDifferentiallyExpressed_EmptyInput_ReturnsEmpty` |
| S2 | ✅ Covered | `FindDifferentiallyExpressed_IdenticalGroups_NotSignificant` |
| S3 | ✅ Covered | `FindDifferentiallyExpressed_SingleReplicate_PValueIsOne` |
| S4 | ✅ Covered | `FindDifferentiallyExpressed_BelowFoldChangeThreshold_NotSignificant` |
| C1 | ✅ Covered | `CalculateFoldChange_ZeroControlMean_ReturnsFinitePositiveValue` |

---

## 6. Assumption Register

**Total assumptions:** 3

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Fold-change pseudocount = 1 (degenerate zero-mean regularization) | CalculateFoldChange, C1 |
| 2 | <2 replicates per group ⇒ p-value = 1 (Welch precondition undefined) | FindDifferentiallyExpressed, S3 |
| 3 | se = 0: equal means ⇒ p = 1, unequal means ⇒ p = 0 | FindDifferentiallyExpressed |

All three are degenerate-input conventions; none changes any non-degenerate (N≥2, finite-variance) output defined by the cited formulas.

---

## 7. Open Questions / Decisions

1. DESeq2 uses a shrunken-GLM/Wald estimator; this unit implements the classical replicate-level estimator (mean log-ratio + Welch t-test + BH), which is the textbook DE pipeline and the directly testable contract for the registry-named methods. The DESeq2 dispersion-shrinkage refinement is out of scope and recorded in the algorithm doc §5.3 "Not implemented".
