# Evidence Artifact: TRANS-DIFF-001

**Test Unit ID:** TRANS-DIFF-001
**Algorithm:** Differential Expression (log2 fold change, Welch's t-test, Benjamini-Hochberg FDR)
**Date Collected:** 2026-06-13

---

## Online Sources

### Love, Huber & Anders (2014) — DESeq2 (Genome Biology, PMC4302049)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC4302049/
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed paper)
**Retrieved how:** WebFetch of the PMC article; prompt requested the verbatim definition of log2 fold change between two conditions, the sign convention, and the multiple-testing procedure used.

**Key Extracted Points:**

1. **LFC is the log2 ratio between treatment and control:** the GLM "fit returns coefficients indicating the overall expression strength of the gene and the log₂ fold change between treatment and control" (verbatim). The coefficient sign therefore points "treatment vs control".
2. **Reporting scale:** "the DESeq2 software reports estimated model coefficients and their estimated standard errors on the log₂ scale" (verbatim) — fold change is reported as log2.
3. **Multiple-testing correction:** "The Wald test P values ... are adjusted for multiple testing using the procedure of Benjamini and Hochberg" (verbatim) — establishes BH as the standard FDR procedure for RNA-seq DE.

### Science Park Study Group — "Differential expression analysis" RNA-seq lesson

**URL:** https://scienceparkstudygroup.github.io/rna-seq-lesson/06-differential-analysis/index.html
**Accessed:** 2026-06-13
**Authority rank:** 3 (curated teaching reference corroborating the primary fold-change definition and the two-criterion DE rule)
**Retrieved how:** WebFetch; prompt requested the fold-change definition, sign convention, and how a gene is declared DE using a fold-change threshold and an adjusted-p-value cutoff.

**Key Extracted Points:**

1. **Definition:** log2 fold change = `log₂(condition A / condition B)`; the numerator is the treatment condition (example given: `log2(Pseudomonas syringae DC3000 / mock)`).
2. **Sign convention:** "a log2 equal to 1 means that gene X has a higher expression (x2, two-fold)" in the treated state; positive = upregulated in the numerator (treatment), negative = downregulated.
3. **DE decision = two simultaneous criteria:** (1) an absolute log2 fold-change cutoff (commonly |1| or |2|) AND (2) an adjusted p-value below a selected alpha (typically 0.01 or 0.001). FDR adjustment is essential when running thousands of simultaneous tests.

### Welch's t-test (Wikipedia, citing Welch 1947)

**URL:** https://en.wikipedia.org/wiki/Welch%27s_t-test
**Accessed:** 2026-06-13
**Authority rank:** 4 (Wikipedia citing the primary B. L. Welch 1947 paper; formulas only)
**Retrieved how:** WebFetch; prompt requested the verbatim t-statistic, the Welch-Satterthwaite degrees-of-freedom formula, and the use of unbiased (N-1) sample variances.

**Key Extracted Points:**

1. **t-statistic (verbatim):** `t = (X̄₁ - X̄₂) / √(s²₁/N₁ + s²₂/N₂)`, where `sᵢ` is the corrected (N-1) sample standard deviation and `Nᵢ` the sample size of group i.
2. **Welch-Satterthwaite df (verbatim):** `ν ≈ (s²₁/N₁ + s²₂/N₂)² / [ s⁴₁/(N²₁(N₁-1)) + s⁴₂/(N²₂(N₂-1)) ]`, with `ν₁=N₁-1`, `ν₂=N₂-1`.
3. **No pooled variance:** the denominator "is not based on a pooled variance estimate", so it handles unequal variances/sizes.

### Student's t-distribution CDF via regularized incomplete beta (Wikipedia)

**URL:** https://en.wikipedia.org/wiki/Student%27s_t-distribution
**Accessed:** 2026-06-13
**Authority rank:** 4 (Wikipedia; standard special-function identity)
**Retrieved how:** WebFetch; prompt requested the t-distribution CDF and two-sided p-value in terms of the regularized incomplete beta function `I_x(a,b)`.

**Key Extracted Points:**

1. **CDF (verbatim):** for `t > 0`, `F(t) = 1 − ½·I_{x(t)}(ν/2, ½)` with `x(t) = ν/(t² + ν)`.
2. **Two-sided tail (verbatim):** `A(t|ν) = 1 − I_{ν/(ν+t²)}(ν/2, ½)` is the probability the statistic falls within `[−t, t]`; therefore the two-sided p-value `P(|T| ≥ t) = I_{ν/(ν+t²)}(ν/2, ½)`.

### Benjamini & Hochberg (1995) procedure — Wikipedia "False discovery rate" + R `p.adjust`

**URL:** https://en.wikipedia.org/wiki/False_discovery_rate ; https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html
**Accessed:** 2026-06-13
**Authority rank:** 2–4 (official-spec-level reference implementation R/stats + Wikipedia citing the primary BH 1995 paper)
**Retrieved how:** WebFetch of the FDR article (procedure statement + citation) and the R `p.adjust` manual (BH reference); WebSearch (`R p.adjust BH source code cummin ...`) returned the verbatim R BH expression from R-help/biostars.

**Key Extracted Points:**

1. **BH step-up procedure (verbatim):** with p-values ordered ascending `P(1) ≤ … ≤ P(m)`, "find the largest k for which P(k) ≤ (k/m)α. Reject the null hypothesis ... for all H(i) for i = 1, …, k."
2. **Citation (verbatim):** "Benjamini Y, Hochberg Y (1995). 'Controlling the false discovery rate: a practical and powerful approach to multiple testing.' Journal of the Royal Statistical Society, Series B. 57(1): 289–300."
3. **R reference-implementation BH adjusted p-value (verbatim code):**
   ```r
   i <- lp:1L
   o <- order(p, decreasing = TRUE)
   ro <- order(o)
   pmin(1, cummin(n/i * p[o]))[ro]
   ```
   i.e. sort p descending, multiply each by `n/rank` (rank = m..1), take running minimum from largest to smallest p, clamp to 1, restore original order. R `p.adjust` cites Benjamini & Hochberg 1995 for `method="BH"`.

---

## Documented Corner Cases and Failure Modes

### From Science Park / DESeq2

1. **Two-criterion gate:** a gene is significant only when BOTH |log2FC| ≥ threshold AND adjusted p-value < alpha; failing either is "not significant".
2. **Many simultaneous tests:** raw p-values must be FDR-adjusted (BH); reporting raw p-values across thousands of genes inflates false positives.

### From Welch's t-test

1. **Fewer than 2 replicates per group:** sample variance (divide by N-1) is undefined for N<2; the t-statistic cannot be formed. Resolution adopted (see Assumptions): such a gene cannot be tested → p-value = 1 (not significant), matching the existing analyzer convention.
2. **Zero pooled standard error (se = 0):** both groups constant. If the two means are equal → no difference → p = 1; if means differ with zero variance → t = ∞ → p = 0. (Degenerate; documented convention.)

### From fold change

1. **Zero mean expression:** `log₂(mean2/mean1)` is undefined when a mean is 0. A pseudocount is added to both means before the ratio (standard practice) so the ratio is finite. Resolution adopted: pseudocount = 1 (see Assumptions).

---

## Test Datasets

### Dataset: log2 fold change two-condition derivation

**Source:** Definition `log₂(mean_treatment / mean_control)` (DESeq2 Love et al. 2014; Science Park RNA-seq lesson). Pseudocount c=1 added to each mean (degenerate-input convention).

| Gene | control means | treatment means | mean1 | mean2 | log2FC = log₂((mean2+1)/(mean1+1)) |
|------|---------------|-----------------|-------|-------|-----|
| UP   | 10,10,10 | 40,40,40 | 10 | 40 | log₂(41/11) = 1.8981204… |
| DOWN | 40,40,40 | 10,10,10 | 40 | 10 | log₂(11/41) = −1.8981204… |
| FLAT | 20,20,20 | 20,20,20 | 20 | 20 | log₂(21/21) = 0 |

Derivation: `log₂(41/11) = log₂(3.7272727…) = 1.8981204…`; the DOWN gene is the exact negative (sign convention: positive = up in treatment).

### Dataset: Welch's t-test statistic derivation

**Source:** Welch t-statistic `t = (X̄₂−X̄₁)/√(s²₁/N₁+s²₂/N₂)` (Welch 1947 via Wikipedia). Two-sided p = `I_{ν/(ν+t²)}(ν/2,½)` (Student t CDF identity).

| Quantity | control = {1,2,3}, treatment = {7,8,9} |
|----------|-----------------------------------------|
| X̄₁, X̄₂ | 2, 8 |
| s²₁, s²₂ | 1, 1 (each: Σ(x−x̄)²/(N−1) = 2/2 = 1) |
| se = √(1/3+1/3) | √(2/3) = 0.8164966 |
| t = (8−2)/0.8164966 | 7.348469 |
| Welch df ν = (1/3+1/3)² / (2·[1/(9·2)]) | (4/9)/(1/9) = 4 |
| two-sided p = I_{4/(4+t²)}(2, 0.5) | 0.0018262607 (cross-checked vs SciPy `ttest_ind(equal_var=False)`) |

The exact df here is ν=4 (equal variances and sizes). The two-sided p-value 0.0018262607 was cross-checked against the SciPy reference implementation `scipy.stats.ttest_ind([7,8,9],[1,2,3],equal_var=False)`, which returns t=7.3484692283 and p=0.0018262607. The MUST assertion checks the exact t and df, p to 1e-6, and that p < alpha.

### Dataset: Benjamini-Hochberg adjusted p-values (derived from the R reference algorithm)

**Source:** R `p.adjust(method="BH")` algorithm `pmin(1, cummin(n/i * p[o]))[ro]` (cites Benjamini & Hochberg 1995). m=4 raw p-values.

| rank i (ascending) | raw p(i) | n/i · p(i) = 4/i · p(i) | cumulative min (from largest rank down) → adjusted |
|---|---|---|---|
| 1 | 0.01 | 4/1·0.01 = 0.04 | 0.04 |
| 2 | 0.02 | 4/2·0.02 = 0.04 | 0.04 |
| 3 | 0.03 | 4/3·0.03 = 0.04 | 0.04 |
| 4 | 0.04 | 4/4·0.04 = 0.04 | 0.04 |

Derivation (R order: descending p, multiply by n/i for i=m..1, cummin, clamp): products are (0.04,0.04,0.04,0.04); running minimum keeps 0.04; all clamp < 1. Hence every adjusted p = 0.04. A second example with separated p-values: raw (0.001, 0.4, 0.5, 0.9) → products from largest down: 4/4·0.9=0.9, 4/3·0.5=0.6667, 4/2·0.4=0.8→cummin keeps 0.6667, 4/1·0.001=0.004; cummin sequence (0.9, 0.6667, 0.6667, 0.004) → restore ascending order → adjusted (0.004, 0.6667, 0.6667, 0.9). Monotone non-decreasing in p-order.

---

## Assumptions

1. **ASSUMPTION: Fold-change pseudocount = 1.** `log₂(mean2/mean1)` is undefined when a mean is 0. No retrieved source mandates a specific pseudocount value for the simple ratio estimator; adding 1 to both means is the standard regularization and the value already used by the existing analyzer. It only affects genes with a near-zero mean and never changes the sign or the relative ordering of non-degenerate ratios. The unregularized definition `log₂(mean2/mean1)` (DESeq2/Science Park) is recovered exactly as means grow large.
2. **ASSUMPTION: Fewer than 2 replicates per group → p-value = 1.** Welch's t-statistic requires the unbiased sample variance (N−1 denominator), undefined for N<2. No source specifies behavior for a single replicate; emitting p=1 (gene not testable, treated as not significant) is the conservative convention and matches the existing analyzer. It never affects a properly replicated (N≥2) input.
3. **ASSUMPTION: se = 0 with equal means → p = 1; with unequal means → p = 0.** The t-statistic is 0/0 (p=1, no evidence of difference) or ±∞ (p=0, perfect separation with zero variance). This degenerate convention follows directly from the limit of the t-statistic and matches the existing analyzer.

---

## Recommendations for Test Coverage

1. **MUST Test:** `CalculateFoldChange` returns log₂((mean2+1)/(mean1+1)); UP gene = +1.8981204…, DOWN = −1.8981204… (exact negative), FLAT = 0 — Evidence: DESeq2 (Love 2014); Science Park lesson.
2. **MUST Test:** sign convention — positive log2FC ⇔ higher in treatment (condition 2) — Evidence: Science Park lesson.
3. **MUST Test:** `FindDifferentiallyExpressed` Welch t-statistic and df for {1,2,3} vs {7,8,9}: t = 7.348469, ν = 4, two-sided p < alpha — Evidence: Welch 1947; Student t CDF identity.
4. **MUST Test:** BH adjusted p-values reproduce the R `p.adjust` algorithm: raw (0.001,0.4,0.5,0.9) → adjusted (0.004,0.6667,0.6667,0.9), monotone non-decreasing — Evidence: R p.adjust / Benjamini-Hochberg 1995.
5. **MUST Test:** DE decision requires BOTH |log2FC| ≥ threshold AND adjusted p < alpha (a gene failing either is not significant) — Evidence: Science Park lesson; DESeq2.
6. **SHOULD Test:** empty input → empty output; identical groups → log2FC 0, p = 1, not significant.
7. **SHOULD Test:** group with <2 replicates → p = 1, not significant (failure mode).
8. **COULD Test:** adjusted p-value is always ≥ raw p-value and ≤ 1 (BH invariant).

---

## References

1. Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15:550. https://pmc.ncbi.nlm.nih.gov/articles/PMC4302049/
2. Benjamini Y, Hochberg Y. 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B 57(1):289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x (procedure/citation retrieved via https://en.wikipedia.org/wiki/False_discovery_rate and R p.adjust manual https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)
3. Welch BL. 1947. The generalization of "Student's" problem when several different population variances are involved. Biometrika 34(1–2):28–35. (formulas via) https://en.wikipedia.org/wiki/Welch%27s_t-test
4. Student's t-distribution — cumulative distribution function via regularized incomplete beta function. https://en.wikipedia.org/wiki/Student%27s_t-distribution
5. Science Park Study Group. Introduction to RNA-seq — 06 Differential expression analysis. https://scienceparkstudygroup.github.io/rna-seq-lesson/06-differential-analysis/index.html

---

## Change History

- **2026-06-13**: Initial documentation.
