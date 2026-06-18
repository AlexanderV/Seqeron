# Evidence Artifact: ONCO-EXPR-001

**Test Unit ID:** ONCO-EXPR-001
**Algorithm:** Tumor Gene Expression Outlier (z-score) and Signature Score
**Date Collected:** 2026-06-15

---

## Online Sources

### cBioPortal — mRNA expression z-score normalization specification

**URL:** https://docs.cbioportal.org/z-score-normalization-script/
**Accessed:** 2026-06-15
**Retrieved how:** WebSearch query `cBioPortal z-score threshold 2.0 mRNA expression altered overexpression underexpression default cutoff`, then WebFetch of the page above.
**Authority rank:** 5 (well-maintained bioinformatics database / platform specification) — formula corroborated by the reference implementation source (rank 3, below).

**Key Extracted Points:**

1. **Z-score formula (verbatim):** "`(r - mu)/sigma` where `r` is the raw expression value, and `mu` and `sigma` are the mean and standard deviation" of the reference base population.
2. **Reference (base) population:** Two methods. Method 1 (diploid): mean/sd from samples where the gene's copy-number value is 0 ("the expression distribution for unaltered copies of the gene"). Method 2 (all samples): mean/sd from "all samples with expression values."
3. **Zero-variance handling:** "Z-Score ← NA when standard deviation = 0"; diploid method reports NA when "the gene has no diploid samples."

### cBioPortal — NormalizeExpressionLevels.java (reference implementation)

**URL:** https://github.com/cBioPortal/cbioportal-core/blob/master/src/main/java/org/mskcc/cbio/portal/scripts/NormalizeExpressionLevels.java
**Accessed:** 2026-06-15
**Retrieved how:** `gh search code NormalizeExpressionLevels --repo cBioPortal/cbioportal-core`, then `gh api repos/cBioPortal/cbioportal-core/contents/.../NormalizeExpressionLevels.java` (base64-decoded) and inspected the `avg()` / `std()` / `getZ()` methods.

**Key Extracted Points:**

1. **Formula (verbatim comment):** `(r - mu)/sigma` … `zScore <- (value - mean)/sd`.
2. **Mean (`avg`):** `avg = (Σ vᵢ) / v.length` — arithmetic mean over the reference values.
3. **Standard deviation divisor (`std`):** `std = Σ(vᵢ-avg)² ; std = std/(double)(v.length-1) ; std = sqrt(std)` — i.e. **sample standard deviation, divisor (n−1)**. This authoritatively resolves the n-vs-(n−1) question that the prose spec leaves unstated.
4. **Zero standard deviation:** `getZ` calls `fatalError("cannot normalize relative to distribution with standard deviation of 0.0.")` (throws `RuntimeException`).

### cBioPortal — FAQ: default z-score threshold

**URL:** https://docs.cbioportal.org/user-guide/faq/
**Accessed:** 2026-06-15
**Retrieved how:** WebSearch (same query as source 1), then WebFetch of the FAQ page.

**Key Extracted Points:**

1. **Default outlier threshold (verbatim):** "By default, samples with expression z-scores >2 or <-2 in any queried genes are considered altered." z > 2 ⇒ overexpressed; z < −2 ⇒ underexpressed.
2. **Statistical basis:** a z-score of 2 = "a difference of two standard deviations from the mean"; ±2 ≈ 95th percentile of a normal distribution.

### Lee E., Chuang H-Y., Kim J-W., Ideker T., Lee D. (2008) — combined z-score (signature/pathway activity score)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC2563693/ (PLoS Comput Biol 4(11):e1000217)
**Accessed:** 2026-06-15
**Retrieved how:** WebSearch `combined z-score gene signature Lee 2008 inferring pathway activity average of z-scores divided by sqrt number of genes`, then WebFetch of the PMC article.

**Key Extracted Points:**

1. **Per-gene standardization (verbatim):** "expression values *gij* are normalized to z-transformed scores *zij* which for each gene *i* have mean μi = 0 and standard deviation σi = 1 over all samples *j*." (i.e. zᵢⱼ = (gᵢⱼ − μᵢ)/σᵢ.)
2. **Combined z-score / activity (verbatim):** "The individual *zij* of each member gene in the gene set are averaged into a combined *z*-score which is designated the activity *aj* (the square root of the number of member genes is used in the denominator to stabilize the variance of the mean)." ⇒ for a gene set of k members, activity **a = (Σᵢ zᵢ) / √k**.

### GSVA package vignette — corroboration of the combined-z-score definition

**URL:** https://bioconductor.org/packages/devel/bioc/vignettes/GSVA/inst/doc/GSVA.html
**Accessed:** 2026-06-15
**Retrieved how:** appeared in WebSearch results for the Lee-2008 query; used to corroborate that the "combined z-score" method (Lee et al. 2008) computes a per-gene z-score then sums over the gene set divided by √k.

**Key Extracted Points:**

1. **Corroboration:** GSVA lists the "combined z-score" (`zscore`) method attributed to Lee et al. (2008), standardizing each gene to a z-score and combining gene-set members by Σz/√k.

---

## Documented Corner Cases and Failure Modes

### From cBioPortal NormalizeExpressionLevels.java / spec

1. **Zero standard deviation reference:** the reference implementation aborts (`fatalError`) — a degenerate (constant) reference gene has no defined z-score. We surface this as a thrown exception.
2. **Empty / single-sample reference:** the (n−1) divisor is undefined for n ≤ 1 (a single reference value has no spread); the sample SD requires at least 2 reference samples.
3. **Threshold boundary:** the rule is strict `> 2` / `< -2`; a z of exactly ±2 is **not** an outlier.

### From Lee et al. (2008)

1. **Single-gene signature (k = 1):** a = z₁/√1 = z₁ (well defined).
2. **Empty signature (k = 0):** Σz/√0 is undefined — invalid input.

---

## Test Datasets

### Dataset: Hand-derived reference cohort (clean n−1 statistics)

**Source:** derived directly from the cBioPortal formula z = (r−μ)/σ with σ the sample SD (divisor n−1) per NormalizeExpressionLevels.java `std()`.

| Parameter | Value |
|-----------|-------|
| Reference cohort (gene G) | {2, 2, 4, 6, 6} |
| Mean μ | (2+2+4+6+6)/5 = 4 |
| Σ(rᵢ−μ)² | 4+4+0+4+4 = 16 |
| Sample variance | 16/(5−1) = 4 |
| Sample SD σ | √4 = 2 |
| z for x = 10 | (10−4)/2 = **3.0** (outlier, over; |z|>2) |
| z for x = 8 | (8−4)/2 = **2.0** (boundary, NOT outlier; rule is strict >2) |
| z for x = 4 | 0.0 |
| z for x = −1 | (−1−4)/2 = **−2.5** (outlier, under; z<−2) |

### Dataset: Hand-derived signature (combined z-score)

**Source:** Lee et al. (2008) a = (Σ zᵢ)/√k.

| Parameter | Value |
|-----------|-------|
| Member-gene z-scores | {3, 1, −1, 1} |
| k | 4 |
| Σ zᵢ | 4 |
| √k | 2 |
| Activity a | 4/2 = **2.0** |
| Single-gene set z = {2.5} | a = 2.5/√1 = **2.5** |

---

## Assumptions

1. **ASSUMPTION: Caller supplies the reference cohort and signature gene set.** Per the unit scope, reference distributions and signatures are caller-supplied (not fabricated). The algorithm computes z and a from caller-provided normalized expression; it does not bundle any specific cohort or signature. This is an API/scope decision, not a correctness-affecting numeric assumption.
2. **ASSUMPTION: Inputs are already on a normalization scale where a z-score is meaningful** (e.g. log-transformed TPM/FPKM/microarray intensity), matching cBioPortal's "raw expression value `r`" being whatever expression matrix is supplied. The z-score formula itself is scale-agnostic; this does not change any computed output for given inputs.

---

## Recommendations for Test Coverage

1. **MUST Test:** per-gene z-score equals (x−μ)/σ with σ the sample SD (n−1), on the hand-derived cohort {2,2,4,6,6} (μ=4, σ=2): x=10→3.0, x=4→0.0, x=−1→−2.5. — Evidence: cBioPortal formula + NormalizeExpressionLevels.java `std()`.
2. **MUST Test:** outlier classification uses strict thresholds >2 / <−2: x=10 over, x=−1 under, x=8 (z=2.0) and x=4 NOT outlier. — Evidence: cBioPortal FAQ "z-scores >2 or <-2 … considered altered".
3. **MUST Test:** combined signature score a = Σz/√k = 2.0 for z={3,1,−1,1}; a=2.5 for single-gene {2.5}. — Evidence: Lee et al. (2008).
4. **MUST Test:** zero-SD reference (constant cohort) throws. — Evidence: NormalizeExpressionLevels.java `fatalError`.
5. **SHOULD Test:** null/empty cohort, reference of size 1 (SD undefined), empty signature — argument validation. — Rationale: (n−1) divisor and √k undefined for these inputs.
6. **COULD Test:** symmetry/monotonicity property (z increases with x for fixed cohort; sign flips). — Rationale: invariant check.

---

## References

1. cBioPortal. 2024. mRNA Expression Z-Scores normalization. cBioPortal documentation. https://docs.cbioportal.org/z-score-normalization-script/
2. cBioPortal core. NormalizeExpressionLevels.java (reference implementation; `avg`, `std` divisor n−1, `getZ`). https://github.com/cBioPortal/cbioportal-core/blob/master/src/main/java/org/mskcc/cbio/portal/scripts/NormalizeExpressionLevels.java
3. cBioPortal. FAQ — default z-score threshold ±2 for altered mRNA expression. https://docs.cbioportal.org/user-guide/faq/
4. Lee E, Chuang H-Y, Kim J-W, Ideker T, Lee D. 2008. Inferring Pathway Activity toward Precise Disease Classification. PLoS Comput Biol 4(11):e1000217. https://doi.org/10.1371/journal.pcbi.1000217 (PMC2563693)
5. Hänzelmann S, Castelo R, Guinney J. 2013. GSVA: gene set variation analysis (combined z-score method, Lee et al. 2008). Bioconductor vignette. https://bioconductor.org/packages/devel/bioc/vignettes/GSVA/inst/doc/GSVA.html

---

## Change History

- **2026-06-15**: Initial documentation.
