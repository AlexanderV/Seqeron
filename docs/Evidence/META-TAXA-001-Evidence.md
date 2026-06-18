# Evidence Artifact: META-TAXA-001

**Test Unit ID:** META-TAXA-001
**Algorithm:** Significant Taxa Detection (Wilcoxon rank-sum / Mann–Whitney U test, normal approximation)
**Date Collected:** 2026-06-13

---

## Online Sources

### Mann–Whitney U test — Wikipedia (cites primary Mann & Whitney 1947)

**URL:** https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test
**Accessed:** 2026-06-13 (retrieved via WebFetch on the URL above)
**Authority rank:** 4 (Wikipedia citing primary; primary used = Mann & Whitney 1947, rank 1)

**Key Extracted Points:**

1. **U statistic definition:** `U1 = R1 − n1(n1+1)/2` and equivalently `U1 = n1·n2 + n1(n1+1)/2 − R1`, where `R1` is the sum of ranks assigned to sample 1 in the pooled ranking.
2. **Complementarity:** `U1 + U2 = n1·n2`.
3. **Normal approximation mean:** `m_U = n1·n2 / 2`.
4. **Normal approximation standard deviation:** `σ_U = sqrt( n1·n2·(n1+n2+1) / 12 )`.
5. **Z-score:** `z = (U − m_U) / σ_U`.
6. **Tie correction:** `σ_ties = sqrt( n1·n2·(n1+n2+1)/12 − n1·n2·Σ(t_k³ − t_k) / (12·n·(n−1)) )`, where `n = n1+n2` and `t_k` is the number of tied observations at distinct value `k`.
7. **Tied ranks:** tied observations receive the midpoint (average) of the ranks they would otherwise occupy (worked illustration: unadjusted ranks `(1,2,3,4,5,6)` for four tied middle values become `(1, 3.5, 3.5, 3.5, 3.5, 6)`).
8. **Worked example (tortoise & hare):** 6 tortoises vs 6 hares; pooled finishing ranks (12 = first). Tortoise rank sum `R_T = 12+6+5+4+3+2 = 32`; `U_T = R_T − n(n+1)/2 = 32 − (6·7)/2 = 32 − 21 = 11`; pairwise-comparison count gives `U_T = 11`, `U_H = 25`, and `U_T + U_H = 36 = 6·6`.
9. **Primary citation:** Mann, Henry B.; Whitney, Donald R. (1947). "On a Test of Whether one of Two Random Variables is Stochastically Larger than the Other." *Annals of Mathematical Statistics* 18(1):50–60.

### SciPy `scipy.stats.mannwhitneyu` — reference implementation documentation

**URL:** https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html
**Accessed:** 2026-06-13 (retrieved via WebFetch on the URL above)
**Authority rank:** 3 (reference implementation in an established scientific library, SciPy)

**Key Extracted Points:**

1. **Statistic complement:** "If `U1` is the statistic corresponding with sample x, then the statistic corresponding with sample y is `U2 = x.shape[axis] * y.shape[axis] - U1`."
2. **Asymptotic mean / sd:** mean `μ = n_x·n_y / 2`; standard deviation computed as `sqrt(nx·ny·(N+1)/12)` with `N = nx + ny` (with tie correction when ties present).
3. **Continuity correction:** "the continuity correction performed by `mannwhitneyu` reduces the distance between the test statistic and the mean `μ = n_x n_y / 2` by 0.5." Default `use_continuity = True` when `method='asymptotic'`.
4. **Worked example (reference output):** x = [19, 22, 16, 29, 24] (n1=5), y = [20, 11, 17, 12] (n2=4) → `U1 = 17.0`, `U2 = 3.0`; asymptotic p-value with continuity correction = `0.11134688653314041`; asymptotic p-value without continuity = `0.0864107329737`.

### Xia Y, Sun J (2017) — Hypothesis testing and statistical analysis of microbiome (PMC6128532)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC6128532/
**Accessed:** 2026-06-13 (retrieved via WebFetch on the URL above)
**Authority rank:** 1 (peer-reviewed review, *Genes & Diseases*)

**Key Extracted Points:**

1. **Domain applicability:** "The non-parametric analogous Wilcoxon rank-sum test (also called Mann–Whitney test) ... was used to identify the statistically significant differences in microbial taxa or OTUs."
2. **Role:** "two-sample t-test and its nonparametric counterpart Wilcoxon rank-sum test were widely used in microbiome studies to comparing continuous variables between two groups." → justifies using the rank-sum test per-taxon across two sample groups for significant-taxa detection.

### Abramowitz & Stegun 7.1.26 erf approximation (used by repository `NormalCDF`)

**URL:** https://www.johndcook.com/blog/python_erf/
**Accessed:** 2026-06-13 (retrieved via WebFetch on the URL above)
**Authority rank:** 2 (documents official A&S Handbook formula 7.1.26)

**Key Extracted Points:**

1. **Formula source:** the error-function approximation `y = 1 − ((((a5·t + a4)·t + a3)·t + a2)·t + a1)·t·exp(−x²)`, `t = 1/(1 + p·x)`, is "A&S formula 7.1.26" (Abramowitz & Stegun, Handbook of Mathematical Functions).
2. **Constants:** `a1 = 0.254829592, a2 = −0.284496736, a3 = 1.421413741, a4 = −1.453152027, a5 = 1.061405429, p = 0.3275911` — identical to those in the repository's `StatisticsHelper.Erf`.
3. **Accuracy:** formula 7.1.26 is a finite-precision polynomial approximation (A&S states maximum error |ε| ≤ 1.5×10⁻⁷); consequently p-values derived through `NormalCDF` match exact normal-CDF p-values only to ≈1×10⁻⁶.

---

## Documented Corner Cases and Failure Modes

### From SciPy mannwhitneyu / Mann–Whitney definition

1. **Ties:** present in real abundance data; require midrank assignment and the tie-corrected σ. Without correction σ is overstated and p-values are conservative.
2. **Identical groups (zero variance / all ties):** when every observation is tied across both groups, the corrected σ → 0 and `z` is undefined; the test cannot reject H₀ → report p = 1.
3. **Small samples:** the normal approximation is asymptotic; for very small n it is approximate (SciPy offers an exact method). This implementation uses the asymptotic approximation per the documented σ/μ formulas.

### From Xia & Sun (2017)

4. **Per-taxon multiple testing:** each taxon is tested independently; raw p-values are not corrected for multiplicity here (correction is the caller's responsibility).

---

## Test Datasets

### Dataset: SciPy mannwhitneyu reference example

**Source:** SciPy `scipy.stats.mannwhitneyu` documentation (reference output, accessed 2026-06-13)

| Parameter | Value |
|-----------|-------|
| group1 (x) | 19, 22, 16, 29, 24 |
| group2 (y) | 20, 11, 17, 12 |
| n1, n2 | 5, 4 |
| U1 (group1) | 17.0 |
| U2 (group2) | 3.0 |
| m_U = n1·n2/2 | 10 |
| σ_U = sqrt(n1·n2·(n1+n2+1)/12) | sqrt(200/12) = 4.08248290463863 |
| z (no continuity, using max U) | (17−10)/4.08248290463863 = 1.7146428199482247 |
| z (continuity, max U) | (17−10−0.5)/4.08248290463863 = 1.5921683328090657 |
| two-tailed p (SciPy, with continuity) | 0.11134688653314041 |
| two-tailed p (SciPy, no continuity) | 0.0864107329737 |

### Dataset: Mann–Whitney tortoise & hare example

**Source:** Wikipedia "Mann–Whitney U test" (accessed 2026-06-13), citing Mann & Whitney (1947)

| Parameter | Value |
|-----------|-------|
| n1 (tortoises), n2 (hares) | 6, 6 |
| tortoise rank sum R_T | 32 |
| U_T = R_T − n1(n1+1)/2 | 32 − 21 = 11 |
| U_H | 25 |
| U_T + U_H | 36 (= n1·n2) |
| σ_U = sqrt(6·6·13/12) | sqrt(39) = 6.244997998398398 |

---

## Assumptions

1. **ASSUMPTION: Continuity correction enabled by default** — The implementation applies the 0.5 continuity correction by default, matching SciPy's documented default (`use_continuity=True` for the asymptotic method). This is a documented reference-implementation default, not invented; exposed as an optional parameter so the un-corrected variant remains reachable.
2. **ASSUMPTION: Two-tailed alternative** — Significance is judged with a two-tailed test (`p = 2·SF(|z|)`), the default for "is taxon differentially abundant" questions. SciPy's example p-values quoted above are two-sided.
3. **ASSUMPTION: Group label set** — Each profile is mapped to one of exactly two groups; taxa absent from a profile contribute abundance 0 to that profile's vector (standard for abundance tables where a taxon may be unobserved in a sample).

---

## Recommendations for Test Coverage

1. **MUST Test:** SciPy reference example reproduces U1=17, U2=3, σ=sqrt(200/12), z (no-cc)=1.71464282, z (cc)=1.59216833, and two-tailed p within Erf tolerance — Evidence: SciPy mannwhitneyu documentation.
2. **MUST Test:** Tortoise/hare reproduces U_T=11, U_H=25, U_T+U_H=36 — Evidence: Wikipedia/Mann & Whitney (1947).
3. **MUST Test:** `U1 + U2 = n1·n2` invariant holds for arbitrary inputs — Evidence: Wikipedia/Mann & Whitney (1947).
4. **MUST Test:** ties receive midranks and tie-corrected σ is used — Evidence: Wikipedia tie-correction formula.
5. **MUST Test:** identical groups (all tied) → p = 1, not significant — Evidence: degenerate σ → 0 case.
6. **SHOULD Test:** significance flag = (p < threshold); a clearly separated pair is flagged significant, an overlapping pair is not — Rationale: contract of the public method.
7. **SHOULD Test:** null/empty inputs and single-group inputs handled per validation contract — Rationale: documented failure modes.
8. **COULD Test:** continuity correction toggle changes p in the documented direction (cc p > no-cc p) — Rationale: SciPy parameter semantics.

---

## References

1. Mann, H.B.; Whitney, D.R. (1947). On a Test of Whether one of Two Random Variables is Stochastically Larger than the Other. *Annals of Mathematical Statistics* 18(1):50–60. https://doi.org/10.1214/aoms/1177730491 (definition/formulas accessed via https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test , 2026-06-13)
2. SciPy developers. scipy.stats.mannwhitneyu — SciPy Manual. https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html (accessed 2026-06-13)
3. Xia, Y.; Sun, J. (2017). Hypothesis testing and statistical analysis of microbiome. *Genes & Diseases* 4(3):138–148. https://pmc.ncbi.nlm.nih.gov/articles/PMC6128532/ (accessed 2026-06-13)
4. Abramowitz, M.; Stegun, I.A. (1964). Handbook of Mathematical Functions, formula 7.1.26 (erf approximation). Documented at https://www.johndcook.com/blog/python_erf/ (accessed 2026-06-13)

---

## Change History

- **2026-06-13**: Initial documentation.
