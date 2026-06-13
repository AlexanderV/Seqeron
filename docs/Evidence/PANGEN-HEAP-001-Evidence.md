# Evidence Artifact: PANGEN-HEAP-001

**Test Unit ID:** PANGEN-HEAP-001
**Algorithm:** Pan-Genome Growth Model (Heaps' law fit + gene presence/absence matrix)
**Date Collected:** 2026-06-13

---

## Online Sources

### Tettelin et al. (2005), PNAS — Streptococcus agalactiae pan-genome

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC1216834/
**Accessed:** 2026-06-13 (retrieved via WebFetch of the PMC full-text page; discovered via WebSearch query "Tettelin 2005 Streptococcus agalactiae pan-genome ... PNAS")
**Authority rank:** 1 (peer-reviewed primary paper, PNAS)

**Key Extracted Points:**

1. **Core-genome decay model:** `F_c = κc · exp[−n/τc] + Ω` with κc = 610 ± 38, τc = 2.16 ± 0.28, Ω = 1,806 ± 16 (extrapolated core size), r² = 0.990.
2. **New (strain-specific) genes model:** `F_s = κs · exp[−n/τs] + tg(θ)` with κs = 476 ± 62, τs = 1.51 ± 0.15, tg(θ) = 33 ± 3.5 (asymptotic new genes per added genome), r² = 0.995.
3. **Worked new-gene numbers:** the second genome added **161 new genes**; the fifth genome added **54 new genes**; the eighth genome still continued to add new genes.
4. **Core fraction:** the core genome is "≈80% of any single genome."
5. **Openness conclusion:** the predicted asymptotic number of new genes per genome is 33 (nonzero, probability < 6×10⁻⁴ that it equals zero) — i.e. an **open** pan-genome (the gene reservoir is vast; unique genes keep being found even after hundreds of genomes).

### Tettelin et al. (2008), Current Opinion in Microbiology — the bacterial pan-genome

**URL:** https://pubmed.ncbi.nlm.nih.gov/19086349/ (record); concepts confirmed via WebSearch query "Tettelin 2008 Comparative genomics bacterial pan-genome ... power law Heaps alpha open closed"
**Accessed:** 2026-06-13 (WebSearch result summary; the power-law form and open/closed rule are restated by the micropan reference implementation below, which cites this paper)
**Authority rank:** 1 (peer-reviewed review, Curr Opin Microbiol 11(5):472–477)

**Key Extracted Points:**

1. **Power-law growth:** the pan-genome new-gene-discovery curve is modelled by a **power law** rather than only exponential decay; openness is judged by the power-law decay exponent.
2. **Core vs dispensable:** pan-genome = core genome (shared by all strains) + dispensable genome (subset-shared + strain-specific genes).
3. **Open/closed framework:** open pan-genomes keep gaining new genes as strains are sequenced; closed pan-genomes have a bounded gene pool.

### micropan `heaps()` — reference implementation (R/powerlaw.R)

**URL:** https://raw.githubusercontent.com/larssnip/micropan/master/R/powerlaw.R (function source); docs https://search.r-project.org/CRAN/refmans/micropan/html/heaps.html
**Accessed:** 2026-06-13 (WebFetch of the raw GitHub source and the CRAN refman page; directory listing via GitHub API contents endpoint)
**Authority rank:** 3 (established bioinformatics reference implementation; implements Tettelin et al. 2008)

**Key Extracted Points:**

1. **Model:** `y.hat <- p[1] * x^(-p[2])` — number of new gene clusters `y` as a function of genome index `x`; `p[1]` = Intercept (K), `p[2]` = decay exponent alpha.
2. **Binarization:** `pan.matrix[which(pan.matrix > 0, arr.ind = T)] <- 1` (presence/absence is binary; counts > 0 mapped to 1).
3. **New-cluster counting per ordering:** `cm <- apply(pan.matrix[sample(nrow(pan.matrix)),], 2, cumsum)` then per added genome `rowSums((cm == 1)[2:ng,] & (cm == 0)[1:(ng-1),])` — a cluster is "new" at row i (i ≥ 2) when its cumulative count is 1 at row i and 0 at row i−1, i.e. it first appears at genome i.
4. **Pooling over permutations:** `x <- rep((2:nrow(pan.matrix)), times = n.perm)`, `y <- as.numeric(nmat)`; genome index starts at **2** (first genome has no predecessor).
5. **Objective (least squares):** `objectFun <- function(p, x, y){ y.hat <- p[1]*x^(-p[2]); J <- sqrt(sum((y - y.hat)^2))/length(x); return(J) }`.
6. **Optimization & bounds:** `optim(p0, objectFun, ..., method = "L-BFGS-B", lower = c(0, 0), upper = c(10000, 2))` — K ∈ [0, 10000], alpha ∈ [0, 2].
7. **Start values:** `p0 <- c(mean(y[which(x == 2)]), 1)` — intercept starts at the mean new-gene count at the second genome, alpha starts at 1.
8. **Return value:** named vector `c("Intercept", "alpha")`.
9. **Open/closed rule (verbatim):** "If 'alpha>1.0' the pan-genome is closed, if 'alpha<1.0' it is open."
10. **n.perm default:** 100 — "The default value of 100 is certainly a minimum."

---

## Documented Corner Cases and Failure Modes

### From micropan `heaps()` / Tettelin 2008

1. **Fewer than 2 genomes:** the new-gene curve `x = 2:ng` is empty when ng < 2 — no fit is defined.
2. **Genome index starts at 2:** the first genome contributes no "new" count (no predecessor); fitting only uses N = 2..G.
3. **Binary presence only:** copy-number > 1 is collapsed to 1 before counting (presence/absence, not abundance).
4. **alpha boundary:** alpha is constrained to [0, 2]; alpha = 1 is the open/closed boundary (alpha < 1 open, alpha > 1 closed); the boundary value itself is treated as not-open (closed) by the strict inequality.

---

## Test Datasets

### Dataset: Tettelin 2005 S. agalactiae new-gene observations (qualitative anchor)

**Source:** Tettelin et al. (2005), PNAS 102(39):13950–13955, https://pmc.ncbi.nlm.nih.gov/articles/PMC1216834/

| Parameter | Value |
|-----------|-------|
| New genes at genome 2 | 161 |
| New genes at genome 5 | 54 |
| Asymptotic new genes per genome tg(θ) | 33 ± 3.5 |
| Core fraction of a single genome | ≈ 80% |
| Pan-genome openness | Open (nonzero asymptote) |

### Dataset: Closed-form power-curve fit (exact, derived)

**Source:** Derived from the micropan model `y = K·x^(−alpha)` (powerlaw.R, point 1/5/6 above). Two points lying exactly on a single power curve determine (K, alpha) uniquely with objective J = 0.

| Genome index x | New gene clusters y |
|---|---|
| 2 | 8 |
| 3 | 4 |

Solving `K·2^(−alpha) = 8`, `K·3^(−alpha) = 4`: ratio `8/4 = (3/2)^alpha` → **alpha = ln(2)/ln(3/2) = 1.70951129135145** (within [0,2]); **K = 8·2^alpha = 26.1640013949735**. alpha > 1 → **closed**. This curve is realized by a fixed-order presence/absence matrix where genome 2 introduces 8 clusters absent in genome 1 and genome 3 introduces 4 clusters absent in genomes 1–2.

### Dataset: Constant new-gene curve (open boundary, exact)

**Source:** Derived from the same model. If every added genome introduces the same number of new clusters (y constant), the best power fit is alpha = 0, K = mean(y).

| Genome index x | New gene clusters y |
|---|---|
| 2 | 1 |
| 3 | 1 |

Best fit `K·x^(−alpha) ≡ 1` → **alpha = 0**, **K = 1**, J = 0. alpha < 1 → **open**.

---

## Assumptions

1. **ASSUMPTION: Bounded coordinate-descent optimizer vs. L-BFGS-B.** micropan calls R's `optim(method="L-BFGS-B")`. This implementation minimizes the identical objective `J = sqrt(Σ(y−K·x^(−alpha))²)/|x|` over the identical box constraints (K ∈ [0,10000], alpha ∈ [0,2]) from the identical start point, using deterministic coordinate descent with geometric step refinement. For data on an exact power curve within the bounds the global minimum is unique and both optimizers reach it; the recovered (K, alpha) is verified to match the analytic solution to < 1e-9 (see derived datasets). The optimization *method* is non-correctness-affecting (the objective, bounds, and start point — which determine the result — are copied verbatim from the source).
2. **ASSUMPTION: Permutation RNG.** micropan uses R's `sample()`; the per-genome ordering is random, so the pooled curve (and thus the fit) is permutation-dependent except for permutation-invariant matrices. This implementation uses a fixed seed and uses the natural input order for the first permutation, so single-permutation fixed-order fits are exactly reproducible. The averaging principle (pool over many orderings) matches the source; the specific RNG stream does not affect correctness for the deterministic test matrices used.

---

## Recommendations for Test Coverage

1. **MUST Test:** new-gene-count curve extraction equals the micropan `(cm==1 & prev cm==0)` first-appearance rule on a fixed-order presence/absence matrix — Evidence: micropan powerlaw.R points 2–4.
2. **MUST Test:** exact (K, alpha) recovery for the closed-form power curve x=[2,3], y=[8,4] → alpha = ln2/ln(3/2) ≈ 1.7095, K ≈ 26.1640, IsOpen = false — Evidence: micropan model + derived dataset.
3. **MUST Test:** constant new-gene curve → alpha = 0, K = 1, IsOpen = true (open) — Evidence: derived dataset + open/closed rule.
4. **MUST Test:** open/closed classification follows alpha < 1 / alpha > 1 strictly — Evidence: micropan rule (point 9).
5. **MUST Test:** fewer than 2 genomes (and empty/null input) → degenerate fit (Intercept = 0, predictor → 0), not an exception — Evidence: corner case 1.
6. **MUST Test:** binarization — a cluster present multiple times / duplicate gene ids counts once — Evidence: micropan binarization (point 2).
7. **SHOULD Test:** dictionary overload clusters then fits and agrees with the matrix overload — Rationale: convenience wrapper must delegate to the canonical matrix path.
8. **COULD Test:** predictor monotonic decreasing in N for alpha > 0 — Rationale: model property `N^(−alpha)` is decreasing.

---

## References

1. Tettelin H, Masignani V, Cieslewicz MJ, Donati C, Medini D, et al. (2005). Genome analysis of multiple pathogenic isolates of *Streptococcus agalactiae*: Implications for the microbial "pan-genome". *PNAS* 102(39):13950–13955. https://doi.org/10.1073/pnas.0506758102 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC1216834/)
2. Tettelin H, Riley D, Cattuto C, Medini D (2008). Comparative genomics: the bacterial pan-genome. *Current Opinion in Microbiology* 11(5):472–477. https://doi.org/10.1016/j.mib.2008.09.006 (record: https://pubmed.ncbi.nlm.nih.gov/19086349/)
3. Snipen L, Liland KH. micropan: Microbial Pan-Genome Analysis — `heaps()` (R/powerlaw.R). https://raw.githubusercontent.com/larssnip/micropan/master/R/powerlaw.R ; docs https://search.r-project.org/CRAN/refmans/micropan/html/heaps.html (accessed 2026-06-13).

---

## Change History

- **2026-06-13**: Initial documentation.
