# Evidence Artifact: ONCO-CCF-001

**Test Unit ID:** ONCO-CCF-001
**Algorithm:** Cancer Cell Fraction (CCF) point estimation and 1D CCF clustering into clones/subclones
**Date Collected:** 2026-06-15

---

## Online Sources

### A practical guide to cancer subclonal reconstruction from DNA sequencing (Tarabichi et al. 2021, *Nature Methods*)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7867630/
**Accessed:** 2026-06-15 (retrieved via WebFetch of the PMC full text after Google/WebSearch query "A practical guide to cancer subclonal reconstruction from DNA sequencing"; redirected from ncbi.nlm.nih.gov/pmc → pmc.ncbi.nlm.nih.gov)
**Authority rank:** 1 (peer-reviewed *Nature Methods* review)

**Key Extracted Points:**

1. **CCF formula (Box 1, verbatim):** "the fraction of cancer cells from the sequenced sample carrying a set of SNVs, *i.e.* CCF = CP / purity. It can be inferred from the VAF (*f*), given a sample purity (*ρ*), local copy-number (*N_T*) and the inferred multiplicity of the mutations m: CCF = f·(ρ·N_T + 2(1−ρ)) / (ρ·m)." (The flattened HTML renders the fraction as "fm ρ(ρN_T+2(1−ρ))"; the multiplicity m divides — see corroboration in the two sources below, which both place m and ρ in the denominator.)
2. **Multiplicity (Box 1, verbatim):** "the number of DNA copies bearing a mutation *m*, which can be estimated from the VAF *f*, sample purity *ρ* and total copy number of the region in the tumor cells (*N_T*) as m = f·(ρ·N_T + 2(1−ρ)) / ρ." For clonal copy-number regions the result is rounded to the nearest non-zero integer.
3. **VAF (Box 1, verbatim):** "the fraction of mutated reads for a given variant, which is a readout for the proportion of DNA mutated in the sequenced tissue."
4. **Clonal cluster rule (SNV clustering section, verbatim):** "The cluster with the highest CP can be deemed clonal, and the remaining clusters can be assigned CPs and associated with a subclonal lineage." (CP = cellular prevalence; CCF = CP/purity, so the cluster with the largest CCF is the clonal cluster.)

### Estimation of cancer cell fractions and clone trees from multi-region sequencing (Zheng et al. 2022, *Bioinformatics* — PICTograph)

**URL:** https://academic.oup.com/bioinformatics/article/38/15/3677/6596597
**Accessed:** 2026-06-15 (retrieved via WebFetch after WebSearch query for the CCF formula)
**Authority rank:** 1 (peer-reviewed *Bioinformatics*)

**Key Extracted Points:**

1. **VAF↔CCF relation (Materials & Methods, verbatim):** "VAF = (m·CCF·p) / (c·p + 2·(1 − p))", where m is multiplicity, CCF is cancer cell fraction, p is purity, c is tumor copy number, and the normal contributes 2·(1−p). Solving for CCF gives CCF = VAF·(c·p + 2(1−p)) / (m·p) — identical to the Box 1 formula above and to McGranahan 2016.
2. **No explicit [0,1] cap** is stated in the generative model; CCF is treated as a continuous fraction.

### Clonal neoantigens elicit T cell immunoreactivity… (McGranahan et al. 2016, *Science* 351:1463–1469)

**URL:** https://www.science.org/cms/asset/fc2a5f47-6e70-4dff-94df-787c854dcd03/pap.pdf (full text 403/paywalled; formula extracted from the indexed WebSearch result for query "McGranahan 2015 … CCF formula observed mutation copy number VAF purity")
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed *Science*)

**Key Extracted Points:**

1. **Observed mutation copy number (verbatim from indexed text):** "n_mut = VAF × (1/p) × [p×CN_t + CN_n×(1−p)]", where VAF is the variant allele frequency, p the tumor purity, CN_t the tumor locus-specific copy number, and CN_n the normal locus-specific copy number. CCF = n_mut / multiplicity, so CCF = VAF·[p·CN_t + CN_n·(1−p)] / (p·m) with CN_n = 2.

### CNAqc — Computation of Cancer Cell Fractions (Caravagna lab, CNAqc vignette)

**URL:** https://caravagnalab.github.io/CNAqc/articles/a4_ccf_computation.html
**Accessed:** 2026-06-15 (retrieved via WebFetch)
**Authority rank:** 3 (reference implementation / tool documentation)

**Key Extracted Points:**

1. **Worked outputs (verbatim):** "VAF = 0.883, mutation_multiplicity = 2, CCF = 0.993" and "VAF = 0.471, mutation_multiplicity = 1, CCF = 1.06" — demonstrates that the raw formula can yield CCF slightly above 1 due to sampling noise, i.e. the formula does not intrinsically bound CCF ≤ 1.

### K-means clustering — Lloyd's algorithm (Wikipedia, citing Lloyd 1982)

**URL:** https://en.wikipedia.org/wiki/K-means_clustering
**Accessed:** 2026-06-15 (retrieved via WebFetch; primary citation Lloyd, S.P. 1982, *IEEE Trans. Inf. Theory* 28(2):129–137)
**Authority rank:** 4 (Wikipedia, used for the cited primary Lloyd 1982)

**Key Extracted Points:**

1. **Assignment step (verbatim):** "Assign each observation to the cluster with the nearest mean (centroid): that with the least squared Euclidean distance."
2. **Update step (verbatim):** "Recalculate means (centroids) for observations assigned to each cluster."
3. **Objective (verbatim):** minimize the within-cluster sum of squares (WCSS), Σ_{i=1..k} Σ_{x∈S_i} ‖x − μ_i‖².

---

## Documented Corner Cases and Failure Modes

### From Tarabichi et al. 2021 (PMC7867630) / CNAqc

1. **CCF > 1 from noise:** the raw formula yields CCF slightly above 1 (CNAqc: 1.06) when VAF is sampled above its expectation; the registry invariant 0 ≤ CCF ≤ 1 is enforced by capping the *reported* CCF at 1 (a mutation present in all cancer cells has CCF = 1, per the McGranahan clonal definition).
2. **Multi-copy loci (ambiguous CCF):** when N_T > 2 the multiplicity m is not 1; an integer multiplicity (≥1, ≤ tumor copy number) must be supplied/estimated, otherwise CCF is ambiguous. The guide notes multiplicity is rounded to the nearest non-zero integer for clonal CN regions.
3. **Unknown purity:** purity ρ appears in the denominator; ρ must be in (0, 1].

---

## Test Datasets

### Dataset: CNAqc worked CCF outputs

**Source:** CNAqc vignette "Computation of Cancer Cell Fractions" (Caravagna lab).

| VAF | multiplicity | (purity, N_T) consistent values | CCF |
|-----|--------------|---------------------------------|-----|
| 0.471 | 1 | ρ=1.0, N_T=2 → 0.471·2/1 = 0.942 (vignette ρ≈0.89, N_T=4 → 1.06) | ~1.06 (uncapped) |

### Dataset: Derived exact CCF cases (from the McGranahan/PMC/Zheng formula CCF = VAF·(ρ·N_T + 2(1−ρ))/(ρ·m))

**Source:** Derivation from the verbatim formula above.

| Case | VAF (f) | purity (ρ) | N_T | m | Derivation | CCF |
|------|---------|-----------|-----|---|------------|-----|
| A clonal diploid | 0.40 | 0.80 | 2 | 1 | 0.40·(1.6+0.4)/0.80 = 0.40·2.0/0.80 | 1.0 |
| B subclonal | 0.20 | 0.80 | 2 | 1 | 0.20·2.0/0.80 = 0.40/0.80 | 0.5 |
| C multi-copy m=2 | 0.50 | 1.0 | 4 | 2 | 0.50·(4.0+0.0)/(1.0·2) = 2.0/2.0 | 1.0 |
| D purity 0.5 | 0.25 | 0.50 | 2 | 1 | 0.25·(1.0+1.0)/0.50 = 0.50/0.50 | 1.0 |
| E noise cap | 0.471 | 1.0 | 2 | 1 | 0.471·2.0/1.0 = 0.942 (uncapped); >1 case capped | 0.942 |
| F overshoot→cap | 0.60 | 0.80 | 2 | 1 | 0.60·2.0/0.80 = 1.50 raw → capped | 1.0 (raw 1.5) |

### Dataset: Derived 1D CCF clustering case (deterministic Lloyd k=2)

**Source:** Derivation from Lloyd 1982 assignment/update steps; clonal cluster = highest centroid (PMC7867630).

| CCF values | k | Expected centroids | Expected assignments | Clonal cluster |
|------------|---|--------------------|--------------------|----------------|
| {1.0, 0.98, 0.96, 0.50, 0.48, 0.52} | 2 | {0.50, 0.98} | low: indices 3,4,5; high: 0,1,2 | high (centroid 0.98) |

---

## Assumptions

1. **ASSUMPTION: CCF reported value is capped to [0,1].** The raw formula can exceed 1 (CNAqc, 1.06). The registry invariant is 0 ≤ CCF ≤ 1; we report min(1, raw) as the bounded CCF (consistent with the McGranahan clonal definition that a mutation in all cancer cells has CCF = 1) while also exposing the uncapped raw value. Justification: invariant + McGranahan clonal definition; no source forbids capping.
2. **ASSUMPTION: 1D clustering algorithm = deterministic Lloyd k-means with quantile seeding.** Sources name CCF clustering broadly (Dirichlet process, variational beta mixtures) but the unit requires a *deterministic, well-defined* 1D method (per task constraints, no fabricated Dirichlet process). Lloyd's k-means (Lloyd 1982) is fully specified; determinism is achieved by sorting values and seeding centroids at evenly-spaced quantiles (no RNG). The clonal-cluster rule (highest centroid) is source-backed (PMC7867630).

---

## Recommendations for Test Coverage

1. **MUST Test:** EstimateCCF on cases A–D returns the exact derived CCF (1.0, 0.5, 1.0, 1.0) within 1e-10 — Evidence: McGranahan 2016 / PMC7867630 / Zheng 2022 formula.
2. **MUST Test:** EstimateCCF caps reported CCF at 1.0 when raw > 1 (case F) and exposes raw (1.5) — Evidence: invariant + CNAqc CCF=1.06.
3. **MUST Test:** EstimateCCF multiplicity must be ≥1 and ≤ tumor copy number; invalid inputs throw — Evidence: multiplicity def (PMC7867630).
4. **MUST Test:** EstimateCCF rejects purity ∉ (0,1], negative/>1 VAF, copy number < 1 — Evidence: formula domain.
5. **MUST Test:** ClusterCCFValues on the derived dataset returns the exact centroids/assignments and identifies the highest-centroid cluster as clonal — Evidence: Lloyd 1982 + PMC7867630.
6. **MUST Test:** ClusterCCFValues is deterministic (identical output across repeated runs / shuffled input grouped by value) — Evidence: deterministic seeding.
7. **SHOULD Test:** ClusterCCFValues with k=1 returns one cluster (mean) — Rationale: boundary.
8. **SHOULD Test:** EstimateCCF/ClusterCCFValues null/empty handling — Rationale: documented failure modes.

---

## References

1. Tarabichi M, Salcedo A, Deshwar AG, et al. (2021). A practical guide to cancer subclonal reconstruction from DNA sequencing. *Nature Methods* 18:144–155. https://pmc.ncbi.nlm.nih.gov/articles/PMC7867630/ (DOI: 10.1038/s41592-020-01013-2)
2. Zheng J, Wang K, et al. (2022). Estimation of cancer cell fractions and clone trees from multi-region sequencing of tumors. *Bioinformatics* 38(15):3677–3683. https://academic.oup.com/bioinformatics/article/38/15/3677/6596597 (DOI: 10.1093/bioinformatics/btac367)
3. McGranahan N, Furness AJS, Rosenthal R, et al. (2016). Clonal neoantigens elicit T cell immunoreactivity and sensitivity to immune checkpoint blockade. *Science* 351(6280):1463–1469. https://www.science.org/doi/10.1126/science.aaf1490 (DOI: 10.1126/science.aaf1490)
4. Caravagna G, et al. CNAqc — Computation of Cancer Cell Fractions. https://caravagnalab.github.io/CNAqc/articles/a4_ccf_computation.html (accessed 2026-06-15)
5. Lloyd SP (1982). Least squares quantization in PCM. *IEEE Trans. Inf. Theory* 28(2):129–137. https://doi.org/10.1109/TIT.1982.1056489 (algorithm steps retrieved via https://en.wikipedia.org/wiki/K-means_clustering)

---

## Change History

- **2026-06-15**: Initial documentation.
