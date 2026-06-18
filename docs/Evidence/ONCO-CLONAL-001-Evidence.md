# Evidence Artifact: ONCO-CLONAL-001

**Test Unit ID:** ONCO-CLONAL-001
**Algorithm:** Clonal vs Subclonal Mutation Classification (cancer cell fraction posterior)
**Date Collected:** 2026-06-14

---

## Online Sources

### Landau et al. (2013) — Evolution and Impact of Subclonal Mutations in Chronic Lymphocytic Leukemia (Cell)

**URL:** https://www.cell.com/fulltext/S0092-8674(13)00071-8 (open-access manuscript PDF retrieved from https://wulab.dfci.harvard.edu/sites/default/files/23415222.pdf)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed primary paper)
**How retrieved:** WebSearch query `Landau 2013 "Evolution and Impact of Subclonal Mutations" PMC full text cancer cell fraction clonal classification ABSOLUTE`; then fetched the Wu Lab PDF and extracted text locally with `pdftotext`.

**Key Extracted Points:**

1. **Clonal/subclonal classification rule (verbatim):** "We classified a mutation as clonal if the CCF harboring it was >0.95 with probability > 0.5, and subclonal otherwise" (Results, near Figure 2A). Restated in the Extended Experimental Procedures: "Mutations were thereafter classified as clonal based on the posterior probability that the CCF exceeded 0.95, and subclonal otherwise."
2. **Expected allele fraction formula (verbatim):** "The expected allele-fraction f of a mutation present in one copy in a fraction c of cancer cells is calculated by … f(c) = αc/(2(1 − α) + αq)", where α = sample purity, c = CCF, q = absolute somatic copy-number at the locus.
3. **Posterior construction (verbatim):** "Then P(c) α Binom(a|N,f(c)), assuming a uniform prior on c", with "c ∈ [0.01,1]"; the CCF posterior "was then obtained by calculating these values over a regular grid of 100 c values and normalizing by dividing them by their sum".
4. **Inputs:** "a somatic mutation observed in a of N sequencing reads on a locus of absolute somatic copy-number q in a sample of purity α."
5. **Subclonal CCF spread (sanity reference):** "For sSNVs designated as subclonal, median CCF was 0.49", confirming subclonal mutations sit well below 1.

### Satas, Zaccaria, El-Kebir et al. (2021) — DeCiFering the Elusive Cancer Cell Fraction (Cell Systems)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC8542635/
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed primary paper)
**How retrieved:** WebSearch query `cancer cell fraction clonal subclonal classification confidence interval 0.95 McGranahan 2015 definition`; then WebFetch of the PMC article.

**Key Extracted Points:**

1. **Multiplicity-general CCF formula (Eq. 1, verbatim):** "c ≈ (1/ρ)·(ρ·N_tot + (1−ρ)·2)/M · v̂", where ρ = tumor purity, N_tot = average total copy number in cancer cells, M = number of SNV copies in mutated cells (multiplicity), v̂ = VAF estimate. This generalises Landau's M = 1 expected-allele-fraction relation to arbitrary multiplicity M (inverting f(c) = ρ·M·c / (2(1−ρ) + ρ·N_tot)).
2. **Clonal/subclonal definition (verbatim):** "SNVs are classified as clonal if they are inferred to be present in all cells in a tumor sample (CCF ≈ 1) or subclonal if they are inferred to be present only in a subpopulation (CCF ≪ 1)."

---

## Documented Corner Cases and Failure Modes

### From Landau et al. (2013)

1. **CCF lower bound > 0:** the grid is c ∈ [0.01, 1] — a detected mutation is present in at least one cancer cell, so CCF cannot be exactly 0.
2. **Probabilistic, not point, classification:** classification is by the posterior probability mass above 0.95, not by the point estimate alone, so a point estimate near 1 with wide uncertainty (shallow coverage) may still be subclonal.

### From Satas et al. (2021)

1. **Multiplicity matters:** for M > 1 (e.g. mutation on both copies after copy-neutral LOH or after duplication), the same VAF corresponds to a lower CCF; ignoring M overestimates CCF.

---

## Test Datasets

### Dataset: Posterior-grid classification cases (derived from Landau 2013 model)

**Source:** Landau et al. (2013), f(c) = αc/(2(1−α)+αq); P(c) ∝ Binom(a|N,f(c)); uniform prior; 100-point grid c∈[0.01,1]; clonal iff P(CCF>0.95) > 0.5. Multiplicity generalisation f(c) = αMc/(2(1−α)+αq) per Satas et al. (2021) Eq. 1. Expected values computed independently of the implementation by direct grid evaluation.

| Case | a | N | q | M | ρ | CCF mean | P(CCF>0.95) | Status |
|------|---|---|---|---|----|----------|-------------|--------|
| A1 | 300 | 300 | 2 | 1 | 1.0 | 0.999486 | 1.000000 | Clonal |
| B2 | 400 | 1000 | 2 | 1 | 0.8 | 0.972455 | 0.864167 | Clonal |
| C1 | 240 | 1000 | 2 | 1 | 0.8 | 0.601297 | 0.000000 | Subclonal |
| D  | 200 | 1000 | 2 | 1 | 1.0 | 0.401198 | ≈0 | Subclonal |
| E  | 100 | 100 | 2 | 2 | 1.0 | 0.994330 | 0.998016 | Clonal |

### Dataset: Point-estimate clonal threshold (IdentifyClonalMutations)

**Source:** Landau et al. (2013) — clonal iff CCF > 0.95.

| CCF input | 0.96 | 0.95 | 1.00 | 0.50 | 0.951 |
|-----------|------|------|------|------|-------|
| Clonal (CCF > 0.95)? | yes | no (boundary, not strict) | yes | no | yes |

Clonal indices: {0, 2, 4}.

---

## Assumptions

1. **ASSUMPTION: Registry `ploidy` parameter is the per-variant local copy number `q`.** The registry stub signature is `ClassifyClonality(variants, purity, ploidy)`. Landau's model uses the **per-locus** absolute somatic copy number q, not a single genome-wide ploidy scalar. The canonical method therefore carries q per variant (`ClonalityVariant.LocalCopyNumber`) and takes `(variants, purity)`. This mirrors the prior ONCO-WGD decision (registry note at line 4405) where a registry scalar was superseded by per-segment data to match the authoritative definition. Non-correctness-affecting (API shape only): the numerical rule and outputs are exactly Landau's.

---

## Recommendations for Test Coverage

1. **MUST Test:** Clonal call for a deep-coverage variant whose posterior mass above 0.95 exceeds 0.5 (cases A1, B2, E). — Evidence: Landau (2013), clonal iff P(CCF>0.95) > 0.5.
2. **MUST Test:** Subclonal call for a variant with CCF well below 1 (cases C1, D). — Evidence: Landau (2013).
3. **MUST Test:** Multiplicity M=2 raises CCF for the same VAF (case E). — Evidence: Satas (2021) Eq. 1.
4. **MUST Test:** Invariant ClonalCount + SubclonalCount = total; ClonalFraction = ClonalCount/total. — Evidence: registry invariant.
5. **MUST Test:** `IdentifyClonalMutations` returns indices with CCF strictly > 0.95 (boundary 0.95 excluded). — Evidence: Landau (2013), CCF > 0.95.
6. **SHOULD Test:** Null inputs, purity ∉ (0,1], invalid read counts / copy number / multiplicity, CCF ∉ [0,1] throw. — Rationale: documented domain and standard validation.
7. **COULD Test:** Empty variant set → empty calls, counts 0, ClonalFraction 0. — Rationale: degenerate-input contract.

---

## References

1. Landau DA, Carter SL, Stojanov P, et al. (2013). Evolution and Impact of Subclonal Mutations in Chronic Lymphocytic Leukemia. *Cell* 152(4):714–726. https://doi.org/10.1016/j.cell.2013.01.019 (PMC: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3575604/)
2. Satas G, Zaccaria S, El-Kebir M, Raphael BJ (2021). DeCiFering the Elusive Cancer Cell Fraction in Tumor Heterogeneity and Evolution. *Cell Systems* 12(10):1004–1018. https://doi.org/10.1016/j.cels.2021.07.006 (PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC8542635/)

---

## Change History

- **2026-06-14**: Initial documentation.
