# Evidence Artifact: ONCO-PURITY-001

**Test Unit ID:** ONCO-PURITY-001
**Algorithm:** Tumor Purity Estimation (from somatic SNV VAF / allele-specific copy number)
**Date Collected:** 2026-06-14

---

## Online Sources

### CNAqc — Quality Control of allele-specific copy numbers, mutations and tumour purity (vignette)

**URL:** https://caravagnalab.github.io/CNAqc/articles/CNAqc.html
**Retrieved by:** WebFetch of the URL above (CNAqc package vignette), 2026-06-14, prompting for the expected-VAF formula and its variable definitions.
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference implementation / official package documentation), backing the peer-reviewed Genome Biology 2024 paper (rank 1).

**Key Extracted Points:**

1. **Expected-VAF formula (verbatim):** `v_m(c) = mπc / [2(1-π) + π(n_A+n_B)]`.
2. **Variable definitions (verbatim phrasing):** `m` = "mutations present in m copies of the tumour genome" (multiplicity); `π` = "tumour purity"; `c` = "mutations present in a percentage 0<c<1 of tumour cells" (cancer cell fraction / clonality); a segment is written `n_A:n_B` for the allele-specific copy numbers `n_A` and `n_B`.
3. **Denominator interpretation (verbatim phrasing):** `2(1−π)` represents "a healthy diploid normal" contribution and `π(n_A+n_B)` is "the proportion of all reads from the tumour"; total denominator = mean number of allele copies per cell across the mixture.

### CNAqc — Computational validation of clonal and subclonal CNAs (Genome Biology 2024, paper)

**URL:** https://link.springer.com/article/10.1186/s13059-024-03170-5 (search-indexed full text; abstract/figure captions returned via web search snippet of this DOI)
**Retrieved by:** WebSearch queries `CNAqc tumor purity expected VAF formula …` and `CNAqc "expected VAF" diploid heterozygous "purity" peak 0.5 …`, 2026-06-14, surfacing the paper's worked numeric examples. (Full-text and bioRxiv PDF returned HTTP 403 to direct fetch; numeric examples taken from the indexed snippets of this DOI.)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, Genome Biology).

**Key Extracted Points:**

1. **General clonal/subclonal form (verbatim from snippet):** for mutations on (possibly two-state) CNAs, `v = (m₁ρ₁ + m₂ρ₂)π / { 2(1−π) + π[(n_{A,1}+n_{B,1})ρ₁ + (n_{A,2}+n_{B,2})ρ₂] }`. For a single clonal state this reduces to the vignette form with c = 1.
2. **Worked example — diploid heterozygous, purity 60%:** "real purity of 60%" for a heterozygous diploid mutation (m=1, n_A+n_B=2) corresponds to expected VAF 30%; the purity tolerance band 55–65% maps to "a VAF range of 27.50–32.5%". (Confirms v = π/2 ⇒ π = 2·VAF.)
3. **Worked example — 2:1 segment, purity = 1:** "For 2:1 segments, there are two peaks of clonal mutations at 33% and 66% VAF" — i.e. m=1 ⇒ 1/3 ≈ 0.333 and m=2 ⇒ 2/3 ≈ 0.667 with n_A+n_B = 3, π = 1.

### FACETS — allele-specific copy number and clonal heterogeneity (NAR 2016)

**URL:** https://academic.oup.com/nar/article/44/16/e131/2460163
**Retrieved by:** WebFetch of the URL above, 2026-06-14, prompting for the purity/copy-number mixing model and full citation.
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, Nucleic Acids Research).

**Key Extracted Points:**

1. **Mixing model (verbatim phrasing):** tumor parental copy numbers adjusted for cellular fraction Φ are `m* = mΦ + (1 − Φ)` and `p* = pΦ + (1 − Φ)`; the model mixes a "normal diploid (1,1)" genotype with an "aberrant (m,p)" genotype at mixing proportion Φ (cellular fraction). This independently confirms the `2(1−π) + π·n_tot` denominator structure used by CNAqc (normal contributes 2 copies weighted 1−π).

### ABSOLUTE — Absolute quantification of somatic DNA alterations (Carter et al. 2012)

**URL / retrieval:** Citation record retrieved via Europe PMC REST API
`https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:22544022&format=json&resultType=core` (WebFetch), 2026-06-14. (PubMed/Nature direct pages returned reCAPTCHA / 403; the EuropePMC core record supplied the verified bibliographic metadata.)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, Nature Biotechnology).

**Key Extracted Points:**

1. **Verified citation:** Carter SL, Cibulskis K, Helman E, McKenna A, Shen H, Zack T, Laird PW, Onofrio RC, Winckler W, Weir BA, Beroukhim R, Pellman D, Levine DA, Lander ES, Meyerson M, Getz G. "Absolute quantification of somatic DNA alterations in human cancer." *Nature Biotechnology* 30(5):413–421, 2012. DOI: 10.1038/nbt.2203.
2. **Context (from search overview):** ABSOLUTE converts allelic fractions of point mutations into per-cancer-cell allele counts (cellular multiplicity) "by correcting for sample purity and local copy-numbers" — the same purity/copy-number correction inverted here to estimate purity.

---

## Documented Corner Cases and Failure Modes

### From CNAqc

1. **Multiplicity ambiguity on amplified segments:** on a 2:1 (n_tot=3) segment clonal mutations form two VAF peaks (1/3 and 2/3) because multiplicity m may be 1 or 2; purity cannot be inferred from VAF alone without knowing m and the copy-number state. Copy-neutral diploid heterozygous (1:1) loci avoid this — m = 1, n_tot = 2 — which is why purity = 2·VAF is the robust closed form there.
2. **Subclonal mutations (c < 1):** for c < 1 the VAF is depressed (v ∝ c); treating a subclonal VAF as clonal underestimates purity. Purity must be estimated from clonal mutations.

### From the checklist by-area definition (ONCO-PURITY-001)

1. **Purity < 0.1 (below detection limit):** very low purity yields VAFs near sequencing noise.
2. **No heterozygous SNPs / no informative variants:** with no usable variant evidence purity is undefined.
3. **High stromal contamination:** equivalent to low purity.

---

## Test Datasets

### Dataset: CNAqc clonal heterozygous diploid worked example

**Source:** CNAqc (Genome Biology 2024), doi:10.1186/s13059-024-03170-5; CNAqc vignette.

| Parameter | Value |
|-----------|-------|
| Segment (n_A:n_B) | 1:1 (copy-neutral diploid, n_tot = 2) |
| Multiplicity m | 1 (heterozygous somatic SNV) |
| Clonality c | 1 (clonal) |
| Expected VAF at purity 0.60 | 0.30 (= π/2) |
| Purity from VAF 0.30 | 0.60 (= 2·VAF) |
| Purity tolerance band 0.55–0.65 | VAF band 0.275–0.325 |

### Dataset: CNAqc 2:1 amplified segment worked example (purity = 1)

**Source:** CNAqc (Genome Biology 2024), doi:10.1186/s13059-024-03170-5.

| Parameter | Value |
|-----------|-------|
| Segment (n_A:n_B) | 2:1 (n_tot = 3) |
| Purity π | 1.0 |
| Clonal VAF, m = 1 | 1/3 ≈ 0.3333… |
| Clonal VAF, m = 2 | 2/3 ≈ 0.6667… |

---

## Assumptions

1. **ASSUMPTION: VAF-only purity estimator uses the copy-neutral diploid heterozygous model.** `EstimatePurityFromVAF` assumes the supplied variants are clonal (c = 1) heterozygous (m = 1) somatic SNVs at copy-neutral diploid (n_tot = 2) loci, giving the closed form ρ = 2·VAF. This is the textbook special case and the band the CNAqc example uses; it is stated explicitly in the API contract. It is a modeling scope choice, not an invented numeric constant — the formula itself is fully source-derived.
2. **ASSUMPTION: aggregation across variants uses the median VAF.** When multiple clonal heterozygous SNVs are supplied, their per-variant purity estimates are combined by the median (robust to subclonal/outlier VAFs). The literature establishes the per-variant relation; the choice of a robust central estimator over the set is a documented, non-correctness-affecting aggregation policy (it does not change the single-variant formula). Recorded as an assumption for transparency.

---

## Recommendations for Test Coverage

1. **MUST Test:** `EstimatePurityFromVAF` on a clonal heterozygous diploid SNV with VAF 0.30 returns purity 0.60. — Evidence: CNAqc worked example (π/2 = 0.30 ⇔ ρ = 2·VAF).
2. **MUST Test:** `EstimatePurityFromVAF` boundary VAF 0.50 → purity 1.0; VAF 0.0 → purity 0.0. — Evidence: ρ = 2·VAF closed form, INV 0 ≤ ρ ≤ 1.
3. **MUST Test:** `EstimatePurity` (allele-specific) inverting v = mπ/[2(1−π)+π·n_tot] recovers π for a clonal SNV on a 2:1 (n_tot=3) segment: π=1, m=1, v=1/3 ⇒ π=1; m=2, v=2/3 ⇒ π=1. — Evidence: CNAqc 2:1 peaks 33%/66%.
4. **MUST Test:** `EstimatePurity` on diploid heterozygous segment (n_tot=2, m=1) with VAF 0.30 returns 0.60 (must agree with VAF-only estimator). — Evidence: CNAqc 60%/30% example.
5. **MUST Test:** invalid inputs (VAF outside [0,1], VAF > 0.5 for the diploid model implying purity > 1, empty variant list, non-positive copy number) throw / are rejected. — Evidence: invariant 0 ≤ ρ ≤ 1; formula domain.
6. **SHOULD Test:** median aggregation across several heterozygous SNVs with mixed VAFs returns the median-derived purity. — Rationale: robustness corner case.
7. **COULD Test:** purity below detection (VAF near 0) yields a small purity near 0 without error. — Rationale: low-purity edge case from checklist.

---

## References

1. Antonello A, Bergamin R, Calonaci N, Househam J, Milite S, Williams MJ, Anselmi F, d'Onofrio A, Sundaram V, Sosinsky A, Cross WCH, Caravagna G (2024). Computational validation of clonal and subclonal copy number alterations from bulk tumor sequencing using CNAqc. *Genome Biology* 25(1):38. https://doi.org/10.1186/s13059-024-03170-5
2. CNAqc package vignette (Caravagna lab). Quality control of allele-specific copy numbers, mutations and tumour purity. https://caravagnalab.github.io/CNAqc/articles/CNAqc.html (accessed 2026-06-14)
3. Carter SL, Cibulskis K, Helman E, et al. (2012). Absolute quantification of somatic DNA alterations in human cancer. *Nature Biotechnology* 30(5):413–421. https://doi.org/10.1038/nbt.2203
4. Shen R, Seshan VE (2016). FACETS: allele-specific copy number and clonal heterogeneity analysis tool for high-throughput DNA sequencing. *Nucleic Acids Research* 44(16):e131. https://doi.org/10.1093/nar/gkw520

---

## Change History

- **2026-06-14**: Initial documentation.
