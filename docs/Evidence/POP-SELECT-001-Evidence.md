# Evidence Artifact: POP-SELECT-001

**Test Unit ID:** POP-SELECT-001
**Algorithm:** Selection Signature Detection (integrated Haplotype Score, iHS)
**Date Collected:** 2026-06-13

---

## Online Sources

### Voight, Kudaravalli, Wen & Pritchard (2006) — A Map of Recent Positive Selection in the Human Genome

**URL:** https://web.stanford.edu/group/pritchardlab/publications/VoightEtAl06.pdf
**Accessed:** 2026-06-13 (PDF fetched and extracted with `pdftotext`; primary source for iHS)
**Authority rank:** 1 (peer-reviewed paper, PLoS Biology 4(3):e72)

**Key Extracted Points:**

1. **iHH definition (p.0447, lines 102–110 of extracted text):** "the area under the EHH curve … This integrated EHH (iHH) (summed over both directions away from the core SNP) will be denoted iHH_A or iHH_D, depending on whether it is computed with respect to the ancestral or derived core allele."
2. **Unstandardized iHS formula (lines 111–118):** "unstandardized iHS = ln(iHH_A / iHH_D)". Numerator = ancestral, denominator = derived.
3. **Sign interpretation (lines 182–185):** "When the rate of EHH decay is similar on the ancestral and derived alleles, iHH_A/iHH_D ≈ 1, and hence the unstandardized iHS is ≈ 0. Large negative values indicate unusually long haplotypes carrying the derived allele; large positive values indicate long haplotypes carrying the ancestral allele."
4. **Standardization (lines 188–214):** "iHS = (ln(iHH_A/iHH_D) − E_p[ln(iHH_A/iHH_D)]) / SD_p[ln(iHH_A/iHH_D)]. The expectation and standard deviation … are estimated from the empirical distribution at SNPs whose derived allele frequency p matches the frequency at the core SNP. The iHS is constructed to have an approximately standard normal distribution and hence the sizes of iHS signals from different SNPs are directly comparable regardless of the allele frequencies at those SNPs."
5. **Integration method & cutoff (Materials and Methods, "Calculation of iHS"):** "The EHH values at successive SNPs are joined by straight lines, and then we compute the total area under each curve, between the nearest points to the left and right of the core SNP where the EHH drops below 0.05." (straight lines = trapezoidal rule; cutoff = 0.05).
6. **Genome-wide scan criterion (Materials and Methods, "Identifying candidate signals of selection"):** "the size of signal in each region was quantified by the proportion of SNPs with |iHS| > 2." Standard window size = 50 SNPs.
7. **MAF filter (Materials and Methods):** "iHS was computed … for every SNP with ancestral state information and with minor allele frequency > 5%."

### Sabeti et al. (2002) — Detecting recent positive selection from haplotype structure

**URL:** https://pubmed.ncbi.nlm.nih.gov/12397357/ (PubMed record, accessed 2026-06-13);
formula confirmed via Szpiech & Hernandez (2014) and the rehh vignette (see below). The reich.hms.harvard.edu PDF mirror returned HTTP 403 and could not be opened this session.
**Authority rank:** 1 (peer-reviewed paper, Nature 419:832–837; originator of EHH)

**Key Extracted Points:**

1. **EHH definition:** EHH is the probability that two randomly chosen chromosomes carrying the core allele are homozygous (identical) over a surrounding chromosomal region; computed in a stepwise manner for each haplotype length ("extended HH").

### Szpiech & Hernandez (2014) — selscan: an efficient multi-threaded program for EHH-based scans

**URL:** https://arxiv.org/pdf/1403.6854
**Accessed:** 2026-06-13 (PDF fetched and extracted with `pdftotext`; reference implementation, MBE 31(10):2824–2827)
**Authority rank:** 3 (established reference implementation; restates the formal definitions)

**Key Extracted Points:**

1. **EHH formula (Eq. 1):** EHH(x_i) = Σ_{h∈C(x_i)} C(n_h,2) / C(n,2), where n_h is the number of observed haplotypes of type h and n the sample size.
2. **Core-conditioned EHH (Eq. 3):** EHH_c(x_i) = Σ_{h∈H_c(x_i)} C(n_h,2) / C(n_c,2), where n_c is the number of chromosomes carrying core allele c. Example: C(x_1) = {11,10,00,01}.
3. **iHH by trapezoidal quadrature (Eq. 4/7):** iHH_c = Σ_downstream ½(EHH_c(x_{i-1}) + EHH_c(x_i))·g(x_{i-1},x_i) + Σ_upstream ½(…)·g(…); g = genetic distance between markers.
4. **Truncation cutoff:** "the sums … are truncated at x_i — the marker at which the EHH … is EHH(x_i) < 0.05."
5. **Sign convention difference (explicit note):** selscan uses unstandardized iHS = ln(iHH_1/iHH_0) (derived/ancestral) and states: "this definition differs slightly from that in Voight et al. (2006), where unstandardized iHS is defined with iHH_1 and iHH_0 swapped." (Confirms Voight = ln(iHH_A/iHH_D).)
6. **Standardization (Eq. 6):** iHS = (ln(iHH_1/iHH_0) − E_p[…]) / SD_p[…], normalized in frequency bins across the genome.

### rehh package vignette (Gautier, Klassmann & Vitalis)

**URL:** https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html
**Accessed:** 2026-06-13 (HTML fetched via WebFetch; established R reference implementation)
**Authority rank:** 3 (reference implementation, CRAN package rehh)

**Key Extracted Points:**

1. **EHH formula (algebraically identical to selscan):** EHH = (1 / (n_a(n_a−1))) · Σ_k n_k(n_k−1), where n_a is the number of core-allele chromosomes and n_k the count of the k-th shared haplotype. Since C(n,2)=n(n−1)/2 the factor of 2 cancels, matching selscan Eq. 3.
2. **iHH:** "The integrated EHH (iHH) is defined as the area under the EHH curve … using the trapezoidal rule."
3. **Default cutoff:** parameter `limehh` default value is 0.05; "EHH is usually computed only for a region [where] it surpasses a given threshold (e.g., EHH > 0.05)."
4. **Worked numeric example (`calc_ehh()` output for SNP F1205400):** IHH_A = 284,429.9; IHH_D = 2,057,107.4 ⇒ unstandardized iHS (Voight) = ln(284429.9/2057107.4) = −1.97857 (verified by computation this session).

---

## Documented Corner Cases and Failure Modes

### From Voight et al. (2006)

1. **Monomorphic / no ancestral state:** iHS is computed only for SNPs that are polymorphic and have ancestral state information and MAF > 5%; otherwise no value is reported.
2. **Long gaps:** "If the region spanned by EHH > 0.05 reached a chromosome end or the start of a gap > 200 kb, then no iHS value was reported." (Gap handling is a data-curation detail; not part of the core score for in-memory inputs.)
3. **Balanced decay:** when EHH decays at similar rates, iHH_A/iHH_D ≈ 1 ⇒ unstandardized iHS ≈ 0.

### From Szpiech & Hernandez (2014)

1. **EHH integration truncation:** integration stops at the marker where EHH first drops below 0.05; the area is summed over both directions.
2. **Sign convention pitfall:** the iHS sign depends on which allele is the numerator; Voight (ancestral/derived) and selscan (derived/ancestral) differ by sign. Implementations must declare their convention.

---

## Test Datasets

### Dataset: rehh worked iHH ratio (SNP F1205400)

**Source:** rehh vignette, `calc_ehh()` output (Gautier et al.)

| Parameter | Value |
|-----------|-------|
| IHH_A (ancestral) | 284429.9 |
| IHH_D (derived) | 2057107.4 |
| Unstandardized iHS (Voight, ln(A/D)) | −1.978569274 (computed) |

### Dataset: Constructed deterministic haplotype panel (derived under-selection vs neutral ancestral)

**Source:** derived by hand from the formulas above (EHH Eq. 3, trapezoidal iHH, ln(iHH_A/iHH_D)).

| Parameter | Value |
|-----------|-------|
| Haplotypes (core index 2, positions 0,10,20,30,40) | 3× `AA1GG` (derived, identical flanks), `TC0TC`,`GA0AG`,`CT0CA` (ancestral, all distinct) |
| EHH_derived at each flank marker | 1.0 (all identical) |
| EHH_ancestral at first flank marker | 0.0 (all distinct) → truncated |
| iHH_D | (1+1)/2·10 + (1+1)/2·10 (right) + same (left) = 40.0 |
| iHH_A | (1+0)/2·10 (right, then EHH<0.05 stop) + (1+0)/2·10 (left) = 10.0 |
| Unstandardized iHS = ln(10/40) | −1.386294361 |
| Derived allele frequency | 0.5 |

### Dataset: EHH unit values (selscan Eq. 3)

**Source:** Szpiech & Hernandez (2014) Eq. 3, worked by hand.

| Input (extended haplotypes) | EHH |
|-----------------------------|-----|
| `11`,`11`,`11`,`10` (n_c=4) | (C(3,2)+C(1,2))/C(4,2) = 3/6 = 0.5 |
| `00`,`00`,`01`,`01` (n_c=4) | (1+1)/6 = 0.333333 |
| single haplotype | 1.0 (trivially homozygous) |
| three all-distinct haplotypes | 0.0 |

---

## Assumptions

1. **ASSUMPTION: Standard-deviation estimator in standardization** — Voight et al. specify "SD of ln(iHH_A/iHH_D)" within a frequency bin but do not state whether it is the population (N) or sample (N−1) standard deviation. The implementation uses the sample standard deviation (N−1), the conventional unbiased estimator for an empirical distribution. This affects only the magnitude scaling of standardized scores, not the sign or relative ordering, and not the unstandardized iHS (the canonical evidence-backed value).

2. **ASSUMPTION: Frequency-bin width** — Voight bins by derived allele frequency "matching" the core SNP but the exact bin width used for the published genome scan is not stated in the retrieved text. The implementation defaults to 20 equal-width bins (width 0.05), matching the rehh `freqbin = 0.05` convention. Bin count is an explicit parameter, so callers can override it.

---

## Recommendations for Test Coverage

1. **MUST Test:** `CalculateEhh` equals Σ C(n_h,2)/C(n_c,2) on the worked unit values (0.5, 0.3333, 1.0, 0.0). — Evidence: selscan Eq. 3 / rehh.
2. **MUST Test:** `CalculateIHS` on the constructed panel returns iHH_A=10, iHH_D=40, unstandardized iHS=ln(0.25)=−1.386294361, derived freq 0.5. — Evidence: Voight iHH (trapezoid + 0.05 cutoff) + ln(iHH_A/iHH_D).
3. **MUST Test:** balanced decay ⇒ unstandardized iHS ≈ 0. — Evidence: Voight lines 182–185.
4. **MUST Test:** sign convention — long derived haplotype ⇒ negative iHS (Voight). — Evidence: Voight lines 182–185.
5. **MUST Test:** `StandardizeIHS` reproduces (x−mean)/sd within a frequency bin; single-element bin ⇒ 0. — Evidence: Voight Eq. (standardization).
6. **MUST Test:** `ScanForSelection` window proportion = (count of |iHS|>2) / window size. — Evidence: Voight "proportion of SNPs with |iHS| > 2".
7. **MUST Test:** edge cases — null inputs, empty haplotypes, monomorphic core (throws), inconsistent lengths (throws), invalid core allele (throws), coreIndex out of range (throws).
8. **SHOULD Test:** EHH of a single chromosome = 1; empty sample = 0. — Rationale: boundary of the combinatorial formula.
9. **COULD Test (property):** ln(iHH_A/iHH_D) = −ln(iHH_D/iHH_A) (Voight vs selscan sign symmetry). — Rationale: O(n·h) algorithm invariant.

---

## References

1. Voight BF, Kudaravalli S, Wen X, Pritchard JK (2006). A Map of Recent Positive Selection in the Human Genome. PLoS Biology 4(3):e72. https://doi.org/10.1371/journal.pbio.0040072 (PDF: https://web.stanford.edu/group/pritchardlab/publications/VoightEtAl06.pdf)
2. Sabeti PC, Reich DE, Higgins JM, et al. (2002). Detecting recent positive selection in the human genome from haplotype structure. Nature 419:832–837. https://pubmed.ncbi.nlm.nih.gov/12397357/
3. Szpiech ZA, Hernandez RD (2014). selscan: an efficient multi-threaded program to perform EHH-based scans for positive selection. Molecular Biology and Evolution 31(10):2824–2827. https://arxiv.org/pdf/1403.6854 (https://doi.org/10.1093/molbev/msu211)
4. Gautier M, Klassmann A, Vitalis R. rehh package vignette. CRAN. https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html (accessed 2026-06-13)

---

## Change History

- **2026-06-13**: Initial documentation.
