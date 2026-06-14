# Evidence Artifact: ONCO-VAF-001

**Test Unit ID:** ONCO-VAF-001
**Algorithm:** Variant Allele Frequency Analysis (empirical VAF, Wilson score confidence interval, purity/ploidy correction)
**Date Collected:** 2026-06-14

---

## Online Sources

### GATK Mutect2 — AlleleFraction / FAQ (allele fraction definition)

**URL:** https://gatk.broadinstitute.org/hc/en-us/articles/360050722212-FAQ-for-Mutect2
**Accessed:** 2026-06-14 (retrieved via WebSearch query "GATK Mutect2 AF allele fraction definition alt reads divided by total reads VCF FORMAT")
**Authority rank:** 3 (reference implementation documentation)

**Key Extracted Points:**

1. **AF FORMAT field:** Mutect2 emits a per-sample `AF` FORMAT field documented as "Allele fractions of alternate alleles in the tumor".
2. **Mutect2 AF is a model estimate, not the raw ratio:** Mutect2's `AF` is *not* the empirical ratio alt/total; it is a posterior/mean-likelihood estimate marginalised over allele fractions. The *empirical* allele fraction (alt-supporting reads / total reads at the locus, derived from the `AD` field) is the simple, model-free quantity. This unit implements the empirical VAF = altReads / totalReads, matching the `AD`-based definition (alt AD / sum AD), which is the standard textbook VAF.

### Wikipedia — Binomial proportion confidence interval (Wilson score interval)

**URL:** https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
**Accessed:** 2026-06-14 (retrieved via WebFetch of the URL above)
**Authority rank:** 4 (Wikipedia citing the primary source Wilson 1927; the cited primary is used)

**Key Extracted Points:**

1. **Wilson score interval (verbatim form extracted from the page):**
   p ≈ (1 / (1 + z²/n)) · ( p̂ + z²/(2n) ± (z/(2n)) · √( 4n·p̂(1−p̂) + z² ) ).
   Equivalently, center = (p̂ + z²/(2n)) / (1 + z²/n) and margin = (z/(1+z²/n)) · √( p̂(1−p̂)/n + z²/(4n²) ).
2. **Attribution:** "The Wilson score interval was developed by E.B. Wilson (1927)." Primary reference: Wilson, E.B. (1927). *Probable inference, the law of succession, and statistical inference.* JASA 22(158):209–212.
3. **95% confidence z value:** the page states z₀.₀₅ = 1.96 for a 95% interval.
4. **Wald (normal-approximation) interval, for contrast:** p ≈ p̂ ± (z/√n)·√(p̂(1−p̂)); criticised for overshoot and zero-width intervals at extreme proportions — the Wilson interval avoids both.

### Tarabichi et al. — Principles of Reconstructing the Subclonal Architecture of Cancers (PMC5538405)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC5538405/
**Accessed:** 2026-06-14 (retrieved via WebFetch of the URL above)
**Authority rank:** 1 (peer-reviewed review, Cold Spring Harb Perspect Med)

**Key Extracted Points:**

1. **VAF / CCF / purity relationship (eq. extracted):** u_i = f_i / ( (1/ρ)·[ρ·n_tot,t,i + (1−ρ)·n_tot,n,i] )⁻¹ — i.e. f_i (VAF) is proportional to ρ·CCF·m divided by the average total copies per cell ρ·n_tot,t + (1−ρ)·n_tot,n, where ρ = purity, m = multiplicity, n_tot,t = tumor total copy number, n_tot,n = normal total copy number (= 2).
2. **Diploid heterozygous clonal case (extracted verbatim meaning):** "For a clonal mutation in a diploid region with no copy number alterations and heterozygous presentation, the expected VAF = purity × 0.5." Worked example: 80% purity ⇒ expected VAF = 0.4, "carried by half the reads that represent tumor DNA".

### CNAqc — Computational validation of CNAs from bulk tumor sequencing (Genome Biology 2024)

**URL:** https://link.springer.com/article/10.1186/s13059-024-03170-5 (formula text retrieved via WebSearch query "CNAqc expected VAF formula multiplicity purity copy number genome biology 2024")
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed, Genome Biology)

**Key Extracted Points:**

1. **Expected (clonal) VAF formula (extracted):** v = (m·π) / ( 2(1−π) + π·n_tot ), where m = mutation multiplicity, π = tumor purity, n_tot = nA+nB = tumor total copy number, normal ploidy = 2.
2. **Worked peak example:** for a 2:1 (n_tot = 3) segment there are clonal-mutation VAF peaks at 33% and 66%. For diploid (n_tot = 2, m = 1): v = π/2.
3. **Inversion used for purity correction:** solving for the per-tumour-copy mutant fraction m·CCF gives m·CCF = VAF·(2(1−π) + π·n_tot)/π. This is the purity/ploidy adjustment computed by `AdjustVAFForPurity`.

---

## Documented Corner Cases and Failure Modes

### From GATK Mutect2 docs

1. **VAF > 1 due to alignment artifacts:** raw alt counts may exceed total at a locus because of misaligned / overlapping reads; an empirical VAF computed naively could exceed 1. The empirical definition requires altReads ≤ totalReads; a count violating this is invalid input.
2. **Zero coverage (totalReads = 0):** VAF is undefined (0/0). No reads ⇒ no observation; division by zero must be guarded.

### From Wikipedia (Binomial CI)

1. **Extreme proportions (p̂ = 0 or 1):** the Wilson interval stays within [0,1] (no overshoot) and has non-zero width, unlike the Wald interval. Lower bound at p̂=0 is ≥ 0; upper bound at p̂=1 is ≤ 1.

### From CNAqc / Tarabichi

1. **Purity = 0:** the purity correction divides by π; π = 0 (pure normal) makes the adjusted fraction undefined.

---

## Test Datasets

### Dataset: Empirical VAF point values

**Source:** definition VAF = altReads/totalReads (GATK AD-based; samtools mpileup read counts).

| altReads | totalReads | VAF |
|----------|-----------|-----|
| 25 | 100 | 0.25 |
| 50 | 100 | 0.50 |
| 5 | 20 | 0.25 |
| 0 | 10 | 0.00 |
| 10 | 10 | 1.00 |

### Dataset: Wilson 95% score interval (z = 1.96)

**Source:** Wikipedia Binomial proportion confidence interval (Wilson 1927); values independently computed from the extracted formula.

| altReads | totalReads | center | lower | upper |
|----------|-----------|--------|-------|-------|
| 25 | 100 | 0.2592487019 | 0.1754509400 | 0.3430464637 |
| 50 | 100 | 0.5000000000 | 0.4038298286 | 0.5961701714 |
| 5 | 20 | 0.2902825314 | 0.1118600528 | 0.4687050100 |
| 0 | 10 | 0.1387700844 | 0.0000000000 | 0.2775401688 |
| 10 | 10 | 0.8612299156 | 0.7224598312 | 1.0000000000 |

### Dataset: Purity/ploidy-corrected VAF (CNAqc inversion m·CCF = VAF·(2(1−π)+π·n_tot)/π)

**Source:** CNAqc (Genome Biology 2024) eq.; Tarabichi 2017 diploid case.

| VAF | purity π | ploidy n_tot | adjusted (m·CCF) |
|-----|----------|--------------|-------------------|
| 0.40 | 0.80 | 2 | 1.00 |
| 0.35 | 0.70 | 2 | 1.00 |
| 0.20 | 0.50 | 2 | 0.80 |
| 0.30 | 0.50 | 4 | 1.80 |

---

## Assumptions

1. **ASSUMPTION: Wilson z = 1.96 for the default 95% interval.** The source states z₀.₀₅ = 1.96 for 95%; we use 1.96 verbatim (not the more precise 1.959964) to match the cited value. This is source-backed, not a free assumption; documented here for traceability.
2. **ASSUMPTION: `AdjustVAFForPurity` normal copy number = 2 and assumes the observed VAF is from a heterozygous-eligible diploid-normal background.** The CNAqc/Tarabichi formulas fix normal ploidy at 2 (autosomal diploid). Sex chromosomes / non-diploid normals are out of scope.

---

## Recommendations for Test Coverage

1. **MUST Test:** Empirical VAF = altReads/totalReads for representative ratios (0.25, 0.50, 1.0, 0.0). — Evidence: GATK AD-based definition.
2. **MUST Test:** Wilson 95% interval center/lower/upper against the extracted formula (exact values above). — Evidence: Wikipedia/Wilson 1927.
3. **MUST Test:** Wilson no-overshoot at p̂=0 (lower = 0) and p̂=1 (upper = 1). — Evidence: Wikipedia (overshoot avoided).
4. **MUST Test:** `AdjustVAFForPurity` diploid heterozygous: VAF 0.4, π 0.8 ⇒ 1.0; VAF 0.2, π 0.5, ploidy 2 ⇒ 0.8. — Evidence: CNAqc / Tarabichi.
5. **MUST Test:** totalReads = 0 ⇒ VAF defined as 0 (no coverage); altReads > totalReads ⇒ ArgumentOutOfRangeException; negative counts ⇒ exception. — Evidence: GATK corner cases.
6. **SHOULD Test:** purity = 0 ⇒ ArgumentOutOfRangeException (division by zero). — Rationale: CNAqc denominator division by π.
7. **COULD Test:** invariant 0 ≤ VAF ≤ 1 and lower ≤ center ≤ upper over a sweep. — Rationale: bounds invariants.

---

## References

1. Wilson, E.B. (1927). Probable inference, the law of succession, and statistical inference. Journal of the American Statistical Association 22(158):209–212. https://doi.org/10.1080/01621459.1927.10502953 (formula retrieved via Wikipedia: https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval)
2. Benjamin D, Sato T, Cibulskis K, et al. / GATK team. FAQ for Mutect2; AlleleFraction. Broad Institute GATK documentation. https://gatk.broadinstitute.org/hc/en-us/articles/360050722212-FAQ-for-Mutect2
3. Tarabichi M, Salcedo A, Deshwar AG, et al. (2017/2021). Principles of Reconstructing the Subclonal Architecture of Cancers. Cold Spring Harb Perspect Med. https://pmc.ncbi.nlm.nih.gov/articles/PMC5538405/
4. Househam J, Cresswell GD, et al. / Caravagna G (2024). Computational validation of clonal and subclonal copy number alterations from bulk tumor sequencing using CNAqc. Genome Biology 25:38. https://doi.org/10.1186/s13059-024-03170-5

---

## Change History

- **2026-06-14**: Initial documentation.
