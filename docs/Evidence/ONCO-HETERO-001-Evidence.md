# Evidence Artifact: ONCO-HETERO-001

**Test Unit ID:** ONCO-HETERO-001
**Algorithm:** Tumor Heterogeneity Analysis (MATH score, Shannon clonal diversity, subclone count, subclonal fraction)
**Date Collected:** 2026-06-15

---

## Online Sources

### Mroz EA, Rocco JW (2013) — MATH, a novel measure of intratumor genetic heterogeneity (primary)

**URL:** https://pubmed.ncbi.nlm.nih.gov/23079694/ (Oral Oncology 49(3):211–215)
**Retrieved by:** WebSearch "Mroz Rocco 2013 MATH score mutant-allele tumor heterogeneity median absolute deviation" → WebFetch of the PubMed page.
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed primary paper)

**Key Extracted Points:**

1. **MATH formula:** Fetched text states verbatim **"MATH = 100 * MAD/median"**. The numerator is the Median Absolute Deviation (MAD) of the mutant-allele fractions (distribution width); the denominator is the median of the mutant-allele fractions (distribution centre); the ratio is scaled by 100.
2. **Definition:** MATH is "the ratio of the width to the center of its distribution of mutant-alleles", computed from NGS mutant-allele fractions across mutated loci within a single tumour. Higher MATH = more intratumour genetic heterogeneity.

### Mroz EA et al. (2015) — Intra-tumor Genetic Heterogeneity and Mortality in Head and Neck Cancer (PLOS Medicine)

**URL:** https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1001786
**Retrieved by:** WebSearch (MATH search above) → WebFetch of the PLOS Medicine article.
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed primary paper)

**Key Extracted Points:**

1. **MAD scaling 1.4826:** Methods state **"The median [of absolute deviations] is then multiplied by a factor of 1.4826, so that the expected MAD of a normally distributed variable is equal to its SD."**
2. **Full procedure:** absolute difference of each MAF from the median → median of those absolute differences → × 1.4826 → divide by median MAF → × 100 (percentage). Confirms **MATH = 100 × MAD/median**.

### maftools `mathScore.R` (reference implementation)

**URL:** https://raw.githubusercontent.com/PoisonAlien/maftools/master/R/mathScore.R
**Retrieved by:** WebSearch "maftools math.score source code mad 1.4826 median vaf 100 mathScore.R" → WebFetch of the raw GitHub source.
**Accessed:** 2026-06-15
**Authority rank:** 3 (established bioinformatics library, Bioconductor)

**Key Extracted Points:**

1. **Exact code:** `abs.med.dev = abs(vaf - median(vaf))`; `pat.mad = median(abs.med.dev) * 100`; `pat.math = pat.mad * 1.4826 / median(vaf)`.
2. **Algebraic equivalence:** MATH = 100 × 1.4826 × median(|VAF − median(VAF)|) / median(VAF) — identical to the primary papers.
3. **VAF normalisation:** if max VAF > 1, divide VAFs by 100 (assumes percent input). This implementation requires VAFs already in [0, 1], so no rescaling.

### Liu Z et al. (2017) — Quantification of within-sample genetic heterogeneity from SNP-array data (Shannon-based ITH)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC5468233/ (BMC Genomics 18:457)
**Retrieved by:** WebSearch "intratumor heterogeneity Shannon diversity index clone fractions" → WebFetch of the PMC article.
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed) / 4 (Shannon primary cited within)

**Key Extracted Points:**

1. **Shannon formula:** Fetched text gives **S = −∑ᵢⁿ pᵢ × ln(pᵢ)** using the **natural logarithm**; "Shannon diversity indices were calculated using the clonal frequencies of each mixture".
2. **Richness:** richness is defined as **the number of clones present** (monoclonal = 1 clone; polyclonal = ≥ 2; mixtures of 1–5 profiles).

### Shannon CE (1948) diversity/entropy definition (corroborating)

**URL:** WebSearch "Shannon 1948 diversity index entropy H = -sum p_i ln p_i definition natural logarithm ecology"
**Retrieved by:** WebSearch (results summary).
**Accessed:** 2026-06-15
**Authority rank:** 4 (Wikipedia-citing-primary tier; corroborates the primary above)

**Key Extracted Points:**

1. **Definition:** H' = −Σ(pᵢ ln pᵢ), R = richness (number of classes), pᵢ = fraction of the i-th class; base e (natural log) usual in ecology; changing the base multiplies H by a constant.
2. **Behavior:** H increases with both richness and evenness; H → 0 as the number of subclones decreases or one dominates.

### Landau DA et al. (2013) — clonal/subclonal CCF threshold (reused from ONCO-CLONAL-001)

**URL:** existing repository Evidence/code citation; threshold P(CCF > 0.95) > 0.5 ⇒ clonal.
**Retrieved by:** repository `OncologyAnalyzer.cs` constant `ClonalCcfThreshold = 0.95` (Landau et al. 2013, Cell 152(4):714–726).
**Accessed:** 2026-06-15
**Authority rank:** 1 (already-validated primary in ONCO-CLONAL-001)

**Key Extracted Points:**

1. **Subclonal definition:** a mutation is subclonal when its CCF is below 0.95 (the clonal CCF threshold). Used here to compute the fraction of subclonal mutations.

---

## Documented Corner Cases and Failure Modes

### From maftools `mathScore.R` / Mroz papers

1. **Zero median MAF:** MATH = 100·MAD/median divides by the median; a median of 0 makes MATH undefined (division by zero). Must be rejected/guarded.
2. **All identical VAFs:** MAD = 0 ⇒ MATH = 0 (no heterogeneity), a valid minimal value.
3. **Single mutation:** median = that value, MAD = 0 ⇒ MATH = 0.

### From Liu et al. (2017) / Shannon

1. **Single clone:** richness 1, p = 1, H = −1·ln 1 = 0 (minimum diversity).
2. **Even split:** k clones of equal size give H = ln k (maximum for that richness).

---

## Test Datasets

### Dataset: MATH worked example (odd count)

**Source:** Mroz & Rocco (2013) formula MATH = 100·1.4826·median(|f−median(f)|)/median(f), hand-derived.

| Parameter | Value |
|-----------|-------|
| VAFs | 0.10, 0.20, 0.30, 0.40, 0.50 |
| median(VAF) | 0.30 |
| abs deviations | 0.20, 0.10, 0.00, 0.10, 0.20 |
| median(abs dev) (raw MAD) | 0.10 |
| scaled MAD (×1.4826) | 0.14826 |
| MATH = 100·0.14826/0.30 | 49.42 |

### Dataset: MATH worked example (even count)

**Source:** same formula, hand-derived.

| Parameter | Value |
|-----------|-------|
| VAFs | 0.20, 0.40, 0.60, 0.80 |
| median(VAF) | 0.50 |
| abs deviations | 0.30, 0.10, 0.10, 0.30 |
| median(abs dev) (raw MAD) | 0.20 |
| scaled MAD (×1.4826) | 0.29652 |
| MATH = 100·0.29652/0.50 | 59.304 |

### Dataset: Shannon diversity worked examples

**Source:** Liu et al. (2017) / Shannon (1948), H = −Σ pᵢ ln pᵢ (natural log), hand-derived.

| Clone fractions | H |
|-----------------|---|
| {1.0} | 0.0 |
| {0.5, 0.5} | −ln 0.5 = 0.6931471805599453 |
| {0.25, 0.25, 0.25, 0.25} | ln 4 = 1.3862943611198906 |

---

## Assumptions

1. **ASSUMPTION: Clone fractions for Shannon = mutation proportions per CCF cluster** — Liu et al. (2017) compute Shannon over "clonal frequencies"; here clone fractions pᵢ are the proportion of mutations assigned to each CCF cluster (ONCO-CCF-001 clustering). This is a standard operationalisation but the exact pᵢ source (cluster CCF vs. cluster size) is a modelling choice; using cluster sizes (counts) is the natural per-mutation diversity reading.
2. **ASSUMPTION: Median definition for even counts** — R's `median` (used by maftools) averages the two central order statistics; replicated here. Sources do not enumerate an even-count rule explicitly, so the standard R behaviour is adopted.

---

## Recommendations for Test Coverage

1. **MUST Test:** MATH on odd/even worked examples equals 49.42 / 59.304 exactly — Evidence: Mroz & Rocco (2013); maftools `mathScore.R`.
2. **MUST Test:** MATH = 0 when all VAFs identical (MAD = 0) — Evidence: formula.
3. **MUST Test:** Shannon H = ln k for k equal clones; H = 0 for a single clone — Evidence: Liu et al. (2017); Shannon (1948).
4. **MUST Test:** subclone count = number of occupied CCF clusters — Evidence: Liu et al. (2017) richness.
5. **MUST Test:** subclonal fraction = #(CCF < 0.95)/n — Evidence: Landau et al. (2013).
6. **SHOULD Test:** zero-median VAF throws; null/empty/out-of-range throw — Rationale: documented division-by-zero failure mode.
7. **COULD Test:** invariant MATH ≥ 0 over varied inputs — Rationale: registry invariant ITH_score ≥ 0.

---

## References

1. Mroz EA, Rocco JW (2013). MATH, a novel measure of intratumor genetic heterogeneity, is high in poor-outcome classes of head and neck squamous cell carcinoma. Oral Oncology 49(3):211–215. https://pubmed.ncbi.nlm.nih.gov/23079694/
2. Mroz EA, Tward AD, Hammon RJ, Ren Y, Rocco JW (2015). Intra-tumor genetic heterogeneity and mortality in head and neck cancer: analysis of data from The Cancer Genome Atlas. PLOS Medicine 12(2):e1001786. https://doi.org/10.1371/journal.pmed.1001786
3. Mayakonda A et al. maftools `mathScore.R`. https://github.com/PoisonAlien/maftools/blob/master/R/mathScore.R
4. Liu Z, Zhang S (2017). Quantification of within-sample genetic heterogeneity from SNP-array data. BMC Genomics 18:457 (PMC5468233). https://pmc.ncbi.nlm.nih.gov/articles/PMC5468233/
5. Shannon CE (1948). A mathematical theory of communication. Bell System Technical Journal 27:379–423. https://en.wikipedia.org/wiki/Diversity_index#Shannon_index
6. Landau DA et al. (2013). Evolution and impact of subclonal mutations in chronic lymphocytic leukemia. Cell 152(4):714–726. https://doi.org/10.1016/j.cell.2013.01.019

---

## Change History

- **2026-06-15**: Initial documentation.
