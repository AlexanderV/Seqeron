# Evidence Artifact: ONCO-ARTIFACT-001

**Test Unit ID:** ONCO-ARTIFACT-001
**Algorithm:** Sequencing Artifact Detection (OxoG / FFPE deamination substitution classification + strand-orientation bias)
**Date Collected:** 2026-06-14

---

## Online Sources

### Chen L. et al. (2017) — "DNA damage is a pervasive cause of sequencing errors, directly confounding variant identification" (Science) + GenomeWeb / Nature Methods summaries

**URL:** https://www.science.org/doi/10.1126/science.aai8690 (primary; paywalled, abstract retrieved); GIV definition retrieved from https://github.com/Ettwiller/Damage-estimator and the Nature Methods news write-up https://www.nature.com/articles/nmeth.4254 (both opened this session) and the comment article https://pmc.ncbi.nlm.nih.gov/articles/PMC7350422/
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper) + 3 (the authors' own reference implementation, `Ettwiller/Damage-estimator`)

**Retrieval method:** WebSearch query `Chen 2017 GIV score formula G_T over C_A ratio damage estimator damaged 1.5 undamaged 1 paired-end`; WebFetch of the Damage-estimator GitHub repo and the PMC7350422 comment article.

**Key Extracted Points:**

1. **OxoG substitution signature:** Oxidative damage such as 8-oxo-dG (8-oxoguanine) "leads to an excess of G-to-T transversion errors when R1 is mapped to a reference genome, whereas R2 reads will show an excess of the reverse complement of G-to-T — i.e., C-to-A errors." So the OxoG artifact class is **G>T (read 1) / C>A (read 2)** — a single substitution class G:C>T:A with a read-orientation imbalance.
2. **GIV (Global Imbalance Value) definition:** The GIV score relies on the imbalance between variants detected in read 1 and read 2 of a paired-end run; a **separate score is calculated for each of the 12 possible mutation types**. For OxoG it is the **global imbalance in the number of G>T variants in R1 compared with R2**, i.e. GIV_G_T = (count of G>T in R1) / (count of G>T in R2).
3. **GIV interpretation / threshold:** "A GIV score of 1 indicates there is no DNA damage and a GIV score above 1.5 is defined as damaged DNA." (Nature Methods summary, retrieved). The comment article PMC7350422 states GIV_G_T = 2 means "the G>T error rate in the 8-oxoG mode is twice the error rate of the non–8-oxoG mode", and GIV_G_T > 5 indicates significant damage — confirming GIV is a ratio with neutral value 1.
4. **Damage-estimator output:** `estimate_damage.pl` output columns include "[1] raw count of variant type [2] variant type (ex. G_T, G to T) ... [6] GIV-score", confirming GIV is computed per substitution type from R1/R2 counts.

### NEB / Ettwiller Damage-estimator (reference implementation)

**URL:** https://github.com/Ettwiller/Damage-estimator
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference implementation by the paper's authors)

**Retrieval method:** WebFetch of the repository landing page and README.

**Key Extracted Points:**

1. **Mechanism:** "Damage estimation is based on the systematic mutation rate difference between the first in pair and the second in pair reads." Confirms the R1-vs-R2 imbalance is the basis.
2. **Reference value:** "If you have followed the standard protocol for acoustic shearing during library preparation you should obtain a GIV score for G_T around 2." Confirms the ratio's neutral baseline is 1 and an elevated value (≈2) signals OxoG damage.

### FFPE cytosine deamination artifact (Do & Dobrovic 2015, *Clinical Chemistry* review, via ScienceDirect)

**URL:** https://www.sciencedirect.com/science/article/pii/S152515781630188X ("Deamination Effects in Formalin-Fixed, Paraffin-Embedded Tissue Samples in the Era of Precision Medicine"); cross-checked with Oxford NAR review https://academic.oup.com/nar/article/51/14/7143/7205768
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed review)

**Retrieval method:** WebSearch query `FFPE formalin fixation deamination cytosine C to T G to A artifact substitution signature sequencing`; WebFetch of the ScienceDirect article.

**Key Extracted Points:**

1. **FFPE deamination substitution signature:** "A base substitution of C to T (C>T) or G to A (G>A, from C>T on the antisense strand), collectively **C:G>T:A**, is one of the most common artifacts, and it is caused by ... cytosine deamination." So the FFPE artifact class is **C>T / G>A**.
2. **Mechanism:** "Deaminated cytosine results in uracil, which pairs with adenine instead of guanine" → C>T; on the complementary strand this reads as G>A.
3. **Distinct from oxidation:** "The two most prevalent artifact types in FFPE-extracted DNA ... are C>T/G>A caused by cytosine deamination and C>A/G>T that mostly results from base oxidation." Confirms the two artifact classes are disjoint by substitution type.

### GATK FisherStrand / StrandBiasTest (Broad Institute) — strand-bias detection

**URL:** https://gatk.broadinstitute.org/hc/en-us/articles/360035532152-Fisher-s-Exact-Test (FS docs, retrieved via WebSearch summary); source code https://raw.githubusercontent.com/broadinstitute/gatk/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/annotator/StrandBiasTest.java and FisherStrand.java
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference implementation, canonical variant-calling toolkit)

**Retrieval method:** WebSearch query for FisherStrand FS; WebFetch of the two raw GitHub source files.

**Key Extracted Points:**

1. **2×2 contingency table:** The FisherStrand test "creates a 2x2 contingency table of #ref alleles on + strand, #ref alleles on - strand, #alt alleles on + strand, and #alt alleles on - strand." Source confirms cell ordering: `array[0]=ref forward, array[1]=ref reverse, array[2]=alt forward, array[3]=alt reverse` (`StrandBiasTest`, ARRAY_DIM=2, ARRAY_SIZE=4).
2. **Statistic:** FS is the **Phred-scaled p-value of a two-sided Fisher's exact test** on that table: `FS = -10 · log10(p)` (`QualityUtils.phredScaleErrorRate(max(pValue, MIN_PVALUE))`, FisherStrand.java).
3. **Floor:** `MIN_PVALUE = 1E-320` caps the smallest p-value before Phred scaling to avoid infinite FS.
4. **Interpretation:** Null hypothesis = no difference in ref/alt distribution across strands (no strand bias). Higher FS ⇒ stronger strand bias ⇒ more likely artifact. (No strand bias ⇒ p ≈ 1 ⇒ FS ≈ 0.)

---

## Documented Corner Cases and Failure Modes

### From GATK FisherStrand

1. **No strand bias:** when ref and alt alleles are evenly distributed across forward/reverse strands, the Fisher p-value ≈ 1 and FS ≈ 0.
2. **Perfect strand segregation:** when all alt reads are on one strand and all ref reads on the other, the p-value is minimal and FS is large.
3. **Empty / single-orientation table:** a table with a zero margin (e.g. no reverse reads at all) yields p = 1 (no evidence of bias) under the two-sided test.

### From Chen 2017 / Damage-estimator

1. **No R2 G>T counts (zero denominator):** GIV is a ratio of R1/R2 counts; an undefined ratio (0 G>T in R2) must be handled (treated as no imbalance evidence when both are 0; as maximal imbalance when only R2 is 0).
2. **Balanced R1/R2:** GIV = 1 ⇒ undamaged.

### From FFPE review

1. **Substitution not C>T/G>A and not G>T/C>A:** is neither a deamination nor an oxidation artifact by substitution class (e.g. A>G transition) ⇒ not flagged by the substitution-class rule.

---

## Test Datasets

### Dataset: OxoG GIV worked example (Chen 2017 / Damage-estimator)

**Source:** Chen et al. 2017; Damage-estimator README (GIV ≈ 2 for standard shearing; >1.5 = damaged; 1 = undamaged).

| Parameter | Value |
|-----------|-------|
| Substitution | G>T |
| R1 G>T count | 200 |
| R2 G>T count | 100 |
| GIV_G_T = 200/100 | 2.0 (damaged; > 1.5) |
| Balanced R1=R2=100 | GIV = 1.0 (undamaged) |

### Dataset: FisherStrand contingency worked example (GATK)

**Source:** GATK FisherStrand / two-sided Fisher exact test on the 2×2 table [ref_fwd, ref_rev, alt_fwd, alt_rev].

| Parameter | Value |
|-----------|-------|
| Balanced table [10,10,10,10] | two-sided p = 1.0 ⇒ FS = 0.000 |
| Segregated table [20,0,0,20] | two-sided p = (computed exactly) ⇒ FS large (> 0) |

(Exact p for [20,0,0,20] is derived in the test by the hypergeometric two-sided Fisher formula; FS = -10·log10(p).)

### Dataset: FFPE / OxoG substitution classification (substitution class only)

**Source:** Do & Dobrovic 2015 review; Chen 2017.

| ref→alt | Class |
|---------|-------|
| C>T | FFPE deamination |
| G>A | FFPE deamination |
| G>T | OxoG oxidation |
| C>A | OxoG oxidation |
| A>G | neither (not an artifact class) |

---

## Assumptions

1. **ASSUMPTION: No BAM parser.** The checklist signature is `FilterArtifacts(variants, bamFile)`. The repository has no BAM reader; per the unit note ("implement the retrievable rule-based artifact classification"), the read-orientation/strand evidence that a BAM would supply (per-strand and per-read-mate alt/ref counts) is passed directly on the variant observation record instead of parsed from a file. This is an API-shape decision, not a correctness-affecting one — the classification rules (substitution class, GIV ratio, Fisher strand p) are unchanged.
2. **ASSUMPTION: GIV neutral/decision thresholds.** GIV = 1 (undamaged) and GIV > 1.5 (damaged) are taken verbatim from the Nature Methods summary of Chen 2017; the underlying ratio (R1 G>T / R2 G>T) is from the paper and Damage-estimator. The 1.5 cutoff is a documented operational threshold, not invented.

---

## Recommendations for Test Coverage

1. **MUST Test:** C>T and G>A substitutions classify as FFPE deamination. — Evidence: Do & Dobrovic 2015 (C:G>T:A).
2. **MUST Test:** G>T and C>A substitutions classify as OxoG oxidation. — Evidence: Chen 2017 (G>T R1 / C>A R2).
3. **MUST Test:** A>G (and other) substitutions classify as not-an-artifact by class. — Evidence: FFPE/oxidation classes are specific.
4. **MUST Test:** GIV_G_T = R1/R2; R1=200,R2=100 ⇒ 2.0 (damaged > 1.5); R1=R2 ⇒ 1.0 (undamaged). — Evidence: Chen 2017 / Damage-estimator.
5. **MUST Test:** FisherStrand FS = -10·log10(two-sided Fisher p); balanced table [10,10,10,10] ⇒ FS = 0.0; segregated table ⇒ FS > 0. — Evidence: GATK FisherStrand.
6. **MUST Test:** `FilterArtifacts` removes flagged artifacts and keeps real variants; result ⊆ input. — Evidence: composition of the above rules.
7. **SHOULD Test:** GIV with zero R2 count handled (no division error). — Rationale: documented corner case.
8. **SHOULD Test:** null / empty inputs throw / return empty. — Rationale: API contract per sibling methods.
9. **COULD Test:** strand-bias p-value monotonicity (more segregation ⇒ higher FS). — Rationale: invariant check.

---

## References

1. Chen L., Liu P., Evans T.C. Jr., Ettwiller L.M. (2017). DNA damage is a pervasive cause of sequencing errors, directly confounding variant identification. Science 355(6326):752–756. https://www.science.org/doi/10.1126/science.aai8690
2. Ettwiller L. Damage-estimator (reference implementation, GIV score). https://github.com/Ettwiller/Damage-estimator
3. Eberle M. (news summary) / Nature Methods (2017). DNA variants or DNA damage? Nat Methods 14:330. https://www.nature.com/articles/nmeth.4254
4. Comment on "DNA damage is a pervasive cause of sequencing errors" (2018). Science / PMC7350422 (GIV_G_T interpretation). https://pmc.ncbi.nlm.nih.gov/articles/PMC7350422/
5. Do H., Dobrovic A. (2015). Sequence Artifacts in DNA from Formalin-Fixed Tissues / Deamination Effects in FFPE Tissue Samples. Clin Chem / ScienceDirect. https://www.sciencedirect.com/science/article/pii/S152515781630188X
6. Broad Institute. GATK FisherStrand (FS) — Fisher's Exact Test for strand bias. https://gatk.broadinstitute.org/hc/en-us/articles/360035532152-Fisher-s-Exact-Test ; source: https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/annotator/StrandBiasTest.java

---

## Change History

- **2026-06-14**: Initial documentation.
