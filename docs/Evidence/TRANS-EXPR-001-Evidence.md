# Evidence Artifact: TRANS-EXPR-001

**Test Unit ID:** TRANS-EXPR-001
**Algorithm:** Expression Quantification (TPM, FPKM/RPKM, Quantile Normalization)
**Date Collected:** 2026-06-13

---

## Online Sources

### Wagner, Kin & Lynch (2012) — TPM introduction (Theory in Biosciences)

**URL:** https://link.springer.com/article/10.1007/s12064-012-0162-3 (DOI: 10.1007/s12064-012-0162-3)
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed paper)
**Retrieved how:** WebSearch query `Wagner Kin Lynch 2012 "Measurement of mRNA abundance using RNA-seq" Theory in Biosciences TPM equation full text`; landing/metadata page opened on Springer Nature Link.

**Key Extracted Points:**

1. **Purpose of TPM:** The paper proposes TPM (transcripts per million) as a modification of RPKM that eliminates inconsistency in measuring RNA abundance among samples. RPKM is intended to measure relative molar RNA concentration (rmc); the average rmc is constant for a set of transcripts (the inverse of the number of transcripts mapped), but RPKM does not respect this invariance property.
2. **Invariant (sum / average):** The average TPM equals 10^6 divided by the number of annotated transcripts, and is therefore a constant; equivalently, TPM values across a sample sum to exactly 10^6.

### "Misuse of RPKM or TPM normalization…" — Zhao, Ye & Stanton (2020), RNA journal (PMC7373998)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7373998/
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed paper)
**Retrieved how:** WebFetch of the PMC article; prompt requested the verbatim TPM and RPKM equations and the sum-to-million invariant.

**Key Extracted Points:**

1. **RPKM formula (verbatim from page):** `RPKM = 10^9 * Reads mapped to the transcript / (Total reads * Transcript length)`.
2. **TPM formula (verbatim from page):** `TPM = 10^6 * (reads mapped to transcript / transcript length) / Sum(reads mapped to transcript / transcript length)`.
3. **TPM–RPKM relationship (verbatim):** `TPM = 10^6 * RPKM / Sum(RPKM)`.
4. **Invariant:** "The average TPM is equal to 10^6 divided by the number of annotated transcripts in a given annotation, and thus is a constant" → TPM sums to one million within a sample.

### Pimentel, H. (2014) — "What the FPKM? A review of RNA-Seq expression units"

**URL:** https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
**Accessed:** 2026-06-13
**Authority rank:** 3 (widely-cited review by a kallisto/sleuth author; used to corroborate the primary formulas, not as sole authority)
**Retrieved how:** WebFetch; prompt requested the exact TPM and FPKM equations with variable definitions.

**Key Extracted Points:**

1. **TPM (verbatim form):** `TPM_i = (X_i / l_i) * (1 / Σ_j (X_j / l_j)) * 10^6`, where `X_i` = read/fragment counts for feature i, `l_i` = (effective) length of feature i, sum over all features j.
2. **FPKM (verbatim form):** `FPKM_i = X_i / ((l_i/10^3) * (N/10^6)) = (X_i / (l_i * N)) * 10^9`, where `N` = total sequenced fragments/reads in the sample.
3. **Relationship:** `TPM_i = (FPKM_i / Σ_j FPKM_j) * 10^6` — TPM is FPKM rescaled to sum to one million.

### Wikipedia — "Quantile normalization" (citing Bolstad et al. 2003)

**URL:** https://en.wikipedia.org/wiki/Quantile_normalization
**Accessed:** 2026-06-13
**Authority rank:** 4 (Wikipedia article citing the primary Bolstad 2003 method); worked example and tie rule used.
**Retrieved how:** WebFetch; prompt requested the full worked-example matrices (input, sorted, rank means, final) and the exact tie-handling rule.

**Key Extracted Points:**

1. **Algorithm:** To quantile-normalize distributions to each other without a reference, sort each column; set each rank to the mean (arithmetic) of the values at that rank across columns; then re-place the rank means at the original positions. The highest value in each column becomes the mean of the highest values, etc.
2. **Worked-example input (3 columns of 4 values):** C1 = (5,2,3,4); C2 = (4,1,4,2); C3 = (3,4,6,8).
3. **Sorted columns:** C1 = (2,3,4,5); C2 = (1,2,4,4); C3 = (3,4,6,8).
4. **Rank means (lowest→highest rank):** (2+1+3)/3 = 2.00; (3+2+4)/3 = 3.00; (4+4+6)/3 = 4.666…; (5+4+8)/3 = 5.666….
5. **Tie rule (verbatim):** "When values are tied in rank, they should instead be assigned the mean of the values corresponding to the ranks they would normally represent if they were different." In C2 the two tied `4`s span ranks iii and iv → each gets (4.666…+5.666…)/2 = 5.166….
6. **Final normalized matrix (verbatim, row-major rows A-D, read column-wise):** C1 = (5.67, 2.00, 3.00, 4.67); C2 = (5.17, 2.00, 5.17, 3.00) — BOTH tied `4`s (rows A and C) become 5.17; C3 = (2.00, 3.00, 4.67, 5.67). Source LaTeX matrix: rows `5.67 5.17 2.00 / 2.00 2.00 3.00 / 3.00 5.17 4.67 / 4.67 3.00 5.67`.

---

## Documented Corner Cases and Failure Modes

### From Zhao, Ye & Stanton (2020) / Wagner (2012)

1. **All-zero counts:** If every transcript has zero reads, Σ(reads/length) = 0 and TPM is `0/0`; the quantity is undefined. Resolution adopted (see Assumptions): emit TPM = 0 for every gene (the denominator collapses, no proportion is defined).
2. **Zero gene length:** `reads/length` requires length > 0; a length of 0 makes RPK undefined. Inputs must have positive length; non-positive length yields FPKM = 0 and is excluded from a meaningful RPK contribution.
3. **Cross-sample comparison caveat:** Both RPKM and TPM are within-sample relative measures; the paper warns they are misused when compared across samples/protocols. This is an interpretation caveat, not an output rule.

### From Wikipedia (quantile normalization)

1. **Tied ranks:** values tied within a column receive the average of the rank means they would otherwise span (see Source point 5).
2. **Empty input / zero genes:** undefined for an empty matrix; no rank means exist.

---

## Test Datasets

### Dataset: TPM three-gene derivation (derived from the cited formula)

**Source:** Computed from the verbatim TPM formula in Zhao/Ye/Stanton (2020) and Pimentel (2014): `TPM_i = (X_i/l_i)/Σ(X_j/l_j) * 10^6`.

| Gene | Count X | Length l | RPK = X/l | TPM |
|------|---------|----------|-----------|-----|
| A | 10 | 2000 | 0.005 | 125000 |
| B | 20 | 4000 | 0.005 | 125000 |
| C | 30 | 1000 | 0.030 | 750000 |
| **Σ** | 60 | — | 0.040 | **1000000** |

Derivation: ΣRPK = 0.04; TPM_A = 0.005/0.04·10^6 = 125000; TPM_C = 0.03/0.04·10^6 = 750000; total = 10^6 (invariant).

### Dataset: FPKM single-gene derivation

**Source:** `FPKM_i = X_i·10^9/(l_i·N)` (Zhao/Ye/Stanton 2020; Pimentel 2014).

| Parameter | Value |
|-----------|-------|
| Count X | 1000 |
| Length l | 2000 |
| Total reads N | 1,000,000 |
| FPKM = 1000·10^9/(2000·10^6) | 500 |

### Dataset: Quantile normalization worked example

**Source:** Wikipedia "Quantile normalization" (Bolstad et al. 2003). Input columns C1=(5,2,3,4), C2=(4,1,4,2), C3=(3,4,6,8).

| Position (row) | C1 in | C2 in | C3 in | C1 out | C2 out | C3 out |
|----------------|-------|-------|-------|--------|--------|--------|
| A (0) | 5 | 4 | 3 | 5.666… | 5.166… | 2.000 |
| B (1) | 2 | 1 | 4 | 2.000 | 2.000 | 3.000 |
| C (2) | 3 | 4 | 6 | 3.000 | 5.166… | 4.666… |
| D (3) | 4 | 2 | 8 | 4.666… | 3.000 | 5.666… |

Rank means: r0=2.0, r1=3.0, r2=14/3=4.666…, r3=17/3=5.666…; tie mean (r2,r3)=31/6=5.166…. Both tied `4`s in C2 (rows A and C) receive 5.166… (verbatim final matrix shows 5.17 for both).

---

## Assumptions

1. **ASSUMPTION: All-zero TPM denominator → 0.** No retrieved source specifies the output when Σ(reads/length)=0 (every count zero). The mathematically defined quantity is 0/0 (undefined). The implementation emits TPM=0 for all genes. This is a degenerate-input convention, not a value from literature; it does not affect any non-degenerate result because every non-zero-denominator case is fully defined by the cited formula.
2. **ASSUMPTION: Effective length = annotated length.** The cited formulas use effective length `l̃_i` (Pimentel 2014); the library's public API takes the annotated transcript length directly (no fragment-length correction). For the within-sample TPM/FPKM ratios this is the standard substitution `l̃_i = l_i` and does not change the formula's structure; it only omits an optional refinement.

---

## Recommendations for Test Coverage

1. **MUST Test:** TPM three-gene example returns (125000, 125000, 750000) exactly and sums to 10^6 — Evidence: TPM formula, Zhao/Ye/Stanton (2020) + Pimentel (2014).
2. **MUST Test:** TPM values always sum to 10^6 for any non-degenerate input (invariant) — Evidence: Wagner (2012); Zhao/Ye/Stanton (2020).
3. **MUST Test:** FPKM single-gene = 500 for X=1000, l=2000, N=10^6 — Evidence: FPKM formula.
4. **MUST Test:** Quantile normalization reproduces the Wikipedia worked-example output matrix exactly (including tie handling on C2) — Evidence: Wikipedia/Bolstad 2003.
5. **SHOULD Test:** Empty inputs yield empty output for all three methods — Rationale: documented degenerate case.
6. **SHOULD Test:** Non-positive length / total reads → FPKM = 0 — Rationale: documented failure mode.
7. **COULD Test:** All-zero counts → all TPM = 0 — Rationale: degenerate denominator convention.

---

## References

1. Wagner GP, Kin K, Lynch VJ. 2012. Measurement of mRNA abundance using RNA-seq data: RPKM measure is inconsistent among samples. Theory in Biosciences 131(4):281–285. https://doi.org/10.1007/s12064-012-0162-3
2. Zhao S, Ye Z, Stanton R. 2020. Misuse of RPKM or TPM normalization when comparing across samples and sequencing protocols. RNA 26(8):903–909. https://pmc.ncbi.nlm.nih.gov/articles/PMC7373998/
3. Mortazavi A, Williams BA, McCue K, Schaeffer L, Wold B. 2008. Mapping and quantifying mammalian transcriptomes by RNA-Seq. Nature Methods 5(7):621–628. https://doi.org/10.1038/nmeth.1226
4. Pimentel H. 2014. What the FPKM? A review of RNA-Seq expression units. https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
5. Bolstad BM, Irizarry RA, Astrand M, Speed TP. 2003. A comparison of normalization methods for high density oligonucleotide array data based on variance and bias. Bioinformatics 19(2):185–193. (via) https://en.wikipedia.org/wiki/Quantile_normalization

---

## Change History

- **2026-06-13**: Initial documentation.
