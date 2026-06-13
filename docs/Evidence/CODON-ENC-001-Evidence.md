# Evidence Artifact: CODON-ENC-001

**Test Unit ID:** CODON-ENC-001
**Algorithm:** Effective Number of Codons (ENC / Nc), Wright 1990
**Date Collected:** 2026-06-13

---

## Online Sources

### Fuglsang A. (2004) "The 'effective number of codons' revisited" (Biochem. Biophys. Res. Commun. 317:957–964)

**URL:** https://eclass.uoa.gr/modules/document/file.php/D473/%CE%92%CE%B9%CE%B2%CE%BB%CE%B9%CE%BF%CE%B3%CF%81%CE%B1%CF%86%CE%AF%CE%B1/DNA%20Composition/Fuglsang_2004.pdf
**Accessed:** 2026-06-13
**How retrieved:** WebSearch query `Wright 1990 "effective number of codons" Nc formula gene 87` → fetched the PDF with WebFetch (binary saved locally) → extracted full text with `pdftotext`.
**Authority rank:** 1 (peer-reviewed paper that reproduces Wright's original equations verbatim)

**Key Extracted Points:**

1. **Codon homozygosity (Eq. 1, verbatim):** `F̂ = ( n·Σ_{i=1..k} p_i² − 1 ) / ( n − 1 )`, "where n is the total count for the amino acid in the gene, and p_i is the codon frequency for the ith synonymous codon for the particular amino acid." (k = number of synonymous codons for that amino acid.)
2. **Per-amino-acid effective codons (Eq. 2):** `N̂c(aa) = 1 / F̂_aa`.
3. **Gene-level aggregation (Eq. 3, verbatim):** `N̂c = 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆`, where `F̂₂` is the average homozygosity for amino acids with degeneracy two (histidine, glutamine, etc.) "and so on".
4. **Missing-amino-acid rule (Eq. 4):** if a gene lacks (e.g.) threonine, then `F̂₄` is the average of the codon homozygosities of the remaining four-fold amino acids: `F̂₄ = (F̂_pro + F̂_gly + F̂_ala + F̂_val) / 4` — i.e. average the estimable F values within the same degeneracy class.
5. **Upper-bound rule:** "There is a chance that N̂c, calculated through use of Eq. (3), will exceed 61. In that case, Wright recommends re-adjusting the result down to 61."
6. **Isoleucine (3-fold) fallback:** "Wright suggested … using `F̂₃ = (F̂₂ + F̂₄)/2` when the isoleucine estimator was not possible to calculate." (Eq. 5a.)
7. **Range:** "In extremely biased genes the effective number of codons can approach 20, while in unbiased genes it will approach 61."
8. **Calculability constraint:** "Nc can only be calculated if there are at least two codons for each amino acid with synonymous codons in the sequence" — homozygosity needs n ≥ 2 (the (n−1) denominator is undefined for n ≤ 1).
9. **Reference organism:** "In his original paper, Wright used Escherichia coli K12 as a reference organism."

### Fuglsang A. (2006) "Estimating the 'effective number of codons': the Wright way…" (Genetics 172(2):1301–1307)

**URL:** https://academic.oup.com/genetics/article/172/2/1301/5923091
**Accessed:** 2026-06-13
**How retrieved:** WebSearch query `effective number of codons Wright degeneracy classes "nine" twofold "five" fourfold "three" sixfold isoleucine Met Trp standard genetic code` → fetched with WebFetch.
**Authority rank:** 1 (peer-reviewed, Genetics)

**Key Extracted Points:**

1. **F formula confirmed:** `F = (n·Σ p_i² − 1)/(n−1)` with `p_i = n_i/n`; Wright uses sampling **without** replacement, which the paper concludes is the superior estimator.
2. **Constraints confirmed:** Nc is constrained between 20 and 61; absent / insufficient amino acids handled by within-class averaging; isoleucine has special 3-fold handling.

### Standard genetic code degeneracy classes (via WebSearch)

**How retrieved:** WebSearch (same query as above) returned the standard-genetic-code degeneracy partition.
**Authority rank:** 2 (NCBI standard genetic code, Table 1)

**Key Extracted Points:**

1. **Degeneracy partition of the 20 amino acids + stop:** "5 quartets (4 codons each), 9 doublets (2 codons each), 3 sextets (6 codons each), 1 triplet (3 codons) and 2 singlets (1 codon each)." The two singlets are Met (ATG) and Trp (TGG); the single triplet is isoleucine (ATT/ATC/ATA). This is exactly the partition assumed by Eq. (3): the constant `2` = the two single-codon amino acids (Met + Trp), and `9, 1, 5, 3` count the two-, three-, four- and six-fold amino acids. Stop codons are excluded from Nc.

---

## Documented Corner Cases and Failure Modes

### From Fuglsang (2004)

1. **Amino acid with n ≤ 1:** F̂ is undefined (denominator n−1 → 0); such an amino acid cannot contribute its own F. Nc requires ≥ 2 codons for each represented amino acid; the within-class average (Eq. 4) is the prescribed substitute when a class member is missing.
2. **Empty degeneracy class:** if no amino acid in a class has an estimable F (e.g. no isoleucine at all), the class average is undefined. For the 3-fold class Wright gives the explicit fallback `F̂₃ = (F̂₂ + F̂₄)/2` (Eq. 5a).
3. **Overshoot past 61:** Eq. (3) can yield N̂c > 61 for nearly-uniform short genes; re-adjust the final value down to 61.
4. **Homozygosity becomes zero:** if all codons of an amino acid are used equally and counts are small, F̂ can be very low, making the apparent per-aa Nc exceed the degeneracy. The gene-level value is still re-adjusted to ≤ 61.

---

## Test Datasets

### Dataset: Fully unbiased gene (asymptotic limit)

**Source:** Fuglsang (2004), §"The effective number of codons" — unbiased gene gives F̂₂ = 0.5, F̂₄ = 0.25 (and by the same reasoning F̂₃ = 1/3, F̂₆ = 1/6).

| Class | F̂ (unbiased) | contribution to Nc |
|-------|---------------|--------------------|
| 2-fold | 0.5 | 9 / 0.5 = 18 |
| 3-fold (Ile) | 1/3 | 1 / (1/3) = 3 |
| 4-fold | 0.25 | 5 / 0.25 = 20 |
| 6-fold | 1/6 | 3 / (1/6) = 18 |
| singlets (Met,Trp) | — | 2 |
| **Total** | | **61** |

### Dataset: Fuglsang (2004) "no bias discrepancy" simulation (Nc = 40.5)

**Source:** Fuglsang (2004), Simulation results / Table 4: "all twofold degenerate aa have Nc = 1.5, isoleucine Nc = 2, the fourfold degenerate aa have Nc = 2.5, and the sixfold degenerate have Nc = 3.5 (these numbers sum up to 40.5 effective codons)."

| Class | per-aa Nc | count | class Nc total |
|-------|-----------|-------|----------------|
| 2-fold | 1.5 | 9 | 13.5 |
| 3-fold | 2.0 | 1 | 2.0 |
| 4-fold | 2.5 | 5 | 12.5 |
| 6-fold | 3.5 | 3 | 10.5 |
| singlets | 1.0 | 2 | 2.0 |
| **Total** | | | **40.5** |

### Dataset: Hand-derived two-fold example (exact F by Eq. 1)

Phe codons only, TTT × 3 and TTC × 1 (n = 4, p = (3/4, 1/4)):
Σp² = 9/16 + 1/16 = 10/16 = 0.625. F̂ = (4·0.625 − 1)/(4 − 1) = (2.5 − 1)/3 = 0.5. Nc(Phe) = 1/0.5 = 2.

Perfectly even Phe, TTT × 2 and TTC × 2 (n = 4): Σp² = 0.5, F̂ = (4·0.5 − 1)/3 = 1/3, Nc(Phe) = 3 (> 2, illustrating the per-aa overshoot corner case).

---

## Assumptions

1. **ASSUMPTION: Lower clamp at 20.** Wright/Fuglsang state Nc *approaches* 20 in extreme bias and explicitly prescribe re-adjusting **down to 61** at the top. They do not prescribe a hard clamp at 20; 20 is the structural minimum (every degeneracy class collapses to one codon ⇒ Nc(aa)=1). Retaining `Math.Max(20, …)` is consistent with the stated range and cannot raise a legitimately-computed value, but it is not an explicit Wright instruction. Treated as a defensive bound, not an algorithmic parameter.

---

## Recommendations for Test Coverage

1. **MUST Test:** Fully unbiased gene → Nc = 61 (use codon counts equal within every amino acid, large enough n that the asymptotic F values are reached, then the value re-adjusts to exactly 61). — Evidence: Fuglsang (2004) unbiased dataset.
2. **MUST Test:** Single-amino-acid two-fold gene with exact hand-derived F (TTT×3,TTC×1 ⇒ F=0.5 ⇒ Nc(Phe)=2; missing classes filled by Eq.4 within-class average — verify against an explicit derivation). — Evidence: Eq. 1, Eq. 2, Eq. 4.
3. **MUST Test:** Maximally biased gene (one codon per amino acid) → Nc = 20. — Evidence: range statement.
4. **MUST Test:** Invariant 20 ≤ Nc ≤ 61 for arbitrary inputs. — Evidence: range statement + Eq. 3 re-adjustment.
5. **MUST Test:** Isoleucine-absent gene uses F̂₃ = (F̂₂+F̂₄)/2 fallback. — Evidence: Eq. 5a.
6. **MUST Test:** null DnaSequence → ArgumentNullException; empty/whitespace string → 0. — Evidence: contract (no Wright rule for empty gene; degenerate input).
7. **SHOULD Test:** Lowercase input normalized; invalid (non-ACGT) codons skipped. — Rationale: parser robustness consistent with sibling methods.
8. **COULD Test:** Overshoot case (near-uniform short gene) re-adjusts to 61. — Rationale: Eq. 3 overshoot corner case.

---

## References

1. Wright, F. (1990). The 'effective number of codons' used in a gene. *Gene* 87(1):23–29. https://doi.org/10.1016/0378-1119(90)90491-9 (original; equations reproduced verbatim in refs 2 & 3 below, which were the retrieved sources).
2. Fuglsang, A. (2004). The 'effective number of codons' revisited. *Biochem. Biophys. Res. Commun.* 317(3):957–964. https://doi.org/10.1016/j.bbrc.2004.03.138 (retrieved PDF: https://eclass.uoa.gr/modules/document/file.php/D473/%CE%92%CE%B9%CE%B2%CE%BB%CE%B9%CE%BF%CE%B3%CF%81%CE%B1%CF%86%CE%AF%CE%B1/DNA%20Composition/Fuglsang_2004.pdf).
3. Fuglsang, A. (2006). Estimating the 'effective number of codons': the Wright way of determining codon homozygosity leads to superior estimates. *Genetics* 172(2):1301–1307. https://academic.oup.com/genetics/article/172/2/1301/5923091

---

## Change History

- **2026-06-13**: Initial documentation.
