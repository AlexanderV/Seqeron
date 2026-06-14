# Evidence Artifact: ONCO-MSI-001

**Test Unit ID:** ONCO-MSI-001
**Algorithm:** Microsatellite Instability (MSI) Detection — fraction-of-unstable-loci scoring and status classification
**Date Collected:** 2026-06-14

---

## Online Sources

### MSIsensor2 — reference implementation README (niu-lab)

**URL:** https://raw.githubusercontent.com/niu-lab/msisensor2/master/README.md (also https://github.com/niu-lab/msisensor2)
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference implementation in an established bioinformatics tool)
**Retrieved by:** WebSearch "MSIsensor microsatellite instability score percentage unstable loci threshold 20% MSI-H" → WebFetch of the GitHub page and the raw README.md.

**Key Extracted Points:**

1. **MSI score definition:** Verbatim — "the msi score (number of msi sites / all valid sites) can be calculated." The score is the fraction of unstable (somatic/"msi") microsatellite loci among all valid evaluated loci, expressed as a percentage.
2. **MSI-High cutoff:** Verbatim — "The recommended msi score cutoff value is 20% (msi high: msi score >= 20%)." The boundary is inclusive: a sample is MSI-High when msi score ≥ 20%.

### MSIsensor — Niu et al. (2014), Bioinformatics 30(7):1015–1016 (Oxford Academic)

**URL:** https://academic.oup.com/bioinformatics/article/30/7/1015/236553
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper)
**Retrieved by:** WebSearch "Niu MSIsensor 2014 Bioinformatics proportion of unstable microsatellites somatic score percentage definition" → WebFetch of the Oxford Academic article page.

**Key Extracted Points:**

1. **Per-site test:** Each microsatellite site is tested by a chi-square test comparing tumor vs normal repeat-length (allele-frequency) distributions; a site is called unstable/somatic when the difference is significant under a default FDR threshold of 0.05.
2. **MSI score = percentage of somatic sites:** The final MSIsensor MSI score is "the percentage of microsatellite sites with a somatic indel" (percentage of unstable sites among evaluated sites).
3. **Practical decision boundary (this dataset):** "Among 71 MSI samples, 70 have an MSI score >3.5. In addition, 165 of 168 MSS samples have a score <3.5" — i.e. MSIsensor's original cohort separated MSI vs MSS near an MSI score of 3.5%.

### Boland et al. (1998) — NCI Workshop on Microsatellite Instability (PubMed)

**URL:** https://pubmed.ncbi.nlm.nih.gov/9823339/ (Boland CR et al., Cancer Res 58(22):5248–5257)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed consensus paper / NCI workshop)
**Retrieved by:** WebSearch "Bethesda panel microsatellite instability MSI-H two or more markers unstable MSI-L MSS NCI workshop 1998" → WebFetch of the PubMed abstract page.

**Key Extracted Points:**

1. **MSI-H:** Verbatim — a tumor is MSI-H "if two or more of the five markers show instability (i.e., have insertion/deletion mutations)" (≥ 2 of 5 unstable markers).
2. **MSI-L:** Verbatim — "if only one of the five markers shows instability" (exactly 1 of 5 unstable).
3. **MSS:** No marker shows instability (0 of 5 unstable). The abstract notes MSS vs MSI-L can only be reliably distinguished with a larger panel.
4. **Reference panel:** "A panel of five microsatellites has been validated and is recommended as a reference panel" (the Bethesda panel; classic markers BAT-25, BAT-26, D2S123, D5S346, D17S250).

### Geiersbach 2014 / British Journal of Cancer 2014;111:813 — revised Bethesda fractions (search snippet, primary paywalled)

**URL:** https://www.nature.com/articles/bjc2014167
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed); full text behind an IDP login redirect — only the search-result snippet was retrievable, NOT the full article.
**Retrieved by:** WebSearch (snippet) "Bethesda panel ... NCI workshop 1998"; WebFetch returned a 303 redirect to an authentication endpoint, so the full text could not be opened.

**Key Extracted Points (from the retrievable snippet only):**

1. **Fraction form of the 2004 revised criteria:** MSI-H = "≥2 out of 5 markers are unstable or ≥40%"; MSI-L = "1 out of 5 ... or ≥20% and <40%"; MSS = "no unstable markers detected." This is consistent with Boland 1998 (≥2/5, 1/5, 0/5) and is used only as a cross-check, not as the sole authority.

---

## Documented Corner Cases and Failure Modes

### From Boland et al. (1998)

1. **MSS vs MSI-L ambiguity:** With the 5-marker panel, distinguishing true MSS (0/5) from MSI-L (1/5) is unreliable; a larger panel is recommended. Classification at the 1-marker level is inherently low-confidence.

### From MSIsensor / MSIsensor2

1. **Insufficient valid loci:** The msi score is `msi sites / all valid sites`; if no microsatellite locus has sufficient coverage there are zero valid sites and the score is undefined (division by zero).
2. **Tumor-only mode:** MSIsensor2 computes the score without a matched normal; the fraction-of-unstable-loci definition is unchanged.

---

## Test Datasets

### Dataset: MSIsensor2 worked cutoff

**Source:** niu-lab/msisensor2 README (accessed 2026-06-14)

| Parameter | Value |
|-----------|-------|
| MSI score formula | unstable loci / valid loci (as %) |
| MSI-H cutoff | ≥ 20% (inclusive) |
| Example: 5 unstable / 25 valid | 20% → MSI-H |
| Example: 4 unstable / 25 valid | 16% → not MSI-H |
| Example: 0 unstable / 25 valid | 0% → MSS-range |

### Dataset: Bethesda 5-marker panel (Boland 1998)

**Source:** Boland et al. (1998), Cancer Res 58:5248–5257

| Unstable / total markers | Status |
|--------------------------|--------|
| 0 / 5 | MSS |
| 1 / 5 | MSI-L |
| 2 / 5 | MSI-H |
| 3 / 5 | MSI-H |
| 5 / 5 | MSI-H |

---

## Assumptions

1. **ASSUMPTION: MSI-L band for the computational fraction score** — MSIsensor2 defines only a binary MSI-H cutoff (≥20%); it does not define an MSI-L band on the continuous score. The categorical Bethesda classification (MSS/MSI-L/MSI-H) is therefore applied to the **discrete marker-count** input (Boland 1998), and the continuous fraction score is classified only as the source-backed binary MSI-H (≥20%) vs not-high. No MSI-L band is invented for the continuous score.

---

## Recommendations for Test Coverage

1. **MUST Test:** Continuous MSI score = unstable/valid as a fraction, and MSI-H classification at the 20% boundary (≥20% High, <20% not High). — Evidence: MSIsensor2 README.
2. **MUST Test:** Bethesda categorical classification by marker count: 0→MSS, 1→MSI-L, ≥2→MSI-H over a 5-marker panel. — Evidence: Boland 1998.
3. **MUST Test:** Score invariant 0 ≤ score ≤ 1 and boundary exactness (20% inclusive). — Evidence: MSIsensor2 README; checklist invariant.
4. **SHOULD Test:** Zero valid loci → score undefined (throws). — Rationale: division by zero in `unstable/valid`.
5. **SHOULD Test:** unstable > valid, negative counts → invalid input (throws). — Rationale: counts must satisfy 0 ≤ unstable ≤ valid.
6. **COULD Test:** Larger panels at the exact 40% Bethesda fraction map to MSI-H. — Rationale: cross-check with BJC 2014 fraction form.

---

## References

1. Niu B, Ye K, Zhang Q, Lu C, Xie M, McLellan MD, Wendl MC, Ding L (2014). MSIsensor: microsatellite instability detection using paired tumor-normal sequence data. Bioinformatics 30(7):1015–1016. https://doi.org/10.1093/bioinformatics/btt755 (retrieved via https://academic.oup.com/bioinformatics/article/30/7/1015/236553)
2. niu-lab/msisensor2 — Microsatellite instability (MSI) detection for tumor only data. README, accessed 2026-06-14. https://github.com/niu-lab/msisensor2
3. Boland CR, Thibodeau SN, Hamilton SR, Sidransky D, Eshleman JR, Burt RW, Meltzer SJ, Rodriguez-Bigas MA, Fodde R, Ranzani GN, Srivastava S (1998). A National Cancer Institute Workshop on Microsatellite Instability for cancer detection and familial predisposition: development of international criteria for the determination of microsatellite instability in colorectal cancer. Cancer Res 58(22):5248–5257. https://pubmed.ncbi.nlm.nih.gov/9823339/
4. British Journal of Cancer (2014) 111:813 (revised Bethesda fraction form) — search-result snippet only; full text not retrievable (login redirect). https://www.nature.com/articles/bjc2014167

---

## Change History

- **2026-06-14**: Initial documentation (ONCO-MSI-001).
