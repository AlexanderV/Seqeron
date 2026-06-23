# Evidence Artifact: ONCO-CHIP-001

**Test Unit ID:** ONCO-CHIP-001
**Algorithm:** Clonal Hematopoiesis (CHIP) Filtering for cfDNA Liquid Biopsy
**Date Collected:** 2026-06-15

---

## Online Sources

### Steensma et al. (2015) — CHIP definition (Blood)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC4624443/ (PMC mirror of *Blood* 126(1):9–16)
**Retrieved how:** WebSearch query `Steensma 2015 Blood clonal hematopoiesis of indeterminate potential CHIP definition VAF 2% somatic mutation`, then WebFetch of the PMC article URL (the ashpublications.org HTML returned HTTP 403; the PMC mirror was opened instead).
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed paper that coined and formally defined the CHIP term)

**Key Extracted Points:**

1. **VAF threshold:** Quoted from the article: "the mutant allele fraction must be ≥2% in the peripheral blood." Rationale quoted: "with deep enough sequencing, a mutation can be found in every individual, and current outcomes data are based on a minimum variant allele fraction of >2%." ⇒ CHIP requires VAF ≥ 0.02.
2. **Somatic-mutation-in-driver-gene requirement:** "detectable somatic clonal mutations in genes recurrently mutated in hematologic malignancies."
3. **Absence-of-malignancy requirement:** CHIP applies to those who "lack a known hematologic malignancy or other clonal disorder," including patients "who do not meet diagnostic criteria for MDS, as well as those with normal peripheral blood counts."
4. **Driver genes named:** DNMT3A, TET2, ASXL1, TP53, SF3B1, JAK2, PPM1D, BCORL1, GNAS, and others ("the 19 genes most commonly mutated in healthy older adults", Figure 2A).

### Genovese et al. (2014) — Blood-cancer risk inferred from blood DNA (NEJM)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC4290021/ (PMC mirror of *N Engl J Med* 371(26):2477–2487)
**Retrieved how:** WebSearch query `Genovese 2014 NEJM clonal hematopoiesis somatic mutations blood cancer DNMT3A TET2 ASXL1 variant allele fraction`, then WebFetch of the PMC article URL.
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed; 12,380-person whole-exome study establishing CH driver genes)

**Key Extracted Points:**

1. **Recurrent CH driver genes:** Quoted: "Four genes (DNMT3A, TET2, ASXL1, and PPM1D) had disproportionately high numbers of somatic mutations." "DNMT3A had the most observed mutations (190), followed by ASXL1 (35) and TET2 (31)." Additional recurrent variants: JAK2 V617F (24 participants) and SF3B1 K700E (9 participants).
2. **Allelic-fraction detection model:** Quoted: "Assuming that a somatic mutation would be present in only a subset of the cells … we predicted that the mutant allele would be present in less than 50% of the sequencing reads." (Somatic CH variants are sub-clonal, VAF < 0.5.)

### Razavi et al. (2019) — Sources of plasma cfDNA variants (Nature Medicine)

**URL:** https://www.nature.com/articles/s41591-019-0652-7 (paywalled HTML redirected to an auth host; facts taken from the indexing summary + the secondary sources below that quote it)
**Retrieved how:** WebSearch query `Razavi 2019 Nature Medicine clonal hematopoiesis confounds cell-free DNA white blood cell matched filtering CHIP`; the publisher HTML 303-redirected to `idp.nature.com`, so the quantitative findings were taken from the search index abstract and corroborated by the Arango-Argoty 2025 paper (which cites Razavi) and the matched-WBC search below.
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed; *Nat. Med.* 25:1928–1937)

**Key Extracted Points:**

1. **CH is a major cfDNA confounder:** 81.6% of cfDNA mutations in controls and 53.2% in cancer patients had features consistent with clonal hematopoiesis.
2. **Matched-WBC design:** high-intensity sequencing of cfDNA AND matched white-blood-cell (WBC) DNA (508 genes; ~2 Mb; >60,000× raw depth); matched cfDNA–WBC sequencing is required for accurate variant interpretation — a cfDNA variant also present in matched WBC is WBC/CH-derived, not tumor.

### Arango-Argoty et al. (2025) — AI model for CH variants in cfDNA (npj Precision Oncology)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC12092662/ (*NPJ Precis Oncol* 9:147, doi:10.1038/s41698-025-00921-w)
**Retrieved how:** WebSearch query `cfDNA variant also present in matched white blood cell buffy coat classified clonal hematopoiesis removed filter liquid biopsy Razavi method`, then WebFetch of the PMC article URL.
**Accessed:** 2026-06-15
**Authority rank:** 3 (peer-reviewed reference method that operationalizes the matched-WBC rule and cites Razavi 2019 as its basis)

**Key Extracted Points:**

1. **Matched-WBC rule (gold standard):** Quoted: "Studies characterizing the presence of variants in blood plasma cfDNA commonly utilize matched sequencing of a nucleated blood cell fraction (e.g. peripheral blood mononuclear cells or white blood cells) to determine variant origin (tumor or CH)." Razavi et al. (2019) is cited as establishing this matched-WBC approach.
2. **Most prevalent CH genes:** Quoted: "all variants in the most prevalent genes in clonal hematopoietic cells (DNMT3A, TET2, and ASXL1) were removed." (Confirms the canonical top-3.)
3. **VAF caveat:** "The VAF of key hematologic cancer driver genes has also been suggested as way to identify CH variants in cfDNA samples, but the exact relationship between VAF and variant origin remains unclear." (⇒ the gene+VAF heuristic flags *candidate* CHIP; the matched-WBC subtraction is the definitive origin test.)

### Bolton et al. (2020) — Strict matched-WBC origin rule (Nature Genetics)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7891089/ (PMC mirror of *Nat Genet* 52(11):1219–1226, doi:10.1038/s41588-020-00710-0)
**Retrieved how:** WebSearch query `Bolton 2020 Nature Genetics clonal hematopoiesis CH-PD origin assignment matched blood somatic tumor`, then WebFetch of the PMC article URL (the Methods origin-assignment rule was extracted from the retrieved text).
**Accessed:** 2026-06-23
**Authority rank:** 1 (peer-reviewed; 24,439-patient paired blood–tumour MSK-IMPACT study)

**Key Extracted Points:**

1. **CH definition (VAF + read floor):** Quoted from the Methods: clonal-hematopoiesis mutations were defined as those with "a variant allele fraction of at least 2% and at least 10 supporting reads." ⇒ a confident WBC/CH call requires WBC VAF ≥ 0.02 AND ≥ 10 supporting reads.
2. **Blood-to-tumour VAF fold rule (origin assignment):** Quoted: variants detected in blood with "a VAF of at least twice that in the tumor or 1.5 times the VAF if the tumor biopsy site was a lymph node were considered" CH/somatic-blood. ⇒ WBC VAF ≥ 2× tumour VAF (1.5× for lymph node) ⇒ WBC/CH origin; otherwise tumour-derived.
3. **Rationale for the fold ratio:** Quoted: the ratio "was chosen … through simulations of leukocyte contamination in the tumor" (minimising sensitivity/specificity loss of CH calls). ⇒ the fold thresholds are a sourced, not invented, value; they are exposed as parameters.

---

## Documented Corner Cases and Failure Modes

### From Steensma et al. (2015)

1. **Sub-2% VAF mutations:** A driver-gene somatic mutation below VAF 2% does NOT meet the CHIP definition (below the formal threshold and the reliable detection limit).
2. **Driver mutation present but disease present:** if diagnostic criteria for a hematologic malignancy are met, it is not CHIP (out of scope here — the algorithm assays non-diagnostic blood/plasma).

### From Razavi et al. (2019) / Arango-Argoty et al. (2025)

3. **Tumor variant masquerading in cfDNA:** without matched-WBC subtraction, CH variants are mis-called as tumor — the dominant false-positive source (81.6%/53.2%).
4. **Variant absent from matched WBC:** a cfDNA variant NOT seen in matched WBC is retained as a candidate tumor variant even if it falls in a CHIP gene.

---

## Test Datasets

### Dataset: Canonical CHIP driver-gene set (Steensma 2015 / Genovese 2014)

**Source:** Steensma et al. (2015) Fig. 2A; Genovese et al. (2014) recurrent genes.

| Parameter | Value |
|-----------|-------|
| Canonical CHIP genes | DNMT3A, TET2, ASXL1, TP53, JAK2, SF3B1, SRSF2, PPM1D |
| CHIP VAF threshold | 0.02 (≥ 2%) |
| Top-3 by prevalence | DNMT3A, TET2, ASXL1 |

### Dataset: Worked classification cases (derived from the cited definition/threshold)

**Source:** derived from Steensma (2015) VAF ≥ 2% + driver-gene rule; Razavi (2019) matched-WBC rule.

| # | Gene | VAF | In matched WBC? | Expected IdentifyCHIPVariants | Expected FilterCHIP (kept?) |
|---|------|-----|-----------------|-------------------------------|-----------------------------|
| 1 | DNMT3A | 0.05 | — | CHIP (driver gene, VAF ≥ 0.02) | (gene/VAF mode) removed |
| 2 | DNMT3A | 0.01 | — | not CHIP (VAF < 0.02) | kept |
| 3 | EGFR | 0.30 | — | not CHIP (not a CHIP gene) | kept |
| 4 | EGFR | 0.30 | yes | n/a (gene rule) | removed (WBC-matched ⇒ CH) |
| 5 | TP53 | 0.40 | no | CHIP candidate (driver gene) | removed (gene+VAF heuristic, rule b) |

### Dataset: Strict matched-WBC origin calls (Bolton 2020 rule)

**Source:** Bolton et al. (2020) — WBC VAF ≥ 0.02 AND ≥ 10 reads AND WBC VAF ≥ φ × tumour VAF (φ = 2.0; 1.5 lymph node).

| # | Tumour VAF | WBC VAF | WBC reads | Fold φ | Expected CallVariantOrigin |
|---|-----------|---------|-----------|--------|----------------------------|
| 1 | 0.10 | 0.30 | 40 | 2.0 | Chip (0.30 ≥ 2%, 40 ≥ 10, 0.30 ≥ 2×0.10) |
| 2 | 0.40 | — (absent) | — | 2.0 | Tumor (no WBC evidence) |
| 3 | 0.10 | 0.30 | 9 | 2.0 | Tumor (9 < 10 reads) |
| 4 | 0.01 | 0.02 | 10 | 2.0 | Chip (all boundaries inclusive: 0.02=2%, 0.02=2×0.01, 10=10) |
| 5 | 0.005 | 0.015 | 50 | 2.0 | Tumor (0.015 < 0.02 WBC-VAF minimum) |
| 6 | 0.30 | 0.40 | 40 | 2.0 | Tumor (0.40 < 2×0.30=0.60) |
| 7 | 0.25 | 0.40 | 40 | 1.5 | Chip (0.40 ≥ 1.5×0.25=0.375); Tumor under default 2.0 |

> **Note on row 5 (origin caveat).** Under a *strict matched-WBC-only* origin call (Razavi 2019), a CH-gene variant **absent** from the matched WBC would be retained as a candidate tumour variant — the matched-WBC presence test, not gene identity, is the definitive origin test, and "the exact relationship between VAF and variant origin remains unclear" (Arango-Argoty 2025). `FilterCHIP` is deliberately **conservative**: rule (b) is a labelled gene+VAF *heuristic fallback* (see `Clonal_Hematopoiesis_Filtering.md` §2.2) that removes a CH-driver-gene variant meeting VAF ≥ τ even with no WBC evidence, exactly as rows 1 and 5 show. The heuristic over-removes relative to the strict matched-WBC definition; callers who want the strict rule pass an empty/custom `chipGenes` panel so only matched-WBC subtraction applies. The earlier "kept" entry for row 5 contradicted both row 1 and the documented contract and was corrected on 2026-06-16.

---

## Assumptions

1. **ASSUMPTION: Canonical default gene set** — When the caller does not supply a CHIP-gene panel, the default is the source-cited set {DNMT3A, TET2, ASXL1, TP53, JAK2, SF3B1, SRSF2, PPM1D} (Steensma 2015 Fig 2A; Genovese 2014). This is a *labeled* canonical set, not an invented value; callers may override it. Gene-symbol comparison is case-insensitive (HGNC symbols are upper-case).
2. **ASSUMPTION: Matched-WBC presence test** — "Present in matched WBC" is decided by `IsVariantDetected` style alt-read evidence at the same locus (≥ 1 alt read), matching the repository's existing MRD detection convention (Wan 2020) and the matched-WBC subtraction principle (Razavi 2019). The exact universal alt-read cutoff is assay-specific and is configurable.

---

## Recommendations for Test Coverage

1. **MUST Test:** `IdentifyCHIPVariants` flags a driver-gene variant at VAF ≥ 0.02 as CHIP and a sub-2% one as not-CHIP — Evidence: Steensma 2015 (VAF ≥ 2%).
2. **MUST Test:** `IdentifyCHIPVariants` does not flag a non-CHIP-gene variant regardless of VAF — Evidence: Steensma 2015 (driver-gene requirement).
3. **MUST Test:** `IdentifyCHIPVariants` honors a caller-supplied gene panel — Evidence: Razavi/Arango-Argoty (caller-supplied panel) + Framework status.
4. **MUST Test:** `FilterCHIP` removes a cfDNA variant present in matched WBC (even a non-CHIP-gene one) — Evidence: Razavi 2019 matched-WBC rule.
5. **MUST Test:** `FilterCHIP` retains a cfDNA variant absent from matched WBC — Evidence: Razavi 2019.
6. **MUST Test:** VAF exactly at 0.02 boundary is CHIP (≥, inclusive) — Evidence: Steensma 2015 ("≥2%").
6b. **MUST Test:** `CallVariantOrigin` calls CHIP when WBC VAF ≥ 2%, ≥ 10 reads, and WBC VAF ≥ 2× tumour VAF — Evidence: Bolton 2020.
6c. **MUST Test:** `CallVariantOrigin` calls tumour when the locus is absent from matched WBC (even a CH driver gene) — Evidence: Bolton 2020 (no WBC evidence ⇒ tumour).
6d. **MUST Test:** `CallVariantOrigin` read floor (9 < 10 ⇒ tumour) and VAF floor (< 2% ⇒ tumour) and fold rule (< 2× ⇒ tumour); lymph-node 1.5× changes the call — Evidence: Bolton 2020.
7. **SHOULD Test:** case-insensitive gene matching — Rationale: HGNC symbols upper-case but inputs vary.
8. **SHOULD Test:** null / empty inputs throw / return empty — Rationale: repository convention.
9. **COULD Test:** order preservation of kept variants — Rationale: deterministic output contract.

---

## References

1. Steensma DP, Bejar R, Jaiswal S, Lindsley RC, Sekeres MA, Hasserjian RP, Ebert BL. (2015). Clonal hematopoiesis of indeterminate potential and its distinction from myelodysplastic syndromes. *Blood* 126(1):9–16. https://doi.org/10.1182/blood-2015-03-631747 (PMC4624443)
2. Genovese G, Kähler AK, Handsaker RE, et al. (2014). Clonal Hematopoiesis and Blood-Cancer Risk Inferred from Blood DNA Sequence. *N Engl J Med* 371(26):2477–2487. https://doi.org/10.1056/NEJMoa1409405 (PMC4290021)
3. Razavi P, Li BT, Brown DN, et al. (2019). High-intensity sequencing reveals the sources of plasma circulating cell-free DNA variants. *Nat Med* 25:1928–1937. https://doi.org/10.1038/s41591-019-0652-7
4. Arango-Argoty G, et al. (2025). An artificial intelligence-based model for prediction of clonal hematopoiesis variants in cell-free DNA samples. *NPJ Precis Oncol* 9:147. https://doi.org/10.1038/s41698-025-00921-w (PMC12092662)
5. Wan JCM, Heider K, Gale D, et al. (2020). ctDNA monitoring using patient-specific sequencing and integration of variant reads. *Sci Transl Med* 12(548):eaaz8084. https://doi.org/10.1126/scitranslmed.aaz8084
6. Bolton KL, Ptashkin RN, Gao T, et al. (2020). Cancer therapy shapes the fitness landscape of clonal hematopoiesis. *Nat Genet* 52(11):1219–1226. https://doi.org/10.1038/s41588-020-00710-0 (PMC7891089)

---

## Change History

- **2026-06-15**: Initial documentation.
- **2026-06-23**: Added Bolton et al. (2020) strict matched-WBC origin rule (WBC VAF ≥ 2%, ≥ 10 reads, ≥ 2×/1.5× tumour VAF) for the new `CallVariantOrigin` method; added the strict-origin dataset and MUST tests 6b–6d.
