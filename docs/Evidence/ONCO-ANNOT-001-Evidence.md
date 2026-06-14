# Evidence Artifact: ONCO-ANNOT-001

**Test Unit ID:** ONCO-ANNOT-001
**Algorithm:** Cancer-Specific Variant Annotation (AMP/ASCO/CAP 2017 four-tier clinical-significance classification)
**Date Collected:** 2026-06-14

---

## Online Sources

### Li MM et al. (2017) — "Standards and Guidelines for the Interpretation and Reporting of Sequence Variants in Cancer" (J Mol Diagn)

**URL:** https://doi.org/10.1016/j.jmoldx.2016.10.002 (DOI confirmed; Elsevier full text paywalled). Full guideline text retrieved and read this session from an authoritative open copy of the same article (The Journal of Molecular Diagnostics, Vol. 19, No. 1, January 2017, "Li et al"): https://ocpe.mcw.edu/sites/default/files/course/2024-03/AMP-ASCO-CAP%20guidelines%20-%20somatic%20variants.pdf
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed joint consensus guideline / standard)

**Retrieval method:** WebSearch query `Li 2017 AMP ASCO CAP standards interpretation reporting sequence variants cancer tier classification J Mol Diagn`; WebFetch of the MCW-hosted PDF, which was downloaded and converted with `pdftotext -layout` and read directly. Citation metadata (volume 19, issue 1, pages 4–23, DOI 10.1016/j.jmoldx.2016.10.002) cross-confirmed via WebSearch.

**Key Extracted Points:**

1. **Four-tier system (verbatim, abstract):** "tier I, variants with strong clinical significance; tier II, variants with potential clinical significance; tier III, variants of unknown clinical significance; and tier IV, variants deemed benign or likely benign."
2. **Tier ↔ evidence-level mapping (verbatim, body):** "tier I, variants with strong clinical significance (level A and B evidence); tier II, variants with potential clinical significance (level C or D evidence); tier III, variants with unknown clinical significance; and tier IV, variants that are benign or likely benign (Figure 2)."
3. **Four evidence levels (Table 3):** Level A — biomarkers predicting response/resistance to FDA-approved therapies for a specific tumor type, or in professional guidelines. Level B — based on well-powered studies with expert consensus. Level C — FDA/guideline therapies for a *different* tumor type (off-label), clinical-trial inclusion criteria, or diagnostic/prognostic from multiple small studies. Level D — plausible therapeutic significance from preclinical studies, or diagnostic/prognostic from small studies / few case reports without consensus.
4. **Population-frequency benign cutoff (verbatim, Population Databases section):** "There is no standardized cutoff for MAF to be used for eliminating polymorphic or benign variants. In the absence of paired normal tissue, the work group recommends using 1% (0.01) as a primary cutoff, which is also commonly used across many clinical laboratories."
5. **Tier IV population criterion (Table 7):** Population database row reads "MAF [≥] 1% in the general population; or high MAF in some ethnic populations." (The "≥" glyph was lost in PDF extraction but is fixed by point 4's 1% primary cutoff for benign.) Tier IV "FDA-approved therapies, PG, investigational therapies" row = "None".
6. **Tier III population criterion (Table 6):** Population database row = "Absent or extremely low MAF"; FDA/guideline therapies = "Cancer genes: none"; publications = "None or no convincing evidence to determine clinical/biological significance". Figure 2 Tier III box: "Not observed at a significant allele frequency in the general or specific subpopulation databases" AND "No convincing published evidence of cancer association".
7. **Figure 2 Tier IV box (verbatim):** "Observed at significant allele frequency in the general or specific subpopulation databases" / "No existing published evidence of cancer association". This makes population frequency the discriminator between Tier III (rare, cancer association) and Tier IV (common, or no association).
8. **Evidence categories:** the levels apply across three categories — Therapeutic, Diagnosis, Prognosis (Table 3 columns).

### Tate JG et al. (2019) — "COSMIC: the Catalogue Of Somatic Mutations In Cancer" (Nucleic Acids Research)

**URL:** https://academic.oup.com/nar/article/47/D1/D941/5146192 ; DOI https://doi.org/10.1093/nar/gky1015
**Accessed:** 2026-06-14
**Authority rank:** 5 (well-maintained bioinformatics database) / 1 (peer-reviewed database paper)

**Retrieval method:** WebSearch query `COSMIC Catalogue of Somatic Mutations in Cancer Tate 2019 Nucleic Acids Research database somatic variant annotation`.

**Key Extracted Points:**

1. **COSMIC is an external curated database:** "the most detailed and comprehensive resource for exploring the effect of somatic mutations in human cancer", v86 (Aug 2018) holding "almost 6 million coding mutations across 1.4 million tumour samples, curated from over 26,000 publications". It cannot be reproduced/hardcoded in this library — a COSMIC lookup must be performed against caller-supplied records.
2. **AMP/ASCO/CAP use of COSMIC:** Li et al. (2017) Tables 4–6 list "Somatic database: COSMIC, My Cancer Genome, TCGA" as an evidence source, with Tier I "Most likely present", Tier II "Likely present", Tier III "Absent or present without association". COSMIC presence is a somatic-database evidence input, supporting (not by itself determining) the tier.

---

## Documented Corner Cases and Failure Modes

### From Li et al. (2017)

1. **Common polymorphism with no clinical evidence:** MAF ≥ 1% in population databases ⇒ Tier IV benign/likely benign (Table 7).
2. **Clinically significant variant that also appears in population databases:** assigned by evidence level (Tier I for Level A/B), not downgraded by frequency — Tables 4/5 list Tier I/II variants as "Absent or extremely low MAF" but the categorization criterion is the evidence level (Figure 2), and the guideline notes well-studied germline-counterpart variants (TP53, PTEN) may be present in databases yet remain significant.
3. **Rare variant with a cancer association but no clinical evidence level:** Tier III (unknown significance) — distinguished from Tier IV by the presence of a cancer association and absence of a significant population frequency (Figure 2; Table 6).
4. **No cancer association and not common:** Figure 2 Tier IV box lists "No existing published evidence of cancer association" as a Tier IV criterion ⇒ benign/likely benign.

### From COSMIC (Tate 2019)

1. **Variant absent from the supplied catalog:** the lookup must return "not found" (null) rather than fabricating an annotation; COSMIC content is external.

---

## Test Datasets

### Dataset: AMP/ASCO/CAP 2017 tier decision matrix (Li et al. 2017, Figure 2 / Tables 3–7)

**Source:** Li MM et al. (2017), J Mol Diagn 19(1):4–23.

| Evidence level | Population MAF | Cancer association | Expected tier |
|----------------|----------------|--------------------|---------------|
| A | (any) | (any) | Tier I (strong) |
| B | (any) | (any) | Tier I (strong) |
| C | (any) | (any) | Tier II (potential) |
| D | (any) | (any) | Tier II (potential) |
| None | ≥ 0.01 | (any) | Tier IV (benign) |
| None | < 0.01 | false | Tier IV (benign — no cancer association) |
| None | < 0.01 | true | Tier III (unknown) |

### Dataset: Canonical example variants (illustrative, named in Li et al. 2017)

**Source:** Li et al. (2017) — BRAF V600E cited as a Level A/B biomarker (FDA-approved therapy / professional guidelines, diagnostic & prognostic for thyroid carcinoma).

| Gene | Protein change | Evidence level | Population MAF | Cancer assoc. | Expected tier |
|------|----------------|----------------|---------------|---------------|---------------|
| BRAF | p.V600E | A | 0.0 | true | Tier I |
| (rare VUS) | p.X | None | 0.0001 | true | Tier III |
| (common SNP) | p.Y | None | 0.25 | false | Tier IV |

### Dataset: COSMIC lookup (caller-supplied catalog)

**Source:** Tate et al. (2019); Li et al. (2017) somatic-database evidence row.

| Catalog entry | Lookup key | Result |
|---------------|-----------|--------|
| {(BRAF, p.V600E) → "COSV56056643"} | (BRAF, p.V600E) | "COSV56056643" |
| {(BRAF, p.V600E) → "COSV56056643"} | (TP53, p.R175H) | null (not in catalog) |

---

## Assumptions

1. **ASSUMPTION: Caller-supplied evidence inputs.** The AMP/ASCO/CAP guideline classifies variants from external curated knowledge (professional guidelines, population databases, somatic databases, literature). This library does not reproduce those resources; the evidence level, population MAF and cancer-association flag are supplied by the caller (who performed the lookups). This is an input-shape decision, not a correctness-affecting one — the tiering decision rule (Figure 2) is applied verbatim to whatever evidence is supplied.
2. **ASSUMPTION: Tier III vs Tier IV discriminator.** When no clinical evidence level (A–D) is present, Figure 2 distinguishes Tier III ("not observed at significant allele frequency" + "no convincing evidence of cancer association") from Tier IV ("observed at significant allele frequency" or "no published evidence of cancer association"). The implementation uses MAF ≥ 1% OR absence of a cancer association ⇒ Tier IV, otherwise Tier III. This is a direct reading of the Figure 2 boxes and Table 6/7, not invented.

---

## Recommendations for Test Coverage

1. **MUST Test:** Level A and Level B evidence ⇒ Tier I (strong clinical significance). — Evidence: Li et al. (2017) Figure 2 ("level A and B evidence").
2. **MUST Test:** Level C and Level D evidence ⇒ Tier II (potential clinical significance). — Evidence: Li et al. (2017) Figure 2 ("level C or D evidence").
3. **MUST Test:** No evidence level, MAF ≥ 1% ⇒ Tier IV (benign / likely benign). — Evidence: Li et al. (2017) Table 7 + 1% primary cutoff.
4. **MUST Test:** No evidence level, MAF < 1%, no cancer association ⇒ Tier IV (no published cancer association). — Evidence: Li et al. (2017) Figure 2 Tier IV box.
5. **MUST Test:** No evidence level, MAF < 1%, cancer association present ⇒ Tier III (unknown significance). — Evidence: Li et al. (2017) Table 6 / Figure 2 Tier III box.
6. **MUST Test:** Clinical evidence (Level A) takes priority over a high MAF (a significant biomarker is still Tier I). — Evidence: Li et al. (2017) Figure 2 (categorization by evidence level).
7. **MUST Test:** `AnnotateCancerVariants` returns one annotation per input variant, in input order, with correct tiers for a mixed batch. — Evidence: composition of the per-variant rule.
8. **MUST Test:** `GetCOSMICAnnotation` returns the catalog value on a hit and null on a miss. — Evidence: Tate et al. (2019); external-catalog lookup.
9. **MUST Test:** MAF exactly at the 1% boundary (0.01) ⇒ Tier IV (cutoff is "≥ 1%"); just below (0.0099) with cancer association ⇒ Tier III. — Evidence: Li et al. (2017) 1% primary cutoff; "≥" semantics.
10. **SHOULD Test:** Invalid MAF (negative, > 1, NaN) throws `ArgumentOutOfRangeException`; null inputs throw `ArgumentNullException`. — Rationale: API contract per sibling methods.
11. **SHOULD Test:** Empty variant collection ⇒ empty annotation list. — Rationale: boundary.

---

## References

1. Li MM, Datto M, Duncavage EJ, Kulkarni S, Lindeman NI, Roy S, Tsimberidou AM, Vnencak-Jones CL, Wolff DJ, Younes A, Nikiforova MN. (2017). Standards and Guidelines for the Interpretation and Reporting of Sequence Variants in Cancer: A Joint Consensus Recommendation of the Association for Molecular Pathology, American Society of Clinical Oncology, and College of American Pathologists. J Mol Diagn 19(1):4–23. https://doi.org/10.1016/j.jmoldx.2016.10.002 (full text read this session from https://ocpe.mcw.edu/sites/default/files/course/2024-03/AMP-ASCO-CAP%20guidelines%20-%20somatic%20variants.pdf)
2. Tate JG, Bamford S, Jubb HC, et al. (2019). COSMIC: the Catalogue Of Somatic Mutations In Cancer. Nucleic Acids Res 47(D1):D941–D947. https://doi.org/10.1093/nar/gky1015 ; https://academic.oup.com/nar/article/47/D1/D941/5146192

---

## Change History

- **2026-06-14**: Initial documentation.
