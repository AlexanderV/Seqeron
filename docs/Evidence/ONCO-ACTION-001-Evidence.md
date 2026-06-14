# Evidence Artifact: ONCO-ACTION-001

**Test Unit ID:** ONCO-ACTION-001
**Algorithm:** Clinical Actionability Assessment (OncoKB Therapeutic Levels of Evidence)
**Date Collected:** 2026-06-15

---

## Online Sources

### OncoKB: A Precision Oncology Knowledge Base (Chakravarty et al. 2017)

**URL:** https://ascopubs.org/doi/10.1200/PO.17.00011
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed paper, JCO Precision Oncology)
**Retrieved how:** WebSearch query "Chakravarty OncoKB Precision Oncology Knowledge Base levels of evidence JCO Precision Oncology 2017" returned the article landing page; direct WebFetch of the DOI URL returned HTTP 403 (publisher paywall), so the level-of-evidence text was obtained from the official OncoKB Levels-of-Evidence PDF and the OncoKB SOP v3 PDF (below), which restate the same Chakravarty 2017 system.

**Key Extracted Points:**

1. **System purpose:** OncoKB "annotates the biologic and oncogenic effects and prognostic and predictive significance of somatic molecular alterations. Potential treatment implications are stratified by the level of evidence that a specific molecular alteration is predictive of drug response on the basis of US Food and Drug Administration labeling, National Comprehensive Cancer Network guidelines, disease-focused expert group recommendations, and scientific literature." (search-result abstract of PO.17.00011).
2. **Levels are an ordered evidence hierarchy:** treatment implications are stratified by an ordered level of evidence (this is the actionability axis ONCO-ACTION-001 classifies).

### OncoKB Therapeutic Levels of Evidence (official PDF, V2)

**URL:** https://www.oncokb.org/content/files/levelOfEvidence/V2/LevelsOfEvidence.pdf
**Accessed:** 2026-06-15
**Authority rank:** 3 (canonical project documentation / official knowledgebase)
**Retrieved how:** WebSearch query "OncoKB levels of evidence Level 1 Level 2 Level 3A 3B Level 4 R1 R2 definitions therapeutic actionable" surfaced this PDF; downloaded via `curl` (WebFetch saved the 3.8 MB binary locally) and text extracted with `pypdf`.

**Key Extracted Points (verbatim level definitions extracted from the PDF text):**

1. **Level 1 (Standard Care):** "FDA-recognized biomarker predictive of response to an FDA-approved drug in this indication".
2. **Level 2 (Standard Care):** "Standard care biomarker recommended by the NCCN or other professional guidelines predictive of response to an FDA-approved drug in this indication".
3. **Level 3A (Investigational):** "Compelling clinical evidence supports the biomarker as being predictive of response to a drug in this indication".
4. **Level 3B (Investigational):** "Standard care or investigational biomarker predictive of response to an FDA-approved or investigational drug in another indication".
5. **Level 4 (Hypothetical):** "Compelling biological evidence supports the biomarker as being predictive of response to a drug".
6. **Level R1 (Standard Care Resistance):** "Standard care biomarker predictive of resistance to an FDA-approved drug in this indication".
7. **Level R2 (Investigational Resistance):** "Compelling clinical evidence supports the biomarker as being predictive of resistance to a drug".

### OncoKB Curation Standard Operating Procedure v3 (official SOP PDF)

**URL:** https://sop.oncokb.org/static/sop/OncoKB_Curation_Standard_Operating_Procedure_v3.pdf
**Accessed:** 2026-06-15
**Authority rank:** 3 (canonical project documentation)
**Retrieved how:** WebSearch (same query as above) surfaced the SOP; downloaded via `curl` (17 MB) and text extracted with `pypdf`; the relevant passages were located by substring search.

**Key Extracted Points (closely paraphrased from extracted SOP text):**

1. **Highest-vs-investigational grouping:** "The highest levels of evidence, Levels 1 and 2, refer to the standard implications for sensitivity to an FDA-approved drug. Additionally, Level R1 refers to the standard implications for resistance to an FDA-approved drug. Levels 3A, 3B and 4 refer to the investigational implications for sensitivity to either an FDA-approved or investigational drug (in the off-label setting, Level 3B) or an investigational drug (Levels 3A and 4). Level R2 includes investigational implications for resistance to either an FDA-approved or investigational drug."
2. **Provenance / consistency:** the level system "was developed to rank the therapeutic implications associated with an alteration found in a patient tumor sample by the relative weight of the evidence (Chakravarty et al., 2017), and are consistent with the Joint Consensus Recommendation by AMP, ASCO and CAP (Li et al., 2017)".
3. **3A above 3B:** "this system was refined to deprioritize the significance of standard care biomarkers when present in indications outside of the FDA-approved/NCCN listed indication" — i.e. Level 3A (compelling clinical evidence in-indication) ranks above Level 3B (standard care in another indication).

### oncokb-annotator (reference implementation, GitHub)

**URL:** https://raw.githubusercontent.com/oncokb/oncokb-annotator/master/README.md
**Accessed:** 2026-06-15
**Authority rank:** 3 (reference implementation maintained by OncoKB)
**Retrieved how:** WebSearch surfaced the repo; `curl` of the raw README and `grep -nE "LEVEL_R1|LEVEL_1|..."` extracted the column-definition rows.

**Key Extracted Points (verbatim from README annotation-column table):**

1. **HIGHEST_LEVEL order:** "Order: LEVEL_R1 > LEVEL_1 > LEVEL_2 > LEVEL_3A > LEVEL_3B > LEVEL_4 > LEVEL_R2".
2. **HIGHEST_SENSITIVE_LEVEL order:** "Order: LEVEL_1 > LEVEL_2 > LEVEL_3A > LEVEL_3B > LEVEL_4".
3. **HIGHEST_RESISTANCE_LEVEL order:** "Order: LEVEL_R1 > LEVEL_R2".
4. **Highest-level definition:** HIGHEST_LEVEL is "The highest level of evidence for therapeutic implications." — i.e. for a variant with several leveled drug associations, the actionable level is the maximum under the order above.

---

## Documented Corner Cases and Failure Modes

### From oncokb-annotator README

1. **No leveled association (VUS-like):** the annotator HIGHEST_LEVEL columns are populated only for variants that carry a therapeutic implication; a variant with no leveled drug association has no highest level (empty). For ONCO-ACTION-001 this maps to "not actionable" (no level assigned).

### From OncoKB SOP v3

1. **Resistance is not the same axis as sensitivity:** R1/R2 (resistance) and the sensitivity levels are separate implication types; the combined HIGHEST_LEVEL order interleaves them (R1 above 1; R2 below 4) but a separate HIGHEST_SENSITIVE_LEVEL / HIGHEST_RESISTANCE_LEVEL is reported. This means a variant can have both a sensitive and a resistance level.
2. **Conflicting data within Level 1/2/R1:** the SOP notes these are "categorized by their inclusion in either the FDA or NCCN guidelines, and therefore conflicting data is not relevant" — i.e. the level is fixed by guideline inclusion, not by literature vote.

---

## Test Datasets

### Dataset: OncoKB level-ordering examples (oncokb-annotator README)

**Source:** oncokb-annotator README, HIGHEST_LEVEL / HIGHEST_SENSITIVE_LEVEL / HIGHEST_RESISTANCE_LEVEL column definitions.

| Parameter | Value |
|-----------|-------|
| Combined order (highest→lowest) | R1 > 1 > 2 > 3A > 3B > 4 > R2 |
| Sensitive order (highest→lowest) | 1 > 2 > 3A > 3B > 4 |
| Resistance order (highest→lowest) | R1 > R2 |
| Variant with {2, 3A} sensitive associations → highest sensitive | Level 2 |
| Variant with {3A, 3B, 4} sensitive associations → highest sensitive | Level 3A |
| Variant with {1, R1} associations → highest combined | R1 (R1 > 1) |
| Variant with {1} only → highest combined / highest sensitive | Level 1 |
| Variant with {4, R2} → highest combined | Level 4 (4 > R2) |
| Variant with no associations → highest | none (not actionable) |

### Dataset: Level definitions (OncoKB Levels-of-Evidence PDF)

**Source:** OncoKB Therapeutic Levels of Evidence PDF (V2).

| Parameter | Value |
|-----------|-------|
| Level 1 category | Standard Care (sensitivity), FDA-recognized, in-indication |
| Level 2 category | Standard Care (sensitivity), NCCN/guideline, in-indication |
| Level 3A category | Investigational (sensitivity), compelling clinical evidence, in-indication |
| Level 3B category | Investigational (sensitivity), standard/investigational in another indication |
| Level 4 category | Hypothetical (sensitivity), compelling biological evidence |
| Level R1 category | Standard Care Resistance, FDA-approved drug, in-indication |
| Level R2 category | Investigational Resistance, compelling clinical evidence |

---

## Assumptions

1. **ASSUMPTION: VUS / no-association → NotActionable** — The OncoKB sources define levels only for variants that carry a therapeutic implication; they do not define a named level for a variant with zero leveled drug associations. We model this as a distinct `NotActionable` outcome (no level). Justification: the annotator leaves HIGHEST_LEVEL empty for such variants, so "no level" is the documented observable behavior, but the explicit name `NotActionable` is ours.
2. **ASSUMPTION: caller supplies the knowledgebase** — Per the unit scope, drug–gene–level associations are caller-supplied evidence inputs; the library does not embed or reproduce the OncoKB curated database (3,000+ alterations across 418 genes). The classifier ranks caller-supplied levels; it does not look up biomarkers. This is a framework boundary, not an algorithm parameter.

---

## Recommendations for Test Coverage

1. **MUST Test:** Each of Levels 1, 2, 3A, 3B, 4, R1, R2 is parsed/recognized and ranks in the exact order R1 > 1 > 2 > 3A > 3B > 4 > R2 — Evidence: oncokb-annotator README HIGHEST_LEVEL order.
2. **MUST Test:** Highest *sensitive* level of {1,2,3A,3B,4} subsets equals the max under 1 > 2 > 3A > 3B > 4 (e.g. {2,3A}→2; {3A,3B,4}→3A) — Evidence: README HIGHEST_SENSITIVE_LEVEL order.
3. **MUST Test:** Highest *resistance* level of {R1,R2} equals R1 — Evidence: README HIGHEST_RESISTANCE_LEVEL order.
4. **MUST Test:** A variant with both sensitive and resistance associations reports the correct max on each axis and the correct combined max (e.g. {1, R1}→combined R1, sensitive 1, resistance R1) — Evidence: README orders + SOP separate-axis statement.
5. **MUST Test:** A variant with no associations is NotActionable (no level) — Evidence: ASSUMPTION (annotator empty HIGHEST_LEVEL).
6. **SHOULD Test:** Null inputs throw `ArgumentNullException`; the per-variant outputs preserve input order — Rationale: library convention (mirrors `AnnotateCancerVariants`).
7. **COULD Test:** Standard-care vs investigational grouping (1,2,R1 standard; 3A,3B,4,R2 investigational) — Rationale: SOP grouping statement.

---

## References

1. Chakravarty D, Gao J, Phillips SM, et al. (2017). OncoKB: A Precision Oncology Knowledge Base. JCO Precision Oncology 2017:1-16. https://doi.org/10.1200/PO.17.00011
2. OncoKB. Therapeutic Levels of Evidence (V2). Memorial Sloan Kettering Cancer Center. https://www.oncokb.org/content/files/levelOfEvidence/V2/LevelsOfEvidence.pdf (accessed 2026-06-15)
3. OncoKB. Curation Standard Operating Procedure v3. https://sop.oncokb.org/static/sop/OncoKB_Curation_Standard_Operating_Procedure_v3.pdf (accessed 2026-06-15)
4. oncokb-annotator. README (annotation columns: HIGHEST_LEVEL / HIGHEST_SENSITIVE_LEVEL / HIGHEST_RESISTANCE_LEVEL). GitHub. https://github.com/oncokb/oncokb-annotator (accessed 2026-06-15)

---

## Change History

- **2026-06-15**: Initial documentation.
