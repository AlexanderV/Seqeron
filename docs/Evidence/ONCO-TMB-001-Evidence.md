# Evidence Artifact: ONCO-TMB-001

**Test Unit ID:** ONCO-TMB-001
**Algorithm:** Tumor Mutational Burden (TMB) — mutations/Mb and TMB-high classification
**Date Collected:** 2026-06-14

---

## Online Sources

### Chalmers et al. (2017), Genome Medicine — "Analysis of 100,000 human cancer genomes reveals the landscape of tumor mutational burden"

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC5395719 (PMC full text; DOI 10.1186/s13073-017-0424-2)
**Retrieved by:** WebSearch "Chalmers 2017 Genome Medicine tumor mutational burden ... formula" → WebFetch of the BioMedCentral article (301 → Springer, paywalled) → WebFetch of the open PMC full text `https://pmc.ncbi.nlm.nih.gov/articles/PMC5395719`.
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed primary paper)

**Key Extracted Points:**

1. **TMB definition / formula (Methods):** "TMB was defined as the number of somatic, coding, base substitution, and indel mutations per megabase of genome examined." (verbatim from the retrieved Methods text). Thus **TMB = (counted mutations) / (sequenced coding region in Mb)**, units mut/Mb.
2. **Coding region size of the assay:** the assay targets "315 genes (1.1 Mb of coding genome)" — i.e. the denominator for the 315-gene FoundationOne panel is **1.1 Mb**.
3. **Mutation counting / synonymous:** "All base substitutions and indels in the coding region of targeted genes, including synonymous alterations, are initially counted before filtering." Synonymous alterations are counted "in order to reduce sampling noise."
4. **Germline / driver filtering (not counted):** "Alterations listed as known somatic alterations in COSMIC and truncations in tumor suppressor genes were not counted"; "Alterations predicted to be germline by the somatic-germline-zygosity algorithm were not counted"; "Known germline alterations in dbSNP were not counted. Germline alterations occurring with two or more counts in the ExAC database were not counted."
5. **Small-panel variance:** percentage deviation of panel TMB vs whole-exome TMB increases "especially with less than 0.5 Mb sequenced" — below ~0.5 Mb the estimate is unstable.

### FDA Approval Summary: Pembrolizumab for TMB-High Solid Tumors (Marcus et al., Clin Cancer Res 2021; PMC8416776)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC8416776/
**Retrieved by:** WebSearch "FDA pembrolizumab tumor mutational burden 10 mutations per megabase FoundationOne CDx" → WebFetch of the PMC article.
**Accessed:** 2026-06-14
**Authority rank:** 2 (official FDA regulatory summary, peer-reviewed journal)

**Key Extracted Points:**

1. **TMB-high cutoff:** the approval defines TMB-H as **"≥10 mutations/megabase (mut/Mb)"** (verbatim, abstract and throughout).
2. **Approval date:** "June 16, 2020."
3. **Companion diagnostic:** the FDA-approved test for TMB determination is **FoundationOne CDx (F1CDx)**.
4. **TMB meaning:** "TMB reflects the overall somatic genomic burden of mutations within a given tumor"; in the WES analyses "TMB was assessed as the number of nonsynonymous single nucleotide variants and indels found in protein coding regions."

### Fancello / Sholl review — "Tumor Mutational Burden (TMB) as a Predictive Biomarker in Solid Tumors" (Friends of Cancer Research harmonization context; PMC7710563)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7710563/
**Retrieved by:** WebSearch "Friends of Cancer Research TMB harmonization project Merino 2020 mutations per megabase definition coding" → WebFetch of the PMC article.
**Accessed:** 2026-06-14
**Authority rank:** 1–2 (peer-reviewed review summarizing the Friends of Cancer Research TMB Harmonization Project, Merino et al. 2020 J Immunother Cancer 8:e000147)

**Key Extracted Points:**

1. **Harmonized definition:** "TMB, defined as the number of somatic mutations per megabase of interrogated genomic sequence." The Harmonization Consortium recommends **reporting TMB in mutations/megabase (mut/Mb)** to keep values comparable across studies.
2. **Coding / non-synonymous:** "MSK-IMPACT and FoundationOne CDx® detect somatic coding mutations (non-synonymous) ... inclusive of frameshift, point mutations, and small insertions and deletions (indels)." (the clinically reported TMB is over coding mutations).
3. **FDA cutoff confirmation:** "The U.S. Food and Drug Administration (FDA) approved pembrolizumab monotherapy for the subgroup of solid tumor patients with TMB ≥10 mut/Mb."
4. **FoundationOne CDx panel size:** "FoundationOne CDx® ... cover ... ~0.8 Mb over 324 genes" (the current F1CDx version; the historical 315-gene panel in Chalmers 2017 is 1.1 Mb).

---

## Documented Corner Cases and Failure Modes

### From Chalmers et al. (2017)

1. **Small denominator instability:** below ~0.5 Mb of sequenced coding region, the percentage deviation of panel-based TMB from whole-exome TMB rises sharply — panel TMB estimates become unreliable.
2. **Filtering before counting:** known driver/germline alterations are removed before the count, so a raw variant list overstates TMB if not filtered (this library counts the caller-supplied somatic-mutation list, leaving filtering to the upstream somatic caller, per ONCO-SOMATIC-001).

### From the FDA Approval Summary (PMC8416776)

1. **Threshold is inclusive (≥):** TMB-H is `TMB ≥ 10`, so exactly 10 mut/Mb is TMB-high (boundary is inclusive).

---

## Test Datasets

### Dataset: FoundationOne 315-gene panel worked example (Chalmers 2017)

**Source:** Chalmers et al. (2017), Genome Medicine — "315 genes (1.1 Mb of coding genome)"; TMB = mutations / Mb.

| Parameter | Value |
|-----------|-------|
| Panel coding size | 1.1 Mb |
| Counted somatic mutations | 11 |
| Expected TMB | 11 / 1.1 = 10.0 mut/Mb |
| Classification (cutoff ≥10) | TMB-High |

### Dataset: Whole-exome example (Friends of Cancer Research / common reporting)

**Source:** Definition TMB = mutations / Mb (Merino/FoCR; Chalmers). Exome coding ≈ 30–40 Mb is the WES denominator; the algorithm takes the denominator as a parameter.

| Parameter | Value |
|-----------|-------|
| Counted mutations | 300 |
| Sequenced region | 30 Mb |
| Expected TMB | 300 / 30 = 10.0 mut/Mb |
| Classification (cutoff ≥10) | TMB-High |

### Dataset: FDA cutoff boundary cases

**Source:** FDA Approval Summary (PMC8416776) — TMB-H ≡ TMB ≥ 10 mut/Mb.

| Counted mutations | Region (Mb) | TMB (mut/Mb) | Classification |
|-------------------|-------------|--------------|----------------|
| 99 | 10 | 9.9 | Not high (< 10) |
| 100 | 10 | 10.0 | High (= 10, inclusive) |
| 150 | 10 | 15.0 | High |
| 0 | 10 | 0.0 | Not high |

---

## Assumptions

1. **ASSUMPTION: Two-tier (High vs Not-High) classification using the single FDA ≥10 cutoff.** The only TMB threshold retrieved from authoritative sources is the FDA/F1CDx TMB-high cutoff of **≥10 mut/Mb** (PMC8416776). The Registry by-area note lists "Low (<6/Mb), Intermediate (6–20/Mb), High (>20/Mb)", but **no authoritative source was retrieved that defines those 6 and 20 boundaries**; they are tumor-type-specific research cut-points, not a harmonized standard. Per the evidence-first / no-fabrication policy, `ClassifyTMB` implements the source-backed **TMB-High = TMB ≥ 10**, TMB-Low = TMB < 10. The unsupported 6/20 boundaries are NOT implemented (would be fabricated constants). Conflict noted in TestSpec §7 and the checklist.

---

## Recommendations for Test Coverage

1. **MUST Test:** `CalculateTMB` returns count / regionMb for worked examples (11 mut / 1.1 Mb = 10.0; 300 / 30 = 10.0; 150 / 10 = 15.0) — Evidence: Chalmers 2017 formula "number of ... mutations per megabase of genome examined".
2. **MUST Test:** `ClassifyTMB` boundary — 9.9 → Low, 10.0 → High (inclusive), 15.0 → High, 0 → Low — Evidence: FDA Approval Summary "TMB ≥10 mut/Mb".
3. **MUST Test:** `CalculateTMB` with `targetRegionMb = 0` throws (division by zero) — Evidence: Registry edge case + formula has Mb in denominator.
4. **SHOULD Test:** small panel (< 0.5 Mb) still computes the ratio (no exception) but is flagged/known-unstable — Rationale: Chalmers 2017 variance note (the value is still mathematically defined; instability is documentation, not an error).
5. **SHOULD Test:** negative mutation count / negative region rejected — Rationale: counts and region size are non-negative by definition.
6. **COULD Test:** invariant — TMB is monotone non-decreasing in mutation count for fixed region, and non-increasing in region size for fixed count — Rationale: division property of the formula.

---

## References

1. Chalmers ZR, Connelly CF, Fabrizio D, et al. (2017). Analysis of 100,000 human cancer genomes reveals the landscape of tumor mutational burden. Genome Medicine 9:34. https://doi.org/10.1186/s13073-017-0424-2 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC5395719)
2. Marcus L, Fashoyin-Aje LA, Donoghue M, et al. (2021). FDA Approval Summary: Pembrolizumab for the Treatment of Tumor Mutational Burden–High Solid Tumors. Clinical Cancer Research 27(17):4685–4689. https://doi.org/10.1158/1078-0432.CCR-21-0327 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC8416776/)
3. Fancello L, Gandini S, Pelicci PG, Mazzarella L (2019); and Friends of Cancer Research TMB Harmonization Project (Merino DM, McShane LM, et al. 2020, J Immunother Cancer 8:e000147). Tumor Mutational Burden as a Predictive Biomarker in Solid Tumors (review). https://pmc.ncbi.nlm.nih.gov/articles/PMC7710563/

---

## Change History

- **2026-06-14**: Initial documentation (ONCO-TMB-001).
