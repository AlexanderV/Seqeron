---
type: concept
title: "Clonal hematopoiesis (CHIP) filtering for cfDNA liquid biopsy"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-CHIP-001-Evidence.md
source_commit: 90f75a142c015ef57f04ebf747b01f8b855634db
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-chip-001-evidence
      evidence: "Test Unit ID: ONCO-CHIP-001 ... Algorithm: Clonal Hematopoiesis (CHIP) Filtering for cfDNA Liquid Biopsy"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:sequencing-artifact-detection
      source: onco-chip-001-evidence
      evidence: "Razavi 2019: CH is the dominant cfDNA confounder (81.6% controls / 53.2% cancer patients). CHIP filtering removes CH-derived false positives before clinical interpretation, the biological-origin sibling of the technical-artifact FilterArtifacts QC filter."
      confidence: high
      status: current
---

# Clonal hematopoiesis (CHIP) filtering for cfDNA liquid biopsy

The **pre-interpretation biological-origin filter** of the Oncology family: in a **liquid biopsy**
(cell-free DNA / cfDNA from plasma), a large fraction of somatic-looking variants are **not tumor** —
they come from **clonal hematopoiesis (CHIP)**, the age-related expansion of blood-cell clones
carrying leukemia-driver mutations. Razavi 2019 measured this as the **dominant false-positive
class**: **81.6%** of cfDNA mutations in healthy controls and **53.2%** in cancer patients have CH
features. This unit identifies and removes those CH-derived calls before variants reach the
clinical-interpretation ONCO layers. Validated under test unit **ONCO-CHIP-001**; the
literature-traced record is [[onco-chip-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the evidence-artifact pattern. Research-grade
([[scientific-rigor|research-grade]]), **not for clinical or diagnostic use**.

## 1. The CHIP definition (Steensma 2015 / Genovese 2014)

**CHIP** = **a somatic mutation in a driver gene recurrently mutated in hematologic malignancies, at
VAF ≥ 2%, in a person with no diagnosed hematologic malignancy or MDS** (Steensma et al. 2015,
*Blood*, which coined the term). The three parts:

- **VAF ≥ 0.02 (inclusive).** Below 2% is not CHIP — both below the formal threshold and the reliable
  detection limit ("with deep enough sequencing a mutation can be found in every individual").
- **Driver gene.** The canonical set (Steensma Fig 2A + Genovese recurrent genes) defaults to
  **{DNMT3A, TET2, ASXL1, TP53, JAK2, SF3B1, SRSF2, PPM1D}**; **top-3 by prevalence DNMT3A / TET2 /
  ASXL1** (Genovese: DNMT3A 190 mutations, ASXL1 35, TET2 31). The panel is a *labeled* canonical
  default, **caller-overridable**.
- **Absence of malignancy.** If diagnostic criteria for a hematologic malignancy are met it is not
  CHIP — out of scope here (the assay reads non-diagnostic plasma).

CH variants are **sub-clonal (VAF < 0.5** in blood, Genovese) — a single expanded clone, not germline.

## 2. Three methods

### `IdentifyCHIPVariants` — candidate flagging (gene + VAF)

Flags a variant as CHIP iff it is **in a CHIP gene AND VAF ≥ 0.02**. Gene comparison is
**case-insensitive** (HGNC symbols are upper-case). This is a *candidate* flag only — Arango-Argoty
2025 cautions that "the exact relationship between VAF and variant origin remains unclear", so gene +
VAF alone does not prove blood origin.

### `FilterCHIP` — removal, with two rules

Removes CH-derived variants from a cfDNA call set:

- **Rule (a) — matched-WBC subtraction (the definitive origin test, Razavi 2019).** A cfDNA variant
  **also present in the patient's matched white-blood-cell (WBC) DNA is WBC/CH-derived, not tumor**,
  and is removed — **even a non-CHIP-gene variant**. Matched cfDNA + WBC sequencing is the
  gold-standard design; "present in matched WBC" is decided by an `IsVariantDetected`-style ≥1-alt-read
  test at the same locus (repo MRD convention, Wan 2020; alt-read cutoff configurable).
- **Rule (b) — labelled gene+VAF heuristic fallback.** When no matched WBC is available, a
  CH-driver-gene variant at VAF ≥ τ is removed anyway. This is **deliberately conservative** — it
  **over-removes** relative to the strict matched-WBC definition. Callers who want the strict rule
  pass an **empty/custom `chipGenes` panel** so only matched-WBC subtraction applies.

A cfDNA variant **absent** from matched WBC is **retained** as a candidate tumor variant even if it
lies in a CHIP gene — presence in WBC, not gene identity, is the definitive origin test.

### `CallVariantOrigin` — strict origin call (Bolton 2020)

The strict matched-WBC origin rule from the 24,439-patient MSK-IMPACT study. Returns **Chip** iff all
three hold, else **Tumor**:

```
WBC VAF ≥ 0.02   AND   WBC supporting reads ≥ 10   AND   WBC VAF ≥ φ × tumour VAF
```

- **φ = 2.0** by default; **φ = 1.5 when the tumour biopsy site is a lymph node** (leukocyte
  infiltration inflates tumour-fraction there). The fold ratio is a **sourced parameter** chosen via
  Bolton's leukocyte-contamination simulations — exposed as a caller knob, not invented.
- **Absent from WBC ⇒ Tumor** (no blood evidence), even for a CH driver gene.

## 3. Worked oracles

**Gene+VAF classification** (gene / VAF / in-WBC → IdentifyCHIP / FilterCHIP kept?):

| Gene | VAF | in WBC | IdentifyCHIP | FilterCHIP |
|------|-----|--------|--------------|------------|
| DNMT3A | 0.05 | — | CHIP | removed |
| DNMT3A | 0.01 | — | not CHIP (VAF<0.02) | kept |
| EGFR | 0.30 | — | not CHIP (not a CHIP gene) | kept |
| EGFR | 0.30 | yes | n/a (gene rule) | removed (WBC-matched) |
| TP53 | 0.40 | no | CHIP candidate | removed (heuristic rule b) |

**Strict origin** (Bolton, φ default 2.0 / 1.5 lymph node): tumour 0.10 / WBC 0.30 / 40 reads → **Chip**;
tumour 0.40 / WBC absent → **Tumor**; WBC 9 reads → **Tumor** (9 < 10); tumour 0.01 / WBC 0.02 / 10 reads
→ **Chip** (all boundaries inclusive: 0.02=2%, 0.02=2×0.01, 10=10); WBC VAF 0.015 → **Tumor** (< 2% floor);
tumour 0.30 / WBC 0.40 → **Tumor** (0.40 < 2×0.30=0.60); tumour 0.25 / WBC 0.40 at φ=1.5 → **Chip**
(0.40 ≥ 0.375) yet **Tumor** under the default 2.0.

## 4. Relationship to the rest of the Oncology family

CHIP filtering is a **QC / origin filter that runs before clinical interpretation**, the
**biological-origin sibling** of the technical-artifact filter [[sequencing-artifact-detection]]
(`FilterArtifacts` removes OxoG / FFPE-deamination / strand-bias *sequencing* artifacts; this removes
*blood-clone* biology). Cleaned tumor variants then flow into the clinical-significance units
[[clinical-actionability-oncokb-levels]] and [[cancer-variant-tier-classification-amp-asco-cap]], and
the clonal-structure layers [[allele-specific-copy-number-ascat]] /
[[cancer-cell-fraction-clonal-clustering]] — a CH variant left unfiltered would corrupt every one of
those downstream calls.

## 5. Corner cases and scope

- **Sub-2% driver mutation** — not CHIP (below threshold + detection limit).
- **Driver mutation with a diagnosed malignancy** — not CHIP (out of scope; assays non-diagnostic plasma).
- **No matched WBC** — only the conservative gene+VAF heuristic (rule b) applies; it over-removes.
- **Variant absent from matched WBC** — retained as candidate tumor even in a CHIP gene.
- **Null / empty inputs** — documented failure modes; kept-variant order preservation is a COULD contract.

Sources are mutually consistent: the VAF-2% + driver-gene definition (Steensma 2015 / Genovese 2014),
the matched-WBC origin gold standard (Razavi 2019 / Arango-Argoty 2025), and the strict fold-ratio
origin rule (Bolton 2020) each cover a disjoint stage and reinforce one another. **Two flagged
assumptions**, both source-consistent: the canonical default gene set and the ≥1-alt-read matched-WBC
presence test. **Not for clinical or diagnostic use.**
</content>
