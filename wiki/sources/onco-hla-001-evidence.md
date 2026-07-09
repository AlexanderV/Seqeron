---
type: source
title: "Evidence: ONCO-HLA-001 (HLA allele nomenclature parsing/validation + allele-specific HLA LOH — LOHHLA)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-HLA-001-Evidence.md
sources:
  - docs/Evidence/ONCO-HLA-001-Evidence.md
source_commit: e7ed7ce32c4de4c278b31b2d397903f64fc2c8c8
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-HLA-001

The validation-evidence artifact for test unit **ONCO-HLA-001** — **HLA allele nomenclature
parsing/validation + allele-specific HLA loss-of-heterozygosity (LOHHLA) classification**. The
**eighteenth ingested unit of the Oncology family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is synthesized in its
own concept, [[hla-nomenclature-and-allele-specific-loh]]; [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (mutually consistent, two disjoint sub-topics, no contradictions):**
  - **WHO HLA Nomenclature — "Naming Alleles" (IPD-IMGT/HLA, hla.alleles.org)** (rank 2, official
    WHO Nomenclature Committee standard) — allele name `HLA-[Gene]*[F1]:[F2][:F3][:F4]`; Field 1 =
    type/allele group, Field 2 = specific HLA protein, Field 3 = synonymous coding difference, Field 4 =
    non-coding (intron / UTR) difference; **two-field ("four-digit") minimum**; optional expression
    suffix `N`/`L`/`S`/`C`/`A`/`Q` (Null / Low / Secreted / Cytoplasm / Aberrant / Questionable).
  - **Marsh et al. (2010)** "Nomenclature for factors of the HLA system, 2010", Tissue Antigens
    75(4):291–455 (rank 1, originating reference for the colon-delimited convention) — up to four
    colon-separated fields; two-field minimum (consistent with hla.alleles.org).
  - **McGranahan et al. (2017)** "Allele-Specific HLA Loss and Immune Escape in Lung Cancer Evolution",
    Cell 171(6):1259–1271, PMC5720478 (rank 1, introduces **LOHHLA**) — per-allele copy number from
    tumor-vs-germline **logR** + **BAF** at allele-distinguishing polymorphic sites; **loss threshold
    copy number < 0.5** (verbatim); **allelic-imbalance guard: paired Student's t-test p < 0.01**
    (verbatim, "to avoid over-calling LOH").
  - **LOHHLA reference implementation (`mskcc/lohhla`, `LOHHLAscript.R`)** (rank 3) — confirms the
    paired t-test (`t.test(...,paired=TRUE)`, `PVal <- PairedTtest$p.value`) and per-allele copy-number
    variables `HLA_type{1,2}copyNum_with{,out}BAF` (the classifier's input).

- **Documented corner cases / failure modes:** nomenclature — single-field name (`HLA-A*02`) invalid,
  trailing suffix must be from N/L/S/C/A/Q, 2–4 fields valid (>4 invalid); LOHHLA — **no-LOH despite
  low copy** when p ≥ 0.01 (over-calling guard), both alleles CN ≥ 0.5 → heterozygous-retained,
  homozygous locus cannot be assessed (no polymorphic sites distinguish the homologs).

- **Datasets (deterministic worked oracles):**
  - **Nomenclature:** `HLA-A*02:01` (valid, 2 fields) · `HLA-B*07:02:01` (3) · `HLA-C*07:02:01:03` (4) ·
    `HLA-A*24:02:01:02L` (4 + `L`); invalid `HLA-A*02` (one field) · `A*02:01` (no prefix) ·
    `HLA-A*02:01:01:01:01` (five fields) · `HLA-A*02:01X` (`X` not a valid suffix).
  - **HLA LOH (LOHHLA thresholds):** (1.8, 0.30, p=0.001) → LOH, lost allele 2 · (0.10, 1.50, 0.0005) →
    LOH, allele 1 · (1.10, 0.90, 0.30) → no (both retained) · (1.60, 0.40, 0.05) → no (p ≥ 0.01 guard) ·
    (1.50, 0.50, 0.001) → no (0.5 not < 0.5) · (1.70, 0.40, 0.01) → no (0.01 not < 0.01). Both
    thresholds strict `<`.

- **Coverage recommendations:** MUST parse each valid name into gene + field tuple + suffix; MUST reject
  malformed names (missing prefix / single field / five fields / invalid suffix / non-numeric field);
  MUST call LOH iff one allele CN < 0.5 AND imbalance p < 0.01, verifying the lost allele and the two
  boundary cases (CN=0.5, p=0.01); SHOULD verify both-retained and non-significant-imbalance → no LOH;
  COULD verify both-alleles-< 0.5 → `HomozygousLoss` label.

## Deviations and assumptions

- **ASSUMPTION — lost-allele tie-break (both alleles < 0.5):** McGranahan 2017 defines a *lost* allele
  by CN < 0.5 but does not specify behaviour when **both** alleles are < 0.5 with significant imbalance
  (biologically a homozygous deletion, not allele-specific LOH). The unit calls LOH only when exactly
  one allele is < 0.5; both < 0.5 → distinct **`HomozygousLoss`** label. Only the label is affected —
  the two source-exact thresholds (0.5, 0.01) are unchanged.

No source contradictions — the WHO/IPD-IMGT/HLA nomenclature standard (naming rules) and the LOHHLA
paper + reference implementation (LOH thresholds) cover disjoint parts of the unit and are mutually
consistent.
