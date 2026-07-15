---
type: source
title: "Evidence: ONCO-ANNOT-001 (cancer variant tier classification — AMP/ASCO/CAP 2017 four-tier)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-ANNOT-001-Evidence.md
sources:
  - docs/Evidence/ONCO-ANNOT-001-Evidence.md
source_commit: c38dde50786144e5976f703468375f68abe6fd1b
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ONCO-ANNOT-001

The validation-evidence artifact for test unit **ONCO-ANNOT-001** — **Cancer-Specific Variant
Annotation** by the **AMP/ASCO/CAP 2017 four-tier clinical-significance classification**
(`AnnotateCancerVariants` + `GetCOSMICAnnotation`). This is the **second ingested unit of the Oncology
family** and one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence
artifact]] pattern. The distinct four-tier decision rule is synthesized in its own concept,
[[cancer-variant-tier-classification-amp-asco-cap]]; the sibling therapeutic-actionability ranking is
[[clinical-actionability-oncokb-levels]]; [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (both mutually consistent, no contradictions):**
  - **Li MM et al. (2017)** "Standards and Guidelines for the Interpretation and Reporting of Sequence
    Variants in Cancer" — the AMP/ASCO/CAP Joint Consensus, *J Mol Diagn* 19(1):4–23 (rank 1, primary;
    Elsevier full text paywalled, full guideline read this session from an authoritative open MCW-hosted
    PDF of the same article). The **four-tier system** (verbatim abstract): Tier I strong clinical
    significance, Tier II potential clinical significance, Tier III unknown clinical significance, Tier IV
    benign/likely benign. The **tier ↔ evidence-level mapping** (Figure 2): Tier I = Level **A/B**, Tier
    II = Level **C/D**, Tier III/IV = no level. The **four evidence levels** (Table 3, across Therapeutic
    / Diagnosis / Prognosis categories): A = FDA-approved therapy / guidelines; B = well-powered studies
    + expert consensus; C = off-label (different tumor type) / trial inclusion / multiple small studies;
    D = preclinical / few case reports. The **1% (0.01) MAF primary benign cutoff** (verbatim): "no
    standardized cutoff ... the work group recommends using 1% (0.01) as a primary cutoff". **Population
    frequency is the Tier III ↔ IV discriminator**: rare + cancer association → Tier III; common (MAF ≥
    1%) or no cancer association → Tier IV (Figure 2 boxes + Tables 6/7).
  - **Tate JG et al. (2019)** "COSMIC: the Catalogue Of Somatic Mutations In Cancer" — *Nucleic Acids
    Res* 47(D1):D941 (rank 1 paper / 5 database). COSMIC is an **external curated database** (v86: ~6M
    coding mutations, 1.4M tumour samples, 26,000+ publications) that **cannot be reproduced/hardcoded**;
    Li 2017 Tables 4–6 list "Somatic database: COSMIC, My Cancer Genome, TCGA" as an evidence source
    (Tier I "Most likely present", Tier II "Likely present"). COSMIC presence **supports** but does not
    by itself **determine** the tier.

- **Documented corner cases / failure modes:**
  - Common polymorphism (MAF ≥ 1%) with no clinical evidence → **Tier IV** (Table 7).
  - Clinically significant variant that also appears in population databases → assigned by **evidence
    level** (Tier I for A/B), **not downgraded** by frequency (well-studied germline counterparts like
    `TP53`, `PTEN` may appear in databases yet stay significant).
  - Rare variant + cancer association + no evidence level → **Tier III** (distinguished from Tier IV by
    the cancer association and absence of significant frequency).
  - No cancer association and not common → **Tier IV** (Figure 2 Tier IV box).
  - Variant **absent** from the supplied COSMIC catalog → lookup returns **null (not found)**, never a
    fabricated annotation.

- **Datasets (documented oracles, Figure 2 / Tables 3–7):**
  - Tier decision matrix: A → Tier I; B → Tier I; C → Tier II; D → Tier II; none + MAF ≥ 0.01 → Tier IV;
    none + MAF < 0.01 + no assoc → Tier IV; none + MAF < 0.01 + assoc → Tier III.
  - Canonical variants: `BRAF p.V600E` (Level A, MAF 0.0, assoc) → **Tier I**; rare VUS (none, 0.0001,
    assoc) → **Tier III**; common SNP (none, 0.25, no assoc) → **Tier IV**.
  - COSMIC lookup: `{(BRAF, p.V600E) → "COSV56056643"}`; key `(BRAF, p.V600E)` → `"COSV56056643"`; key
    `(TP53, p.R175H)` → **null** (not in catalog).
  - Boundary: MAF exactly **0.01** → Tier IV (cutoff is "≥ 1%"); **0.0099** + assoc → Tier III.

## Deviations and assumptions

- **ASSUMPTION 1 — caller-supplied evidence inputs.** The AMP/ASCO/CAP guideline classifies from external
  curated knowledge (professional guidelines, population/somatic databases, literature). The library does
  **not** reproduce those resources; the evidence level, population MAF, and cancer-association flag are
  supplied by the caller. This is an input-shape decision, not a correctness-affecting one — the Figure 2
  tiering rule is applied verbatim to whatever evidence is supplied.
- **ASSUMPTION 2 — Tier III vs Tier IV discriminator.** When no evidence level (A–D) is present, the
  implementation uses **MAF ≥ 1% OR absence of a cancer association ⇒ Tier IV, otherwise Tier III** — a
  direct reading of the Figure 2 boxes and Tables 6/7 (the "≥" glyph lost in PDF extraction is fixed by
  the guideline's 1% primary cutoff for benign), not invented.
- **Coverage recommendations:** MUST-test Level A/B → Tier I; Level C/D → Tier II; no level + MAF ≥ 1% →
  Tier IV; no level + MAF < 1% + no assoc → Tier IV; no level + MAF < 1% + assoc → Tier III; Level A
  overrides high MAF (still Tier I); `AnnotateCancerVariants` one-annotation-per-variant in input order
  for a mixed batch; `GetCOSMICAnnotation` hit-value / miss-null; MAF boundary 0.01 → Tier IV, 0.0099 +
  assoc → Tier III. SHOULD-test invalid MAF (negative / > 1 / NaN) → `ArgumentOutOfRangeException`, null →
  `ArgumentNullException`, empty collection → empty list.

No source contradictions — Li et al. 2017 and Tate et al. 2019 are mutually consistent; the AMP/ASCO/CAP
consensus is the same one the OncoKB SOP (sibling [[clinical-actionability-oncokb-levels]]) declares
itself consistent with.
