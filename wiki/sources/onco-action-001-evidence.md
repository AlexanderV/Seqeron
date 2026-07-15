---
type: source
title: "Evidence: ONCO-ACTION-001 (clinical actionability — OncoKB therapeutic levels of evidence)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-ACTION-001-Evidence.md
sources:
  - docs/Evidence/ONCO-ACTION-001-Evidence.md
source_commit: f43bbbd5379ef59f0f01dabd2d4a7378e8cb23da
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ONCO-ACTION-001

The validation-evidence artifact for test unit **ONCO-ACTION-001** — **Clinical Actionability
Assessment** by the **OncoKB Therapeutic Levels of Evidence**. This is the **first ingested unit of the
Oncology family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The distinct level-ranking method is
synthesized in its own concept, [[clinical-actionability-oncokb-levels]]; [[test-unit-registry]] tracks
the unit.

## What this file records

- **Online sources (all mutually consistent, no contradictions):**
  - **Chakravarty D et al. (2017)** OncoKB: A Precision Oncology Knowledge Base, *JCO Precision
    Oncology* 2017:1–16 (rank 1, primary; DOI landing paywalled HTTP 403, so level text taken from the
    two official OncoKB PDFs below which restate the same system) — OncoKB stratifies a somatic
    alteration's **treatment implications by an ordered level of evidence** that it is predictive of
    drug response, weighted by FDA labeling / NCCN guidelines / expert-group recommendations /
    literature. This ordered actionability axis is exactly what ONCO-ACTION-001 classifies.
  - **OncoKB Therapeutic Levels of Evidence PDF (V2)** (rank 3, canonical project doc) — the **verbatim
    seven level definitions**: Level 1 (FDA-recognized biomarker, in-indication, Standard Care), Level 2
    (NCCN/guideline biomarker, in-indication, Standard Care), Level 3A (compelling **clinical** evidence,
    in-indication, Investigational), Level 3B (standard/investigational biomarker in **another**
    indication, Investigational), Level 4 (compelling **biological** evidence, Hypothetical), R1
    (standard-care resistance to FDA-approved drug, in-indication), R2 (compelling clinical evidence of
    resistance, Investigational).
  - **OncoKB Curation SOP v3 PDF** (rank 3, canonical project doc) — the **grouping**: Levels 1/2/R1 are
    the standard (FDA/NCCN) implications; Levels 3A/3B/4/R2 are investigational. Confirms **3A ranks
    above 3B** (system refined to deprioritize standard-care biomarkers used outside the approved
    indication). Provenance: developed in Chakravarty 2017, **consistent with the AMP/ASCO/CAP Joint
    Consensus (Li et al. 2017)**. Levels 1/2/R1 are fixed by guideline inclusion, so conflicting
    literature is not relevant to them.
  - **oncokb-annotator README** (rank 3, reference implementation) — the **three ordering axes** used as
    the primary oracle:
    - **Combined `HIGHEST_LEVEL`:** `R1 > 1 > 2 > 3A > 3B > 4 > R2` (resistance R1 interleaves *above*
      sensitivity Level 1; R2 sits *below* Level 4).
    - **`HIGHEST_SENSITIVE_LEVEL`:** `1 > 2 > 3A > 3B > 4`.
    - **`HIGHEST_RESISTANCE_LEVEL`:** `R1 > R2`.
    - For a variant with several leveled drug associations, the actionable level is the **maximum** under
      the applicable order.

- **Documented corner cases / failure modes:** a variant with **no leveled drug association** leaves the
  annotator's `HIGHEST_LEVEL` **empty** → modeled as a distinct **`NotActionable`** outcome (no level).
  **Resistance and sensitivity are separate axes** — a variant can carry both a sensitive and a
  resistance level; the combined order interleaves them but the two axis-specific highest levels are
  reported independently. Levels 1/2/R1 are **fixed by FDA/NCCN inclusion**, not by a literature vote,
  so "conflicting data" does not move them.

- **Datasets (documented oracles, from the annotator README ordering):**
  - `{2, 3A}` sensitive associations → highest sensitive = **Level 2**.
  - `{3A, 3B, 4}` sensitive → highest sensitive = **Level 3A**.
  - `{1, R1}` → highest **combined = R1** (R1 > 1), highest sensitive = 1, highest resistance = R1.
  - `{1}` only → combined / sensitive = **Level 1**.
  - `{4, R2}` → highest combined = **Level 4** (4 > R2).
  - `{R1, R2}` resistance → highest resistance = **R1**.
  - No associations → highest = **none (NotActionable)**.

## Deviations and assumptions

- **ASSUMPTION 1 — VUS / no-association → `NotActionable`.** OncoKB defines levels only for variants
  carrying a therapeutic implication; it names no level for zero leveled associations. The unit models
  this as a distinct `NotActionable` outcome. Justification: the annotator leaves `HIGHEST_LEVEL` empty
  for such variants ("no level" is the documented observable), but the explicit **name** `NotActionable`
  is the library's, not OncoKB's.
- **ASSUMPTION 2 — caller supplies the knowledgebase.** Per the unit scope, drug–gene–level
  associations are **caller-supplied evidence inputs**; the library does **not** embed or reproduce the
  OncoKB curated database (3,000+ alterations across 418 genes). The classifier **ranks** caller-supplied
  levels; it does **not** look up biomarkers. This is a framework boundary, not an algorithm parameter.
- **Coverage recommendations:** MUST-test each of the seven levels ranks in the exact combined order;
  each axis-specific highest level; a mixed sensitive+resistance variant reports the correct max per
  axis and combined; no-association → `NotActionable`. SHOULD-test null → `ArgumentNullException` and
  input-order-preserving per-variant outputs (mirrors `AnnotateCancerVariants`). COULD-test the
  standard-care (1/2/R1) vs investigational (3A/3B/4/R2) grouping.

No source contradictions — Chakravarty 2017, the OncoKB Levels-of-Evidence PDF, the OncoKB SOP v3, and
the oncokb-annotator README are mutually consistent (the SOP is explicitly consistent with the
AMP/ASCO/CAP Li et al. 2017 consensus).
