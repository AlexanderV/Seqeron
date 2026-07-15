---
type: source
title: "Evidence: ONCO-CNA-003 (homozygous / deep deletion detection — total-CN-0 call + tumour-suppressor mapping)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-CNA-003-Evidence.md
sources:
  - docs/Evidence/ONCO-CNA-003-Evidence.md
source_commit: 819918712f8e6a3fddb0f4a534fb6f69bc24cf5b
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ONCO-CNA-003

The validation-evidence artifact for test unit **ONCO-CNA-003** — **homozygous (deep) deletion
detection**: identify segments whose classified integer copy number is **exactly 0** (a homozygous /
deep deletion), then **map** each such segment's arm to a panel of known **tumour-suppressor** genes.
The **deletion counterpart of ONCO-CNA-002** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is synthesized in its
own concept, [[homozygous-deletion-detection]]; [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (five, mutually consistent, no contradictions):**
  - **cBioPortal — Discrete Copy Number file format + FAQ** (rank 5, de-facto discrete-CNA convention) —
    the verbatim discrete scale: `−2` = **Deep Deletion**, "a deep loss, possibly a **homozygous
    deletion**"; `−1` = Shallow / single-copy (heterozygous) loss; `0` diploid; `1` low-level gain;
    `2` high-level amplification. The value **−2 ↔ homozygous deletion** — the deepest discrete loss.
  - **Cheng et al. 2017, Nature Communications 8:1221** (rank 1, peer-reviewed) — **homozygous deletion
    = total copy number 0**: "regions having zero copies of **both alleles** in the tumour cells." A
    hemizygous (single-copy / heterozygous) deletion loses only **one** allele; a homozygous deletion
    requires "two independent hits." Homozygous deletions are rare and recurrently target **tumour
    suppressors** (RB1, CDKN2A, PTEN historically delineated this way).
  - **CNVkit `cnvlib/call.py` `absolute_threshold`** (rank 3, reused from ONCO-CNA-001) — integer copy
    number **0 ⇒ `DeepDeletion`** state (log2 ≤ −1.1 by the default thresholds). A homozygous deletion
    is therefore **exactly** the existing CN-0 / `DeepDeletion` classification — **no new numeric
    threshold is invented**; it reuses [[copy-number-alteration-classification]].
  - **NCBI Gene — tumour-suppressor cytogenetic locations** (rank 5, curated DB) — the registry panel
    and arms: **TP53 17p13.1 (17p) · RB1 13q14.2 (13q) · CDKN2A 9p21.3 (9p) · PTEN 10q23.31 (10q) ·
    BRCA1 17q21.31 (17q) · BRCA2 13q13.1 (13q)**. The chromosome+arm prefix is what a homozygous-
    deletion segment's arm label is matched against; **13q carries both RB1 and BRCA2**.

- **Documented corner cases / failure modes:** a single-copy loss (−1, heterozygous, CN ≥ 1) must NOT
  be reported as a homozygous deletion — only total CN exactly 0 qualifies; discrete calls are putative
  and purity/ploidy differences cause false +/− (interpretation caveat only, does not move the CN-0
  definition); one remaining copy = hemizygous, not homozygous.

- **Datasets (worked, deterministic):**
  - Discrete CNA scale (cBioPortal): only **−2** → homozygous deletion; −1 / 0 / 1 / 2 → not.
  - Integer CN ⇒ state (CNVkit, via ONCO-CNA-001): log2 −2.0 (≤ −1.1) → CN 0 → `DeepDeletion` →
    **homozygous**; log2 −0.5 → CN 1 → `Loss` → not; log2 0.0 → CN 2 → `Neutral` → not.
  - Tumour-suppressor arms (NCBI Gene): 17p→TP53, 13q→RB1 **and** BRCA2, 9p→CDKN2A, 10q→PTEN, 17q→BRCA1.

- **Coverage recommendations (8 items):** MUST — CN-0 segment reported as homozygous deletion; CN-1
  single-copy loss NOT reported; neutral/gain/amp NOT reported; `IdentifyDeletedTumorSuppressors` maps
  17p→TP53, 13q→RB1+BRCA2, 9p→CDKN2A, 10q→PTEN, 17q→BRCA1; custom thresholds move the CN-0 boundary.
  SHOULD — boundary log2 exactly at the deletion cutoff (−1.1) is CN 0 (≤, inclusive) → homozygous;
  order-preserving report-each-once (mirror `DetectFocalAmplifications`). COULD — empty→empty,
  null→ArgumentNullException, invalid arm length / End≤Start→ArgumentException (mirror ONCO-CNA-002).

## Deviations and assumptions

- **ASSUMPTION — homozygous deletion identified at the integer-CN level via the existing CN-0
  (`DeepDeletion`) classification.** cBioPortal (−2 = Deep Deletion) and Cheng et al. (total CN 0)
  converge on the same call the repository already realizes (CNVkit, ONCO-CNA-001). A segment is a
  homozygous deletion **iff its classified integer copy number is 0** — no new threshold.
- **ASSUMPTION — curated tumour-suppressor panel is caller-supplied / fixed** (TP53, RB1, CDKN2A, PTEN,
  BRCA1, BRCA2). Arm membership is source-backed (NCBI Gene); the *choice* of panel is a registry
  curated list (analogous to the ONCO-CNA-002 oncogene panel), non-correctness-affecting for the
  deletion-detection logic itself — it only labels which arm maps to which gene name(s).

No source contradictions — cBioPortal (−2 = Deep Deletion = homozygous), Cheng et al. 2017 (total CN 0,
both alleles lost), CNVkit (`DeepDeletion` = integer CN 0), and NCBI Gene (tumour-suppressor cytobands)
corroborate one another. Cross-referenced: GISTIC2 (Mermel 2011) `t_del` and the ONCO-CNA-001/002 chain.
