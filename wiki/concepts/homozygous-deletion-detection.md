---
type: concept
title: "Homozygous / deep deletion detection (total-CN-0 call + tumour-suppressor mapping)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-CNA-003-Evidence.md
source_commit: 819918712f8e6a3fddb0f4a534fb6f69bc24cf5b
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-cna-003-evidence
      evidence: "Test Unit ID: ONCO-CNA-003 ... Algorithm: Homozygous (Deep) Deletion Detection"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:copy-number-alteration-classification
      source: onco-cna-003-evidence
      evidence: "A homozygous deletion is exactly the existing CN-0 (DeepDeletion) classification: 'integer CN 0 ⇒ DeepDeletion' (CNVkit, ONCO-CNA-001); no new numeric threshold is invented."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:focal-amplification-detection
      source: onco-cna-003-evidence
      evidence: "The deletion counterpart of ONCO-CNA-002: an order-preserving segment filter (CN-0 homozygous deletions) plus an arm→gene panel mapping — IdentifyDeletedTumorSuppressors mirrors IdentifyAmplifiedOncogenes."
      confidence: high
      status: current
---

# Homozygous / deep deletion detection (total-CN-0 call + tumour-suppressor mapping)

The **oncology deep-deletion layer**: given per-segment copy-number data it (a) keeps only segments that
are **homozygous deletions** — total (integer) copy number **exactly 0** — then (b) **maps** each such
segment's arm label to a panel of known **tumour-suppressor** genes. `IdentifyDeletedTumorSuppressors`
is the deletion mirror of ONCO-CNA-002's `IdentifyAmplifiedOncogenes`. Validated under test unit
**ONCO-CNA-003**; the literature-traced record is [[onco-cna-003-evidence]], [[test-unit-registry]]
tracks the unit, and [[algorithm-validation-evidence]] describes the evidence-artifact pattern.

**How it differs from its neighbours (the deletion sibling of the CNA family):**

- [[copy-number-alteration-classification]] (ONCO-CNA-001) already bins a log2 ratio into the five
  discrete states, and **integer CN 0 = `DeepDeletion`** is one of them. This unit does **not** invent a
  new threshold: a segment is a homozygous deletion **iff its classified integer copy number is 0**. It
  **reuses** CNA-001's CN-0 call and adds the tumour-suppressor mapping on top.
- [[focal-amplification-detection]] (ONCO-CNA-002) is the **amplification** counterpart — an
  order-preserving segment filter + arm→oncogene panel. This unit is the **mirror on the loss side**:
  filter to CN-0 homozygous deletions, map arm→tumour-suppressor. (One structural difference: CNA-002's
  focal test also applies a **length** gate — segment < 98% of arm; the deletion call here is defined by
  **depth** — total CN 0 — per its cBioPortal / Cheng-et-al. sources, not by a length ratio.)
- [[allele-specific-copy-number-ascat]] (ONCO-ASCAT-001) derives allele-specific nA/nB and could express
  "one allele lost" (hemizygous) vs "both lost" (homozygous) directly; this unit works from the
  **total-CN** classification only — homozygous = total CN 0.

## 1. The homozygous-deletion predicate

A segment is a **homozygous (deep) deletion** iff its classified integer copy number is 0:

```
homozygous deletion:  round/threshold(log2) → integer CN == 0     # == DeepDeletion state (ONCO-CNA-001)
```

Source convergence on the CN-0 definition:

| Source | Definition of homozygous deletion |
|--------|-----------------------------------|
| cBioPortal (file format / FAQ) | discrete value **−2** = "Deep Deletion … possibly a homozygous deletion" (deepest loss) |
| Cheng et al. 2017 (Nat Commun) | **total copy number 0** — "zero copies of **both alleles**"; two independent hits |
| CNVkit `absolute_threshold` (via ONCO-CNA-001) | integer CN **0** ⇒ `DeepDeletion` (log2 ≤ −1.1, default thresholds) |

**Not a homozygous deletion:** a single-copy loss (cBioPortal **−1** "Shallow Deletion", heterozygous /
hemizygous, one allele remaining, CN ≥ 1) is a *loss*, never a homozygous deletion. Neutral / gain /
amplification segments are excluded a fortiori.

## 2. Tumour-suppressor mapping (`IdentifyDeletedTumorSuppressors`)

Each homozygous-deletion segment carries an **arm label** (chromosome + p/q arm, e.g. `17p`). The mapper
matches that arm prefix against a small registry of recurrently deleted **tumour suppressors** and their
cytogenetic locations (NCBI Gene):

| Tumour suppressor | Arm | Cytoband |
|-------------------|-----|----------|
| TP53 | 17p | 17p13.1 |
| RB1 | 13q | 13q14.2 |
| CDKN2A | 9p | 9p21.3 |
| PTEN | 10q | 10q23.31 |
| BRCA1 | 17q | 17q21.31 |
| BRCA2 | 13q | 13q13.1 |

A single arm can carry **multiple** tumour suppressors (**13q → both RB1 and BRCA2**). Only homozygous
deletions feed the mapper — a loss / neutral / gain segment never yields a deleted tumour suppressor.
The arm→gene panel is the algorithm's built-in curated registry (analogous to CNA-002's oncogene panel),
not a caller-supplied knowledgebase.

## Worked dataset (default thresholds −1.1/−0.25/0.2/0.7)

| Segment | Arm | log2 | Integer CN | State | Homozygous deletion? | Maps to |
|---------|-----|------|-----------|-------|----------------------|---------|
| A | 17p | −2.0 | 0 | DeepDeletion | **yes** | TP53 |
| B | 13q | −2.0 | 0 | DeepDeletion | **yes** | RB1 + BRCA2 |
| C | 10q | −0.5 | 1 | Loss | no (single-copy) | — |
| D | 9p | 0.0 | 2 | Neutral | no | — |

## Corner cases and failure modes

- **Shallow vs deep:** a single-copy / heterozygous loss (−1, CN ≥ 1) is NOT a homozygous deletion —
  only total CN 0 qualifies (segment C).
- **Boundary at the deletion cutoff (log2 = −1.1):** inclusive `log2 <= thresh` → CN 0 → homozygous
  (mirrors CNA-001's boundary-inclusive lower-bin rule).
- **Custom thresholds move the CN-0 boundary:** raising the deletion cutoff can reclassify a previously
  CN-1 segment as CN 0 (CNVkit `absolute_threshold` thresholds are parameters).
- **Order-preserving filter:** each qualifying segment is reported once, input order preserved (mirror
  `DetectFocalAmplifications`).
- **Validation:** empty → empty; null → `ArgumentNullException`; invalid arm length / `End ≤ Start` →
  `ArgumentException` (mirror sibling ONCO-CNA-002 validation).

## Assumptions and scope

- **ASSUMPTION — homozygous deletion identified at the integer-CN level via the existing CN-0
  (`DeepDeletion`) classification.** cBioPortal (−2) and Cheng et al. (total CN 0) converge on the call
  the repository already realizes (CNVkit, ONCO-CNA-001). **No new numeric threshold is invented.**
- **ASSUMPTION — curated tumour-suppressor panel is caller-supplied / fixed** (TP53, RB1, CDKN2A, PTEN,
  BRCA1, BRCA2). Arm membership is source-backed (NCBI Gene); the *choice* of panel is a registry list,
  non-correctness-affecting for the detection logic — it only labels arm→gene name(s).

A [[scientific-rigor|research-grade]] correctness reference — **not for clinical or diagnostic use.**
No source contradictions: cBioPortal (−2 = Deep Deletion = homozygous), Cheng et al. 2017 (total CN 0,
both alleles lost), CNVkit (`DeepDeletion` = integer CN 0), and NCBI Gene (tumour-suppressor cytobands)
corroborate one another across the whole predicate.
