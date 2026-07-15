---
type: source
title: "Evidence: PROTMOTIF-COMMON-001 (Common motif finding — built-in PROSITE-pattern dictionary scan)"
tags: [validation, protein, motif]
doc_path: docs/Evidence/PROTMOTIF-COMMON-001-Evidence.md
sources:
  - docs/Evidence/PROTMOTIF-COMMON-001-Evidence.md
source_commit: 12b13e4ecc31636e0c27a2c4b0098bf11d6cc054
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PROTMOTIF-COMMON-001

The validation-evidence artifact for test unit **PROTMOTIF-COMMON-001** — **common motif
finding** (`ProteinMotifFinder.FindCommonMotifs`): scan an amino-acid sequence against a fixed
built-in dictionary (`CommonMotifs`) of canonical PROSITE patterns and aggregate every hit, each
carrying the matching entry's accession/name. It is one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the model, catalog, contract,
invariants and worked oracles are synthesized in [[common-protein-motifs]]. See
[[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources (all authority rank 2 — official database / spec):**
  - **PROSITE PS00001** `ASN_GLYCOSYLATION` — pattern `N-{P}-[ST]-{P}` (N-glycosylation site;
    `{P}` = any residue except Proline).
  - **PROSITE PS00005** `PKC_PHOSPHO_SITE` — `[ST]-x-[RK]` (protein kinase C phosphorylation site).
  - **PROSITE PS00006** `CK2_PHOSPHO_SITE` — `[ST]-x(2)-[DE]` (casein kinase II phosphorylation site).
  - **PROSITE PS00016** `RGD` — `R-G-D` (cell-attachment / RGD sequence).
  - **PROSITE PS00017** `ATP_GTP_A` — `[AG]-x(4)-G-K-[ST]` (ATP/GTP-binding P-loop, motif A).
  - **ScanProsite / PROSITE User Manual** — pattern syntax: IUPAC one-letter code; `x` = any
    residue; `[ ]` = allowed set; `{ }` = excluded set; `-` = element separator; `x(3)` = fixed
    repetition; `x(2,4)` = variable repetition; default reporting is "greedy, overlaps, no
    includes" (overlaps reported unless one match is fully contained in another); hits reported as
    `[start]-[stop]` in **1-based inclusive** coordinates.
- **Documented corner cases:** Proline at an excluded `{P}` position rejects an N-glycosylation
  match (`N-P-[ST]` is not a hit); overlapping occurrences of a pattern are both reported.
- **Test dataset (synthetic windows, 0-based inclusive Start..End as `MotifMatch` records):**
  `AAAANFTAAAA` → PS00001 `4..7 NFTA`; `AAAANPSAAAAANPTAAA` → PS00001 no match; `AAAAASARKAAA`
  → PS00005 `5..7 SAR`; `AAAASAAEASDEDAAA` → PS00006 `4..7 SAAE` + `9..12 SDED`;
  `AAAAAGXXXXGKSAAAA` → PS00017 `5..12 GXXXXGKS`; `AARGDKK` → PS00016 `2..4 RGD`;
  `RGDRGD` → PS00016 `0..2 RGD` + `3..5 RGD` (two non-overlapping).

## Deviations and assumptions

No algorithm deviations. **One assumption**, an API-shape convention: PROSITE/ScanProsite report
**1-based** inclusive coordinates while the repository `MotifMatch` records **0-based**
`Start`/`End` (matching sibling units PROTMOTIF-FIND-001 and PROTMOTIF-PATTERN-001). This is not a
correctness-affecting parameter — matched substring content and relative positions are identical;
only the coordinate origin differs. Tests assert 0-based positions per the repository convention.

## Recommended coverage

MUST: find each canonical PROSITE motif (PS00001/00005/00006/00016/00017) at the exact 0-based
position with the exact substring; Proline exclusion `{P}` rejects `N-P-[ST]` windows; multiple
distinct pattern types aggregated from one whole-dictionary scan; overlapping occurrences both
reported; null/empty → empty. SHOULD: substring invariant `Sequence == protein.Substring(Start,
End-Start+1)`; `MotifName`/`Pattern` carry the matching entry's identity. COULD: determinism.
