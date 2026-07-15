---
type: source
title: "Evidence: ONCO-FUSION-002 (HGNC gene-fusion designation + directional known-fusion lookup)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-FUSION-002-Evidence.md
sources:
  - docs/Evidence/ONCO-FUSION-002-Evidence.md
source_commit: ee968906158cef08ef55972fdafb4aad150e427a
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-FUSION-002

The validation-evidence artifact for test unit **ONCO-FUSION-002** — **Known Fusion Database
Lookup** (HGNC gene-fusion designation + caller-supplied known-fusion match). The **fifteenth
ingested unit of the Oncology family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is synthesized
in [[gene-fusion-nomenclature-known-fusion-lookup]]; it is the **annotation/naming** sibling of the
read-evidence caller [[gene-fusion-detection-read-evidence]] (ONCO-FUSION-001).
[[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (the HGNC gene-fusion nomenclature standard, mutually consistent):**
  - **Bruford et al. (2021)** — *HGNC recommendations for the designation of gene fusions*
    (*Leukemia* 35:3040–3043, PMC8550944, authority rank 2, official nomenclature + peer-reviewed):
    (1) **`::` double-colon separator** — unique, replacing the ambiguous hyphen/slash; (2) **5′
    partner always listed first**, before the `::`, irrespective of chromosomal location or gene
    orientation; (3) partners use **HGNC approved symbols**; (4) **read-through transcripts keep the
    hyphen** (`INS-IGF2`), `::` reserved for true fusions; (5) worked example — for
    `t(9;22)(q34.1;q11.2)`, **BCR** (chr 22) = 5′, **ABL1** (chr 9) = 3′ →
    `GetFusionAnnotation("BCR","ABL1") == "BCR::ABL1"`.
  - **Recommendations for future extensions to the HGNC gene fusion nomenclature** (*Leukemia*
    2021, PMC8632684, rank 2 consortium consensus) — endorses `::` as the cross-resource standard
    delimiter.

- **Documented corner cases / failure modes:**
  - **Direction matters** — the 5′-first rule makes `A::B` and `B::A` two *different* fusions (e.g.
    a reciprocal fusion); a directional lookup must not treat them as equal.
  - **Hyphen ≠ double colon** — `INS-IGF2` (read-through, hyphen) must not be confused with a fusion
    designation; the unit emits `::` for true fusions only.

- **Datasets (deterministic, derived from the cited nomenclature rules — no fusion DB is bundled):**
  - **HGNC worked example** — 5′ = BCR (chr 22), 3′ = ABL1 (chr 9), canonical `BCR::ABL1`,
    reciprocal (different fusion) `ABL1::BCR`.
  - **Caller-supplied known-fusion list** (illustrative; membership is caller-supplied, NOT
    fabricated/bundled) — `BCR::ABL1` → "Chronic myeloid leukemia driver"; `EML4::ALK` → "NSCLC
    driver, ALK TKI target".

- **Coverage recommendations (7 items):** MUST — `GetFusionAnnotation("BCR","ABL1") == "BCR::ABL1"`
  (5′ first, `::`); direction matters (`"ABL1::BCR" ≠ "BCR::ABL1"`); `MatchKnownFusions` returns the
  annotation for a designation present in the set, keyed `5′::3′`; no-match when only the reciprocal
  `B::A` is in the set; null/empty/null-argument rejected. SHOULD — case-insensitive symbol match
  (`bcr`/`abl1`). COULD — annotation round-trips through `MatchKnownFusions` on a `FusionCall`
  produced by `DetectFusions` (integration with ONCO-FUSION-001).

## Deviations and assumptions

- **ASSUMPTION — known-fusion membership is caller-supplied.** `MatchKnownFusions` takes the
  known-fusion set as a parameter; the library bundles no Mitelman/COSMIC/ChimerDB content. Only the
  **designation format** and the **directional 5′/3′ keying** are evidence-defined (Bruford et al.
  2021); set contents are the caller's responsibility. This makes the unit a **Framework**
  algorithm (format/keying source-backed, data supplied).
- **ASSUMPTION — symbol case-insensitivity.** Lookup compares symbols case-insensitively
  (ordinal-ignore-case) while preserving directionality. Not contradicted by the source; flagged
  because Bruford et al. specify approved (uppercase) symbols without addressing case folding — it
  affects matching only, not the formal `::`/order rule.

No source contradictions — Bruford et al. (2021) and the consortium extension letter agree on the
`::` separator and the 5′-first directional rule.
