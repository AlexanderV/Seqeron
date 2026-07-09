---
type: concept
title: "Gene-fusion HGNC designation + directional known-fusion lookup"
tags: [oncology, nomenclature, algorithm]
sources:
  - docs/Evidence/ONCO-FUSION-002-Evidence.md
source_commit: ee968906158cef08ef55972fdafb4aad150e427a
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:gene-fusion-detection-read-evidence
      source: onco-fusion-002-evidence
      evidence: "COULD Test: annotation string round-trips through MatchKnownFusions on a FusionCall produced by DetectFusions — integration with ONCO-FUSION-001 (ONCO-FUSION-002-Evidence.md Coverage item 7)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-fusion-002-evidence
      evidence: "Test Unit ID: ONCO-FUSION-002, Algorithm: Known Fusion Database Lookup (HGNC gene-fusion designation + caller-supplied known-fusion match)"
      confidence: high
      status: current
---

# Gene-fusion HGNC designation + directional known-fusion lookup

The **fifteenth ingested Oncology unit** (ONCO-FUSION-002) and the wiki's **fusion
annotation / naming** method — distinct from the read-evidence caller
[[gene-fusion-detection-read-evidence]] (ONCO-FUSION-001), which decides *whether* a fusion is
detected. This unit does not touch reads at all: it **formats** a fusion's canonical HGNC
designation and **matches** a fusion against a caller-supplied known-fusion set. Validated under
test unit **ONCO-FUSION-002** ([[onco-fusion-002-evidence]]); [[test-unit-registry]] tracks the
unit and [[algorithm-validation-evidence]] describes the artifact pattern.

## HGNC designation format (`GetFusionAnnotation`)

The canonical designation is defined by the **HGNC recommendation** (Bruford et al. 2021, *Leukemia*):

- **Double-colon separator** `::` — a unique, instantly recognizable separator (`BCR::ABL1`),
  chosen to replace the previously ambiguous hyphen (`-`) and forward slash (`/`).
- **Directional 5′-first order** — the **5′ partner is always written first**, before the `::`,
  irrespective of chromosomal location or gene orientation; the 3′ partner follows.
- **HGNC approved gene symbols** designate the partners.

```
GetFusionAnnotation(gene5p, gene3p) = gene5p + "::" + gene3p
```

Worked HGNC example: for `t(9;22)(q34.1;q11.2)`, **BCR** (chr 22) is the 5′ gene and **ABL1**
(chr 9) is the 3′ gene, so `GetFusionAnnotation("BCR","ABL1") == "BCR::ABL1"`.

## Directional known-fusion match (`MatchKnownFusions`)

`MatchKnownFusions` looks a fusion up in a **caller-supplied** known-fusion set, keyed by the
directional `5′::3′` designation, and returns the caller's annotation label (e.g. *"Chronic
myeloid leukemia driver"* for `BCR::ABL1`, *"NSCLC driver, ALK TKI target"* for `EML4::ALK`).

This is a **Framework algorithm**: only the **designation format** and the **directional 5′/3′
keying** are evidence-defined (Bruford et al. 2021). The **set contents are the caller's
responsibility** — the library does **not** bundle Mitelman, COSMIC, or ChimerDB content
(licensing + curation are out of scope). Symbol matching is **case-insensitive**
(ordinal-ignore-case) while preserving directionality.

## Invariants and corner cases

- **Direction matters.** Because the 5′ gene is always first, `A::B` and `B::A` describe two
  *different* fusions (e.g. a reciprocal fusion). A directional lookup must not treat them as
  equal — `GetFusionAnnotation("ABL1","BCR") == "ABL1::BCR" ≠ "BCR::ABL1"`, and a set containing
  only the reciprocal `B::A` yields **no match** for `A::B`.
- **Hyphen ≠ double colon.** Read-through transcripts keep the **hyphen** (`INS-IGF2`); `::` is
  reserved for true fusion genes. This unit produces `::` designations only.
- **Case-insensitive symbols.** `bcr`/`abl1` match `BCR`/`ABL1` (assumption: HGNC symbols are
  formally uppercase, but real inputs vary in case; folding affects matching only, not the
  format/order rule).
- **Input validation.** Null/empty symbols and null arguments are rejected (mirrors sibling
  methods' contract).

## Scope and relation to the read-evidence caller

A [[research-grade-limitations|research-grade]] method, **not for clinical use**. It is the
**naming / annotation** layer that sits downstream of the [[gene-fusion-detection-read-evidence|
read-evidence fusion caller]]: a `FusionCall` produced by `DetectFusions` (ONCO-FUSION-001) can be
annotated via `GetFusionAnnotation` and round-tripped through `MatchKnownFusions` (the unit's
integration COULD-test). It is orthogonal to the copy-number / clonal-structure ONCO units and to
the clinical-interpretation units ([[clinical-actionability-oncokb-levels]],
[[cancer-variant-tier-classification-amp-asco-cap]]) that would consume a matched driver fusion.
Distinct also from the premature-stop / transcript-reconstruction scope deferred to
ONCO-FUSION-003. No source contradictions.
