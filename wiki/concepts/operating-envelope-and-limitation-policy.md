---
type: concept
title: "Operating envelope and LimitationPolicy"
tags: [rigor, governance, limitations]
sources:
  - docs/Validation/LIMITATIONS.md
source_commit: 45545719fbdd7689c20bb680104862f6098adf32
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:scientific-rigor
      source: limitations
      evidence: "LimitationPolicy guards each algorithm's validated scope; LIMITATIONS.md is the honest scope it documents — the runtime-enforcement mechanism named in scientific-rigor."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:validation-and-testing
      source: limitations
      evidence: "Each limitation row is BY-DESIGN and traces to a ✅ CLEAN per-unit report under docs/Validation/reports/; the envelope is the published output of the validation campaign."
      confidence: high
      status: current
---

# Operating envelope and LimitationPolicy

The **validated operating envelope** is the boundary between what each Seqeron algorithm is validated to
compute correctly and what it deliberately does not do. It has two faces: a human-readable catalog of
open limitations (`docs/Validation/LIMITATIONS.md`, summarized in [[limitations]]) and a **runtime
guard**, `Seqeron.Genomics.Core.LimitationPolicy`, that refuses or flags a call that steps outside the
envelope. This is the concrete mechanism named — but not detailed — in [[scientific-rigor]].

## The runtime guard: three modes

`LimitationPolicy` has three modes, least → most permissive: **`Strict` < `Moderate` < `Permissive`**
(default **`Moderate`**). A guarded call throws `SeqeronLimitationException` — naming the limitation,
what it relates to, and how to obtain the result another way — when the effective mode is *more
restrictive* than that limitation's **minimum access mode**. Set it globally via
`LimitationPolicy.DefaultMode`, or scope a region with `using (LimitationPolicy.Use(mode)) { … }`.

- **`Strict`** — only the ideal *and complete* result; throws on every guarded branch.
- **`Moderate`** (default) — throws on **non-ideal-output** branches; **allows** the
  **correct-but-incomplete / narrower-contract** branches.
- **`Permissive`** — allows everything (historical best-effort).

Each guarded unit declares a **minimum access mode** (single source of truth: `LimitationCatalog`):

- **`Permissive` min** (non-ideal output — blocked in Strict *and* Moderate): PARSE-FASTQ-001 (encoding
  undetermined), CHROM-CENT-001 (`Sf1OrSf2Dimeric`), DISORDER-REGION-001 (uncalibrated confidence),
  MIRNA-TARGET-001 (partial context++), MIRNA-CLEAVAGE-001 (approximate 3p/star span).
- **`Moderate` min** (correct-but-incomplete — blocked only in Strict): ONCO-MHC-001 (SMM/BIMAS matrix
  score), ONCO-IMMUNE-001 (ABIS/matrix deconvolution & ESTIMATE purity), META-BIN-001 (domain-level
  CheckM), PROBE-DESIGN-001 (qualitative MGB rules).

Some limitations have **no runtime guard** (documented only): the irreducible RNA-STRUCT-001 pair (the
result is exact for the stated energy model / grammar and the shortfall is undetectable per call) and
MIRNA-PRECURSOR-001 (read-stacking is not implemented, so nothing is returned to gate).

## The limitation taxonomy

The envelope classifies every open limitation into one of three kinds, which also governs whether it can
ever be closed:

- **Irreducible** — no algorithm can close it (physics / information theory). Never reopens.
- **Data-blocked** — needs a gated / non-redistributable / never-measured model, matrix, or database.
  Reopens when that data arrives; several units already accept it via a caller-supplied loader.
- **Scope** — a deliberate out-of-scope boundary; use the named reference tool or supply the input.

Every listed limitation is `BY-DESIGN` with its unit `✅ CLEAN` for its stated contract — these are
honest scope boundaries, never defects. The full per-unit enumeration lives in [[limitations]].

## Why it matters

This is the governance artifact that makes [[scientific-rigor]] enforceable rather than aspirational: the
assistant cannot silently return a result outside an algorithm's validated scope. It is the per-unit,
runtime-checked face of the project-level [[research-grade-limitations]] disclaimer, and the published
scope output of the [[validation-and-testing]] campaign (each row traces to a `✅ CLEAN` report). The
`bio-rigor` skill enforces the same envelope at the discipline layer.
