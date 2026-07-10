---
type: source
title: "LIMITATIONS.md — the validated operating envelope and open scope boundaries"
tags: [validation, governance, limitations]
doc_path: docs/Validation/LIMITATIONS.md
source_commit: 45545719fbdd7689c20bb680104862f6098adf32
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# LIMITATIONS.md — the validated operating envelope and open scope boundaries

The project's operating-envelope document (last reviewed 2026-06-26): the single human-readable
catalog of what Seqeron.Genomics **currently does NOT do**. It is the reference behind the runtime
`LimitationPolicy` and the `bio-rigor` discipline — the mechanism and taxonomy are written up as the
concept [[operating-envelope-and-limitation-policy]]. This page summarizes the document; see that
concept for how the envelope is enforced.

## What the document lists (and deliberately excludes)

It lists **only genuine open limitations**. Anything already implemented — opt-in alternatives,
convention-parity modes, resolved defects — is excluded and recorded instead in the per-unit reports
under `docs/Validation/reports/{UNIT}.md`. Crucially, **every** limitation is `BY-DESIGN` and its unit
is `✅ CLEAN` for its stated contract: these are honest scope boundaries, never defects. Pure
**performance** optimisations (e.g. a seeded index where an exhaustive scan already returns the same
result) are explicitly not limitations and are not listed.

## Three kinds of limitation

- **Irreducible** — no algorithm or model can close it (physics / information theory).
- **Data-blocked** — needs a trained model / matrix / database that is gated, non-redistributable, or
  never measured / published. Reopens only when that data becomes available (several units already
  accept it via a caller-supplied loader).
- **Scope** — a deliberate out-of-scope boundary; use the named reference tool, or supply the input.

## Enumerated limitations (by unit)

- **Irreducible.** PARSE-FASTQ-001 (a FASTQ wholly confined to the Phred+33/+64 overlap range ASCII
  64–74 is information-theoretically ambiguous; `DetectEncoding` flags it `Ambiguous`, defaults +33);
  RNA-STRUCT-001 twice (tertiary-stabilised knots like BWYV are not representable by *any*
  nearest-neighbour energy model; pseudoknot classes outside the csr-PK grammar need an
  never-measured loop–loop interaction energy).
- **Data-blocked.** ONCO-MHC-001 (proprietary NetMHCpan-4.1 ANN / MHCflurry; no redistributable
  SMM/BIMAS matrix); ONCO-IMMUNE-001 (CIBERSORT LM22 matrix is no-redistribution / registration-gated;
  ESTIMATE absolute-purity transform is Affymetrix/ABSOLUTE-calibrated); META-BIN-001 (gated
  `checkm_data` DB + TIGRFAM markers → only domain-level completeness/contamination ships);
  CHROM-CENT-001 (separating SF1 from SF2 + chromosome-specific HOR identity needs the SF-resolved
  consensus-monomer library, published only in no-licence supplements/HMM repos); DISORDER-REGION-001
  (no predictor publishes a calibrated per-region disorder **confidence**; the boundaries themselves
  follow the validated TOP-IDP threshold).
- **Scope.** PROBE-DESIGN-001 (quantitative MGB ΔTm is empirical/proprietary, no closed form; use a
  chemistry-specific tool); MIRNA-TARGET-001 (full context++ needs caller-supplied 3'UTR set /
  alignment+tree+sigmoid params / SPS/ORF features); MIRNA-CLEAVAGE-001 (miRBase 3p/star boundaries
  encode sequencing-read cut sites, not a deterministic fold + fixed-overhang rule; supply MIMAT
  coordinates or read pileups); MIRNA-PRECURSOR-001 (miRDeep2 read-stacking log-odds exists only in a
  gated supplement + GPL source; use miRDeep2 with the caller's reads).

## Research vs clinical framing (the disclaimer)

The document's "How to read this" section states the two-audience stance directly: for **research /
pipeline** use, every row is a normal scope boundary of a from-first-principles library — the algorithm
is validated correct for its stated contract. For **clinical / decision-grade** use, the §2 oncology
rows (ONCO-MHC-001, ONCO-IMMUNE-001) mark layers that require an external validated predictor and
clinical sign-off — the library computes the *rule*, not the trained model behind it. This is the
document that operationalizes the project's [[research-grade-limitations|not-for-clinical-use]] status
at the per-unit level, and the honest-scope backing for [[scientific-rigor]].

## Where this fits

- [[operating-envelope-and-limitation-policy]] — the concept: the `LimitationPolicy` three-mode guard,
  the minimum-access-mode table, and the taxonomy this document catalogs.
- [[scientific-rigor]] — runtime honesty; `LimitationPolicy` guards each algorithm to the envelope this
  document defines.
- [[research-grade-limitations]] — the project-level beta / not-for-clinical disclaimer this document
  makes concrete unit-by-unit.
- [[validation-and-testing]] — the per-unit validation campaign; each row here traces to a `✅ CLEAN`
  report under `docs/Validation/reports/`.
