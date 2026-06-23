# Seqeron.Genomics — Validated Limitations & Operating Envelope

**Library:** Seqeron.Genomics (mission-critical)   **Last reviewed:** 2026-06-23

This file lists the library's **current, sourced limitations** — the only places where an algorithm is
a faithful but *simplified / subset* realisation of a fuller published method, or where it consumes a
caller-supplied input rather than computing it upstream. Each is `BY-DESIGN` and validated correct for
its stated contract; the purpose is to make the honest operating envelope visible in one place.

---

## 1. Algorithmic simplification (a fuller method exists; the implemented one is correct for its scope)

| Unit | Not implemented | Note |
|------|-----------------|------|
| RNA-STRUCT-001 | Pseudoknot classes outside the pknotsRG canonical csr-PK grammar (kissing hairpins, triple-crossing / chained "complex" helix interactions, non-canonical bulged or unequal-length helices); tertiary-stabilised knots as the MFE structure | Classes Reeder & Giegerich (2004) explicitly exclude from csr-PK; and tertiary-stabilised knots (e.g. BWYV / PDB 437D) are not recoverable by *any* nearest-neighbour thermodynamic model — an energy-model floor, not an algorithm gap. |

## 2. Threshold / aggregation layer (classifies a caller-supplied input; does not predict it upstream)

| Unit | Not modelled / caller-supplied |
|------|--------------------------------|
| ONCO-MHC-001 | The peptide–MHC affinity / %Rank prediction (NetMHCpan-style learned model) is caller-supplied: the library classifies a supplied IC50 / %Rank into strong/weak binder but does not bundle a trained predictor (no redistributable, cross-verifiable model available). |

---

## How to read this

- For **research / pipeline** use, these are the normal scope boundaries of a from-first-principles library.
- For **clinical / decision-grade** use, treat §2 as a threshold layer that requires an external validated
  predictor and clinical sign-off — the library computes the rule, not the trained model behind it.
- Each row traces to its per-unit validation report under `docs/Validation/reports/` and a `BY-DESIGN`
  entry in the findings register.
