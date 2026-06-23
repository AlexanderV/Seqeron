# Seqeron.Genomics — Validated Limitations & Operating Envelope

**Date:** 2026-06-16   **Library:** Seqeron.Genomics (mission-critical)
**Basis:** the independent Phase-1 + Phase-2 validation campaign (234 ☑ Registry units;
see [VALIDATION_LEDGER.md](VALIDATION_LEDGER.md), [FINDINGS_REGISTER.md](FINDINGS_REGISTER.md),
and `docs/Validation/reports/`).

This document consolidates **documented, sourced limitations** — places where an algorithm is a
faithful but *simplified / subset* realisation of a fuller published method, or where it consumes a
caller-supplied input rather than computing it upstream. **None of these are defects:** every item
below was validated as correct *for its stated contract* and is `BY-DESIGN` in the findings register.
The point of this file is to make the library's honest operating envelope visible in one place so
downstream users are not surprised.

Genuine algorithm defects found during validation were **all fixed in-session** (13 in Phase 2);
they are recorded in the ledger, not here. Deferred *enhancements* that would lift some of these
limitations are tracked separately as the "Deferred BIG fixes" backlog in the ledger and the
`C. NOT-POSSIBLE (radical)` section of the findings register.

---

## 1. Algorithmic simplifications (a fuller method exists; the simpler one is correct for its scope)

| Unit(s) | Not implemented | Note |
|---------|-----------------|------|
| RNA-STRUCT-001 | Pseudoknot classes outside the pknotsRG canonical csr-PK grammar (kissing hairpins, triple-crossing / chained "complex" helix interactions, non-canonical bulged or unequal-length helices); tertiary-stabilised knots as the MFE structure | Classes Reeder & Giegerich (2004) explicitly exclude from csr-PK; and tertiary-stabilised knots (e.g. BWYV / PDB 437D) are not recoverable by *any* nearest-neighbour thermodynamic model — an energy-model floor, not an algorithm gap. (csr-PK incl. nested/multiple/over-arching knots is handled — RNA-PKPREDICT-001 / RNA-PKRECURSIVE-001.) |

## 2. "Threshold / aggregation / framework" layers — they classify or combine caller-supplied inputs, they do not predict upstream

These units are correct implementations of a published **rule or formula**, but they sit *downstream*
of a model/measurement the caller must provide. They are **decision-support computations, not
validated clinical-grade predictors.**

| Unit | Not modelled / caller-supplied |
|------|--------------------------------|
| ONCO-MHC-001 | The affinity / %Rank prediction (NetMHCpan-style learned model) — caller-supplied. |

## ~~3. Convention divergences~~ — retired (all addressed via opt-in modes or by-design configurability)

The former §3 "Convention divergences" rows were retired on 2026-06-23. None represented a defect;
each was either made caller-selectable via an opt-in mode (defaults unchanged) or was already
by-design configurability. For the record:

- **Variant 0-based vs VCF 1-based** — already handled at serialization (`ToVcfLines` emits
  `Position + 1`). A convenience `Variant.VcfPosition` (1-based) accessor was added for callers that
  consume the record directly (VCF v4.3 §1.4.1). No default change.
- **GC / skew / composition percentage vs Biopython fraction [0,1]** — added an opt-in `fraction`
  parameter to `GcSkewCalculator.AnalyzeGcContent` and `SequenceStatistics.CalculateGcContentProfile`
  (and a `CalculateGcFraction` already exists on `SequenceExtensions`), so callers can request
  Biopython's [0,1] convention. Default percentage output unchanged.
- **Default thresholds differing from a named tool's default** — by-design: every such threshold is
  already exposed as a configurable parameter and the named tool's value is sourced; the caller can
  set it. Not a code gap.
- **DUST `k=3` / LZ base-2 clamp** — by-design: `k=3` and the clamp are exactly sourced and `wordSize`
  is already a caller parameter; `k≠3` was only a *test-assertion-scope* note (exact values asserted
  at k=3), not a behaviour gap. Nothing to fix.
- **Non-ACGT / ambiguous bases skipped or "Other"** — added an opt-in
  `SequenceExtensions.CalculateGcFraction(GcAmbiguityMode)` (Remove/Ignore/Weighted) that reproduces
  Biopython `gc_fraction`'s IUPAC handling (S counts as GC, W as length, weighted `_gc_values`). The
  existing default (A/T/G/C/U only) is unchanged.

---

## How to read this

- For **research / pipeline** use, these limitations are the normal scope boundaries of a from-first-
  principles library and are individually sourced.
- For **clinical or decision-grade** use, treat §2 especially as **heuristics / threshold layers that
  require an external validated predictor and clinical sign-off** — the library computes the rule, not
  the trained model behind it.
- Each row traces to a per-unit validation report (`git show <validate-commit>:docs/Validation/reports/<ID>.md`)
  and a `BY-DESIGN` / PASS-WITH-NOTES entry in the findings register.
