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
| ALIGN-MULTI-001 | Full consistency-based MSA refinement (T-Coffee consistency library / objective) | A distinct objective class from the sum-of-pairs score the library optimizes. |
| RNA-STRUCT-001 | Pseudoknotted (crossing-helix) structure prediction | The O(n³) nearest-neighbour recurrences are pseudoknot-free by construction; pseudoknots need a different algorithm class. |
| PRIMER-TM-001 | A validated primer-design scoring objective | The convenience primer score is heuristic; there is no single canonical primer-score ground truth to assert against (NN-thermodynamics ΔH/ΔS/Tm itself is validated under SEQ-THERMO-001). |

## 2. "Threshold / aggregation / framework" layers — they classify or combine caller-supplied inputs, they do not predict upstream

These units are correct implementations of a published **rule or formula**, but they sit *downstream*
of a model/measurement the caller must provide. They are **decision-support computations, not
validated clinical-grade predictors.**

| Unit | Not modelled / caller-supplied |
|------|--------------------------------|
| ONCO-SIG-002/003/004 | Automatic NMF rank / model-stability selection (k is caller-specified); de-novo extraction uses the Frobenius, not the Poisson/KL, objective; extracted signatures are not auto-matched to COSMIC references. |
| ONCO-MHC-001 | The affinity / %Rank prediction (NetMHCpan-style learned model) — caller-supplied. |
| ONCO-PURITY/PLOIDY/CCF/CLONAL | Allele-specific CN segments, multiplicity and VAF — caller-supplied, not derived upstream. |
| ONCO-CHIP-001 | Strict matched-WBC origin calling — absent matched-WBC data, origin falls back to a gene+VAF heuristic (over-removes). |
| ONCO-MRD-001 | INVAR fragment-length (size) weighting, patient-specific outlier suppression and locus-noise filtering; the per-locus background error model is caller-supplied, not estimated from control plasma. |

## 3. Convention divergences (documented, internally consistent)

| Area | Convention used | Differs from |
|------|-----------------|--------------|
| Variant units (VARIANT-*, SV-*, ONCO somatic) | 0-based internal `Position` | VCF 1-based (1-based only on serialization, e.g. `ToVcfLines`). |
| GC / skew / composition outputs | Percentage (×100) in several units | Biopython fraction [0,1]. |
| Default thresholds | Configurable defaults that may differ from a named tool's published default (e.g. MSI 20% (MSIsensor2) vs MSIsensor ~3.5; ConsensusThreshold 0.5 vs Biopython 0.7; CNA inner cutoffs) | The named tool — but the value is sourced and the parameter is exposed. |
| DUST / compression complexity | `k=3` (DUST) and the LZ base-2 clamp are exactly sourced | `k≠3` is a documented, non-exact-asserted extrapolation. |
| Non-ACGT / ambiguous bases | Skipped or tracked as "Other" in several composition/thermo units | Tools that throw or partially count (e.g. Biopython S/W handling). |

---

## How to read this

- For **research / pipeline** use, these limitations are the normal scope boundaries of a from-first-
  principles library and are individually sourced.
- For **clinical or decision-grade** use, treat §2 especially as **heuristics / threshold layers that
  require an external validated predictor and clinical sign-off** — the library computes the rule, not
  the trained model behind it.
- Each row traces to a per-unit validation report (`git show <validate-commit>:docs/Validation/reports/<ID>.md`)
  and a `BY-DESIGN` / PASS-WITH-NOTES entry in the findings register.
