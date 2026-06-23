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

| Unit(s) | Implemented | Fuller method not implemented | Note |
|---------|-------------|-------------------------------|------|
| ALIGN-MULTI-001 | Star (center-star) MSA, guide-tree progressive MSA (Feng–Doolittle), **and** iterative refinement (`MultipleAlignIterative`: MUSCLE-style tree-dependent restricted partitioning, Edgar 2004 — re-split the alignment at each guide-tree edge, realign the two profiles, keep on non-decreasing SP) | Full consistency-based refinement à la T-Coffee (consistency library / objective) | Iterative refinement removes the single-pass "once a gap, always a gap" limitation (early gap-placement errors are now corrected; refined SP is provably ≥ the progressive seed). Residual: the refinement is SP-guided and restricted to guide-tree edge partitions, not a full consistency-based optimizer. |
| RNA-STRUCT-001 | MFE-optimal structure via Zuker–Stiegler DP traceback (`CalculateMfeStructure`/`PredictStructureMfe`); greedy path retained as a fast heuristic | Pseudoknotted optima | The optimal structure is now returned by DP traceback (its energy equals the scalar MFE). **Residual:** the O(n³) recurrences are pseudoknot-free by construction, so crossing-helix (pseudoknotted) optima are not predicted — a different algorithm class. |
| PRIMER-TM-001 | Honest heuristic / NN-thermodynamics | — | NN-thermodynamics (ΔH/ΔS/Tm) **is** validated (SEQ-THERMO-001); the primer-design convenience scoring is heuristic. |

## 2. "Threshold / aggregation / framework" layers — they classify or combine caller-supplied inputs, they do not predict upstream

These units are correct implementations of a published **rule or formula**, but they sit *downstream*
of a model/measurement the caller must provide. They are **decision-support computations, not
validated clinical-grade predictors.**

| Unit | Computes (validated) | Caller must supply / not modelled |
|------|----------------------|-----------------------------------|
| ONCO-SIG-002/003/004 | NNLS signature **refitting**, bootstrap CIs, aetiology mapping | **De-novo signature extraction (NMF)** is not implemented — reference signatures are an input. |
| ONCO-MHC-001 | Strong/weak-binder **classification** by IC50 / %Rank cutoffs | The **affinity / %Rank prediction** (NetMHCpan-style learned model) — supplied by caller. |
| ONCO-HRD-001 | HRD score = LOH + TAI + LST, ≥42 cutoff | The three component scores are **inputs** (per-segment derivation is ONCO-LOH/CNA). |
| ONCO-PURITY/PLOIDY/CCF/CLONAL | Purity/ploidy/CCF formulas + clonal rule | Allele-specific CN segments, multiplicity, VAF — supplied. |
| EPIGEN-AGE-001 | Horvath `anti.trafo` + linear predictor | The **353-CpG coefficient table** is a caller-supplied input (framework, no fabricated coefficients). |
| ONCO-CHIP-001 | CHIP filter (gene panel + ≥2% VAF + WBC subtraction) | Origin call uses a **gene+VAF heuristic** where matched-WBC data is absent (over-removes vs strict matched-WBC origin). |
| ONCO-MRD-001 | ≥2-of-N positivity, IMAF, Poisson LoD | IMAF is read-pooled, **without** INVAR-style background subtraction / tumour-AF weighting. |

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
