---
type: source
title: "Validation report: DISORDER-PRED-001 (intrinsic disorder prediction — normalized TOP-IDP sliding window, DisorderPredictor.PredictDisorder)"
tags: [validation, analysis]
doc_path: docs/Validation/reports/DISORDER-PRED-001.md
sources:
  - docs/Validation/reports/DISORDER-PRED-001.md
source_commit: 920bd895c63f51453f3d267e06abb5c61d0b6fc3
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: DISORDER-PRED-001

The two-stage **validation write-up** for test unit **DISORDER-PRED-001** — the **core per-residue
intrinsic-disorder predictor** (`DisorderPredictor.PredictDisorder`), which scores each residue as the
sliding-window mean of the **normalized TOP-IDP propensity scale** and thresholds at the published
0.542 cutoff — validated 2026-06-24. This is the *report* artifact that feeds one row of the
[[validation-ledger]]; it records the validator's independent **verdict** on both the algorithm
description (Stage A) and the shipped code (Stage B), and the wider campaign is
[[validation-and-testing]]. The scale, cutoff, Dunker classification, invariants and oracles are
synthesized on the concept [[intrinsic-disorder-prediction-top-idp]] (the shared `PredictDisorder`
anchor the whole protein-disorder family reads from); [[test-unit-registry]] defines the unit.
Distinct from [[disorder-pred-001-evidence]] — the pre-implementation evidence artifact sourced from
`docs/Evidence/` — this page is the independent two-stage re-validation verdict. Sibling reports
[[disorder-lc-001-report]] (SEG low-complexity) and [[disorder-morf-001-report]] (MoRF) cover
different units of the same family.

## Verdict

**Stage A: PASS · Stage B: PASS · End state: CLEAN.** No algorithm defect, no code defect, **no code
changed this session**. Full `DisorderPredictor` test family **113/113 green**, of which the focused
`DisorderPredictor_DisorderPrediction_Tests` is **22/22**. No divergences at either stage.

## Canonical methods & source under test

- `DisorderPredictor.PredictDisorder` (`DisorderPredictor.cs:190`) →
  `CalculatePerResidueScores` (`:227`, centered window truncated at termini) →
  `CalculateDisorderScore` (`:255`, normalize `(prop − TopIdpMin)/TopIdpRange` and average); threshold
  applied as `score >= disorderThreshold` (`:242`). Bundled utility `CalculateHydropathy`
  (Kyte & Doolittle 1982). Constants: `TopIdpMin = −0.884` (W), `TopIdpMax = 0.987` (P),
  `TopIdpRange = 1.871`, `TopIdpCutoff = 0.542`, default window **21**.

## Method scope (honest scoping, from the report)

Per-residue disorder is an explicit **composition heuristic**, not an energy-based (IUPred) or ML
predictor. The report confirms the source XML doc comments scope this honestly — single-feature,
AUC ≈ 0.65–0.72 vs IUPred2A 0.75–0.80 and flDPnn 0.85–0.90 — and recommend dedicated tools for
publication-grade work. **No false IUPred-grade claim.** A [[research-grade-limitations|research-grade]]
implementation.

## Stage A — description (algorithm faithfulness)

Sources re-opened live this session: **Campen et al. 2008 (TOP-IDP, PMC2676888)** via WebFetch of the
PMC full text — Table 2 values confirmed verbatim, cutoff `0.542`, prediction equation
`I = −(<TOP-IDP> − 0.542)` (positive ⟹ ordered, so `<TOP-IDP> > 0.542` ⟹ disordered), scales
normalized to min 0 / max 1, window 21; and **Wikipedia "Intrinsically disordered proteins" + Dunker
2001 (PMID 11381529)** confirming the disjoint 8+8+4 classification (order {W,C,F,I,Y,V,L,N},
disorder {A,R,G,Q,S,P,E,K}, ambiguous {H,M,T,D}) — matching the code exactly.

- **Formula check.** Per-residue `S(aa) = (TOP-IDP(aa) − TopIdpMin)/TopIdpRange` with
  `TopIdpMin = −0.884`, `TopIdpRange = 1.871`; window score = mean of `S(aa)` over the truncated
  centered window; `score >= 0.542` ⟹ disordered. The code applies 0.542 to the **normalized** averaged
  score, consistent with Campen's all-scales-normalized convention for the published cutoff.
- **All 20 TOP-IDP values match** the implementation exactly (anchors W/F/Y/I/E/P/K/S verified directly
  against the fetched Table 2): W −0.884 … P 0.987.
- **Hand-computed worked normalized scores** confirm spec M8b and the test assertions: W → 0.0;
  I → 0.398/1.871 = **0.21272** (< 0.542, ordered); E → 1.620/1.871 = **0.86585** (≥ 0.542, disordered);
  P → 1.871/1.871 = **1.0** (disordered).
- **Edge-case semantics** (all defined/sourced): empty → zeroed result (INV-7); unknown residues
  skipped in the mean (poly-X → 0.0); window truncated at termini; single residue → window of itself;
  case via `ToUpperInvariant`; scores inherently ∈ [0,1].
- **Findings: None.** The cosmetic ranking-string fix from the prior validation (S before K:
  `…Q,S,K,E,P`) is already present in both `DisorderPredictor.cs:84` and the Evidence artifact — no new
  divergence.

## Stage B — implementation

Code path realises **exactly** the validated normalized-TOP-IDP windowed mean with the published 0.542
cutoff — not an approximation. Cross-verification was **run, not merely traced**: all
`DisorderPredictor` tests green — `DisorderPredictor_DisorderPrediction_Tests` 22/22; the full family
(incl. propensity/classification, MoRF, LowComplexity, region) 113/113. Tests assert exact
externally-confirmed values: 20 propensity values (M8), normalized W=0 / I=0.2127 / E=0.8660 / P=1.0
(M8b/S1), poly-I content 0.0, poly-E/P content 1.0 (M4/M5/M6), hydropathy I=4.5 / W=−0.9 / E=−3.5 (C4).

- **Variant/delegate consistency:** `GetDisorderPropensity`, `IsDisorderPromoting` and the three
  classification-set properties all read the same constant tables; disjoint / cover-20 verified (the
  DISORDER-PROPENSITY-001 file).
- **Numerical robustness:** division guarded by `count > 0`; empty handled; scores bounded [0,1] by
  construction.
- **Test-quality audit:** assertions check exact sourced numbers with tight tolerances, deterministic,
  cover the Stage-A edge cases — not tautologies.
- **Findings: None.**

## Findings & follow-ups

- **No algorithm defect and no code defect (State CLEAN).** No code changed this session.
- All 20 TOP-IDP values, cutoff 0.542, window 21, the normalization, and the Dunker 8/8/4
  classification were independently re-confirmed against the fetched Campen 2008 PMC text and
  Wikipedia/Dunker 2001. Full `DisorderPredictor` test set 113/113 green. **No follow-ups.**
