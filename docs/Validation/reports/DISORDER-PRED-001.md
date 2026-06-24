# Validation Report: DISORDER-PRED-001 — Protein Intrinsic Disorder Prediction

- **Validated:** 2026-06-24   **Area:** ProteinPred
- **Canonical method(s):** `DisorderPredictor.PredictDisorder` (+ internal `CalculateDisorderScore`), `CalculateHydropathy`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **State:** CLEAN

## Method scope (honest scoping)

Per-residue disorder is a **composition heuristic**, not an energy-based (IUPred) or
ML predictor. Score = mean of the **normalized TOP-IDP scale** (Campen et al. 2008,
PMC2676888) over a sliding window (default 21, centered, truncated at termini),
thresholded at the published cutoff **0.542** (`score >= 0.542` ⟹ disordered). The
source XML doc comments scope this honestly (single-feature, AUC ≈ 0.65–0.72 vs
IUPred2A 0.75–0.80, flDPnn 0.85–0.90) and recommend dedicated tools for
publication-grade work. No false IUPred-grade claim.

## Stage A — Description

### Sources opened (this session, via web)
- **Campen et al. (2008) TOP-IDP, PMC2676888** (WebFetch of the PMC full text).
  Table 2 values confirmed verbatim; cutoff = 0.542; prediction equation
  `I = -(<TOP-IDP> - 0.542)` (positive ⟹ ordered, so `<TOP-IDP> > 0.542` ⟹
  disordered); scales normalized to min 0 / max 1; window 21.
- **Wikipedia "Intrinsically disordered proteins" + Dunker (2001) PMID 11381529**
  (WebSearch). Confirms classification: order-promoting {W,C,F,I,Y,V,L,N},
  disorder-promoting {A,R,G,Q,S,P,E,K}, ambiguous/neutral {H,M,T,D} — 8+8+4=20,
  disjoint. Matches code exactly.

### Formula check
- Per-residue normalized propensity `S(aa) = (TOP-IDP(aa) − TopIdpMin)/TopIdpRange`
  with `TopIdpMin = −0.884 (W)`, `TopIdpMax = 0.987 (P)`, `TopIdpRange = 1.871`.
  Matches Campen "normalized to min 0, max 1".
- Window score = mean of `S(aa)` over the truncated centered window; `score >= 0.542`
  ⟹ disordered. Direction matches Campen's `I = -(<TOP-IDP>-0.542)` (the code applies
  0.542 to the **normalized** averaged score, consistent with the all-scales-normalized
  convention used for the published cutoff).

### TOP-IDP values — external (Campen Table 2, fetched) vs code constant table
W −0.884, F −0.697, Y −0.510, I −0.486, M −0.397, L −0.326, V −0.121, N 0.007,
C 0.020, T 0.059, A 0.060, G 0.166, R 0.180, D 0.192, H 0.303, Q 0.318, S 0.341,
K 0.586, E 0.736, P 0.987 — **all 20 match the implementation exactly** (key anchors
W/F/Y/I/E/P/K/S verified directly against the fetched Table 2).

### Hand-computed worked examples (normalized score)
- W: (−0.884+0.884)/1.871 = **0.0**
- I: (−0.486+0.884)/1.871 = 0.398/1.871 = **0.21272** (< 0.542 → ordered)
- E: (0.736+0.884)/1.871 = 1.620/1.871 = **0.86585** (≥ 0.542 → disordered)
- P: (0.987+0.884)/1.871 = 1.871/1.871 = **1.0** (≥ 0.542 → disordered)
All confirm spec M8b and the test assertions.

### Edge-case semantics
Empty → zeroed result (INV-7); unknown residues skipped in the mean (poly-X → 0.0);
window truncated at termini; single residue → window of itself; case via
ToUpperInvariant; scores inherently ∈[0,1]. All defined and sourced/standard.

### Findings / divergences (Stage A)
None. The cosmetic ranking-string fix from the prior validation (S before K:
`…Q,S,K,E,P`) is already present in both `DisorderPredictor.cs:84` and
`docs/Evidence/DISORDER-PRED-001-Evidence.md:188`. No new divergence.

## Stage B — Implementation

- **Code path:** `DisorderPredictor.cs` — `PredictDisorder` (:190) →
  `CalculatePerResidueScores` (:227, centered truncated window) →
  `CalculateDisorderScore` (:255, normalize `(prop−TopIdpMin)/TopIdpRange` and
  average); threshold `score >= disorderThreshold` (:242). Constants:
  `TopIdpMin=−0.884`, `TopIdpMax=0.987`, `TopIdpRange=1.871`, `TopIdpCutoff=0.542`.
- **Formula realised:** exactly the validated normalized-TOP-IDP windowed mean with
  published cutoff — not an approximation.
- **Cross-verification (run, not traced):** all DisorderPredictor tests green —
  `DisorderPredictor_DisorderPrediction_Tests` 22/22; the full DisorderPredictor
  family (incl. propensity/classification, MoRF, LowComplexity, region) 113/113.
  Tests assert exact externally-confirmed values: 20 propensity values (M8),
  normalized W=0/I=0.2127/E=0.8660/P=1.0 (M8b/S1), poly-I content 0.0, poly-E/P
  content 1.0 (M4/M5/M6), hydropathy I=4.5/W=−0.9/E=−3.5 (C4).
- **Variant/delegate consistency:** `GetDisorderPropensity`, `IsDisorderPromoting`,
  and the three classification-set properties all read the same constant tables;
  disjoint/cover-20 verified (DISORDER-PROPENSITY-001 file).
- **Numerical robustness:** division guarded by `count > 0`; empty handled; scores
  bounded [0,1] by construction.
- **Test quality:** assertions check exact sourced numbers with tight tolerances,
  deterministic, cover the Stage-A edge cases. Not tautologies.

### Findings / defects (Stage B)
None.

## Verdict & follow-ups
- Stage A PASS, Stage B PASS, **CLEAN**. No code changed this session.
- All 20 TOP-IDP values, cutoff 0.542, window 21, normalization, and the Dunker
  8/8/4 classification independently re-confirmed against the fetched Campen 2008
  PMC text and Wikipedia/Dunker 2001. Full DisorderPredictor test set 113/113 green.
