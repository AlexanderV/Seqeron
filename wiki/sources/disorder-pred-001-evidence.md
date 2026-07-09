---
type: source
title: "Evidence: DISORDER-PRED-001 (intrinsic-disorder prediction — TOP-IDP sliding window)"
tags: [validation, analysis]
doc_path: docs/Evidence/DISORDER-PRED-001-Evidence.md
sources:
  - docs/Evidence/DISORDER-PRED-001-Evidence.md
source_commit: 05fff695e889b79023301d7319afbc8a24e0bec4
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: DISORDER-PRED-001

The validation-evidence artifact for test unit **DISORDER-PRED-001** — **intrinsic-disorder
prediction**, the `PredictDisorder` sliding-window average of the TOP-IDP amino-acid propensity
scale. This is the **third ingested unit of the protein disorder / features family** (after
[[disorder-lc-001-evidence|DISORDER-LC-001]] and [[disorder-morf-001-evidence|DISORDER-MORF-001]])
and one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern. The algorithm itself — the **shared `PredictDisorder` anchor** on which
[[morf-prediction-dip-in-disorder|MoRF prediction]] and disordered-region detection sit — is written
up on the concept page [[intrinsic-disorder-prediction-top-idp]]; this file records the source trace
and worked oracles. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources (four, with authority ranks):**
  - **Campen et al. 2008** "TOP-IDP-Scale: A New Amino Acid Scale Measuring Propensity for Intrinsic
    Disorder" (*Protein Pept Lett* 15(9):956–963, PMC2676888, PMID 18991772, **rank 1 / primary**) —
    the **exact 20-value TOP-IDP scale** (Table 2, order→disorder: W −0.884, F −0.697, Y −0.510,
    I −0.486, M −0.397, L −0.326, V −0.121, N 0.007, C 0.020, T 0.059, A 0.060, G 0.166, R 0.180,
    D 0.192, H 0.303, Q 0.318, S 0.341, K 0.586, E 0.736, P 0.987), the scale derivation (517 scales
    surveyed, simulated annealing, ARV 0.761 = 11% over the best prior scale), the **prediction
    cutoff 0.542** (maximum-likelihood), and the **21-residue** evaluation window.
  - **Wikipedia** "Intrinsically disordered proteins" (rank 3) — IDPs lack fixed 3-D structure and
    exist as dynamic ensembles; disorder-promoting AAs are the polar/charged Q,S,P,E,K (± R,G) and A;
    order-promoting are the bulky hydrophobics V,L,I,M,F,W,Y plus C,N; Uversky 2000 charge–hydropathy
    (CH) plot; MoRFs (Mohan 2006); ~33% of eukaryotic proteins have long (>30 aa) disordered regions.
  - **Wikipedia** "Hydrophilicity plot (Kyte–Doolittle)" (rank 3) — the full 20-value Kyte &
    Doolittle 1982 hydropathy scale (I 4.5 … R −4.5) and the sliding-window averaging method, used by
    the `CalculateHydropathy` utility.
  - **Wikipedia** "Amino acid (properties table)" (rank 3) — cross-validates the Kyte–Doolittle
    values; charge at pH 7 (R/K +1, H ~+0.1, D/E −1); why Pro (rigid cyclic side chain, α-helix
    breaker) and Gly (no side chain, maximal flexibility) promote disorder.
- **Dunker et al. 2001 classification (confirmed by Campen):** disorder-promoting = {A, R, G, Q, S,
  P, E, K}; order-promoting = {W, C, F, I, Y, V, L, N}; ambiguous/borderline = {D, H, M, T}. The
  three sets are disjoint and cover all 20 residues.
- **Datasets:** the TOP-IDP propensity table (all 20 values + per-AA class), the Kyte–Doolittle
  hydropathy scale (all 20), charge-at-pH-7, and the Dunker disorder/order/ambiguous partition —
  each with its primary source.
- **Documented corner cases / failure modes:** input **shorter than the window** → boundary effects
  dominate; scores at **ordered↔disordered boundaries** unreliable (window averaging); **non-standard
  residues** (X/B/Z) absent from the scales; single-residue input (window 1) has no averaging
  context; homopolymers give uniform per-residue scores with no window benefit.
- **Recommended coverage (17 items):** MUST — all-20 TOP-IDP values match Campen Table 2 (M8);
  disorder-promoting {A,R,G,Q,S,P,E,K} (M9) and order-promoting {W,C,F,I,Y,V,L,N} (M10) sets;
  ambiguous {D,H,M,T} NOT disorder-promoting (M10b); hydrophobic poly-I low score (M4);
  charged/polar poly-E / poly-P high score (M5/M6); empty→zero-init (M1); residue-prediction count =
  length (M2); all scores ∈ [0,1] (M3); case-insensitivity (M7); correct residue positions (M13).
  SHOULD — the three public residue-set properties, `CalculateHydropathy` = mean Kyte–Doolittle
  (C4), the ambiguous set {D,H,M,T} (C3), three sets disjoint & complete (C5). COULD — O(n)
  performance claim.

## Deviations and assumptions

**Assumptions: None** in the evidence artifact itself — every parameter is traced to peer-reviewed
sources: disorder propensity + cutoff 0.542 from Campen et al. 2008 (Table 2 / maximum-likelihood),
disorder/order/ambiguous sets from Dunker et al. 2001, hydropathy from Kyte & Doolittle 1982. No
source contradictions.

**Implementation-side note (from the algorithm doc, not the evidence file):** the repository
`DisorderPredictor` is an explicitly **simplified, single-feature TOP-IDP heuristic** — it adds no
evolutionary profiles, predicted structure, or trained-model features and is not competitive with
IUPred2A / MobiDB-lite; non-canonical residues are skipped in the window average (a window with no
recognized residues scores `0.0`); edge windows are **clipped** to the sequence bounds rather than
padded. These are captured on [[intrinsic-disorder-prediction-top-idp]]. The evidence file uses the
`0.542` published cutoff; the sibling [[morf-prediction-dip-in-disorder|MoRF]] unit uses the `0.5`
order/disorder threshold (PMC2570644) — different published thresholds for different purposes, not a
contradiction.
</content>
</invoke>
