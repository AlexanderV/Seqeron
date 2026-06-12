# Validation Report: DISORDER-PRED-001 — Protein Intrinsic Disorder Prediction

- **Validated:** 2026-06-12   **Area:** ProteinPred
- **Canonical method(s):** `DisorderPredictor.PredictDisorder`, `GetDisorderPropensity`, `IsDisorderPromoting`, `CalculateHydropathy`, `{Disorder,Order,Ambiguous}PromotingAminoAcids`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **State:** CLEAN

## Method scope (honest scoping)

This is a **composition/propensity heuristic**, not an energy-based (IUPred) or ML
predictor. Per-residue disorder = average of the **normalized TOP-IDP scale**
(Campen et al. 2008) over a sliding window (default 21), thresholded at the
published normalized cutoff **0.542**. The source's own XML doc comments state this
explicitly (single-feature toolkit, AUC ≈ 0.65–0.72 vs IUPred2A ≈ 0.75–0.80,
flDPnn ≈ 0.85–0.90) and recommend dedicated tools for publication-grade work. The
advertised accuracy is honestly scoped; no false IUPred-grade claim.

## Stage A — Description

### Sources opened
- **Campen et al. (2008), TOP-IDP, PMC2676888** (fetched). Confirms: (a) all 20
  Table-2 scale values; (b) the scale is **normalized to [0,1]** for prediction;
  (c) the prediction index `I = -(<TOP-IDP> - 0.542)` on the **normalized averaged**
  score (so `<TOP-IDP> > 0.542` ⟹ disordered); (d) window size 21; (e) window value
  is plotted at the center residue.
- **Wikipedia "Intrinsically disordered proteins"** + **Dunker (2001) PMID 11381529**
  (search). Confirms classification: order-promoting {W,C,F,I,Y,V,L,N},
  disorder-promoting {A,R,G,Q,S,P,E,K}, neutral/ambiguous {D,H,M,T} — 8+8+4=20,
  disjoint.
- **Kyte & Doolittle (1982)** hydropathy scale (Evidence doc) — all 20 values match.
- **Uversky et al. (2000)** charge-hydropathy rationale — low hydropathy + high net
  charge ⟹ disorder; consistent with poly-I ordered / poly-E,P disordered.

### Formula check
- Per-residue normalized propensity `S(aa) = (TOP-IDP(aa) + 0.884) / 1.871`
  (min = W = −0.884, max = P = 0.987, range 1.871). Matches Campen "normalized to
  min 0, max 1".
- Window score = mean of `S(aa)` over the truncated centered window. Threshold
  `score >= 0.542` ⟹ disordered. Matches Campen's normalized-cutoff convention.

### TOP-IDP values (all 20) — external (Campen Table 2) vs code
W −0.884, F −0.697, Y −0.510, I −0.486, M −0.397, L −0.326, V −0.121, N 0.007,
C 0.020, T 0.059, A 0.060, G 0.166, R 0.180, D 0.192, H 0.303, Q 0.318, S 0.341,
K 0.586, E 0.736, P 0.987 — **all match the implementation exactly.**

### Hand-computed worked examples (normalized score)
- W: (−0.884+0.884)/1.871 = **0.0**
- I: (−0.486+0.884)/1.871 = 0.398/1.871 = **0.2127** (< 0.542 → ordered)
- E: (0.736+0.884)/1.871 = 1.620/1.871 = **0.8660** (≥ 0.542 → disordered)
- P: (0.987+0.884)/1.871 = 1.871/1.871 = **1.0** (≥ 0.542 → disordered)
All confirm the spec's M8b expected values.

### Findings / divergences (Stage A)
- **Cosmetic:** the ranking string in the Evidence doc (line 188) and the source
  comment listed `…Q,K,S,E,P`, but since S=0.341 < K=0.586 the correct order is
  `…Q,S,K,E,P`. Values were always correct; only the human-readable ordering string
  was wrong. Source comment fixed (see Files changed). Not a computational defect.

## Stage B — Implementation

- **Code path:** `DisorderPredictor.cs:180` `PredictDisorder` →
  `CalculatePerResidueScores` (245) → `CalculateDisorderScore` (245) normalizes &
  averages; threshold at `:232`.
- **Formula realised:** normalization `(prop − TopIdpMin)/TopIdpRange` with
  TopIdpMin=−0.884, TopIdpRange=1.871, cutoff 0.542 — exactly the validated method.
- **Edge cases:** empty → zeroed result (`:186`, INV-7); unknown residues skipped in
  averaging (poly-X → 0.0, S3/S6); window truncated at termini (S2); single residue
  → window of itself (S1, poly-P → 1.0); case via `ToUpperInvariant` (M7/S7);
  scores inherently ∈[0,1] (INV-2/M3). `>=` boundary convention is the safe choice.
- **Cross-verification:** tests assert exact externally-confirmed values — all 20
  propensity values (M8), normalized worked examples W=0, I=0.2127, E=0.8660, P=1.0
  (M8b), poly-I content 0.0, poly-E/P content 1.0. 124 disorder tests pass.
- **Variant consistency:** `GetDisorderPropensity`, `IsDisorderPromoting`,
  classification set properties all read the same constant tables; disjoint/cover-20
  verified (C5).
- **Test quality:** assertions check exact sourced numbers with tight tolerances,
  deterministic, cover all Stage-A edge cases. Real tests, not tautologies.

### Findings / defects (Stage B)
None affecting computation. Only the cosmetic comment typo (fixed).

## Verdict & follow-ups
- Stage A PASS, Stage B PASS, **CLEAN**.
- Fix applied: corrected ranking-order comment in `DisorderPredictor.cs:82`
  (`Q,K,S,E,P` → `Q,S,K,E,P`); values unchanged. Build + full suite green
  (4486 passed). The same stale ordering remains in
  `docs/Evidence/DISORDER-PRED-001-Evidence.md:188` (documentation only, non-load-bearing).
