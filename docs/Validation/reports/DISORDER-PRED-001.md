# Validation Report: DISORDER-PRED-001 — Protein Intrinsic Disorder Prediction (TOP-IDP)

- **Validated:** 2026-09-28 (review-2026-09, batch B15; supersedes 2026-06-24)   **Area:** ProteinPred
- **Canonical method(s):** `DisorderPredictor.PredictDisorder` (+ `PredictDisorderRegions`, private `CalculatePerResidueScores` / `CalculateDisorderScore`), `DisorderPredictor.CalculateHydropathy`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES
- **State:** CLEAN

## Stage A — Description

### Sources opened (this session)
- **Campen et al. (2008) PMC2676888** — direct fetch blocked by egress proxy; content obtained via
  WebSearch snippets of the PMC page: "the average global TOP-IDP value and average window-by-window
  TOP-IDP values are calculated based on the normalized TOP-IDP scale"; "normalized to have the minimal
  value of zero and the maximal value of 1"; index I = −(⟨TOP-IDP⟩ − 0.542), positive ⟹ ordered,
  negative ⟹ disordered; cut-off 0.542 from maximum likelihood. (Table 2 values were fetched verbatim
  in the 2026-06 session; not re-fetchable now.)
- **localCIDER 0.1.21** (PyPI sdist, `localcider/backend/sequence.py::fraction_disorder_promoting`,
  citing Campen 2008): order list W,F,Y,I,M,L,V,N,C / disorder list T,A,G,R,D,H,Q,K,S,E,P — the same
  TOP-IDP ranking split as the code's scale.
- **Biopython 1.88** `Bio.SeqUtils.ProtParamData.kd` and `ProteinAnalysis.gravy()` — Kyte-Doolittle table
  and GRAVY definition (identical to the 20 values in code).

### Formula check
S(aa) = (TOP-IDP(aa) + 0.884)/1.871 ∈ [0,1]; residue score = mean S over the centered window (21,
truncated at termini); disordered iff score ≥ 0.542. Matches the paper (cutoff applied to the averaged
normalized scale). Hydropathy = mean KD over standard residues (GRAVY).

### Notes (divergences, documented, not defects)
1. Tie at exactly 0.542: Campen's I = 0 is unclassified; code resolves `>=` toward disorder (measure-zero).
2. Terminal windows are truncated (paper silent on termini handling).
3. Even `windowSize` w is centered with w/2 each side → spans w+1 residues (paper uses odd 21). Now documented.

## Stage B — Implementation

- **Code path:** `DisorderPredictor.cs` `PredictDisorder` → `ComputeDisorder` → `CalculatePerResidueScores`
  → `CalculateDisorderScore`; `CalculateHydropathy` → `SequenceStatistics.CalculateHydrophobicity`.
- **Cross-verification (independent Python re-implementation of the Campen formula, scratch `topidp_ref.py`):**

| Input | Quantity | Reference | Code |
|---|---|---|---|
| α-synuclein P37840, defaults | score[0] / [10] / [69] / [70] / [120] / [139] | 0.4721345 / 0.5422361 / 0.5406836 / 0.5272709 / 0.6138301 / 0.6496769 | identical (1e-12) |
| same | content / mean | 102/140 = 0.7285714 / 0.5766837 | identical |
| same | regions (0-based incl., ≥5) | [10–43], [47–66], [94–139] | identical |
| I10P10E5, w=7 | 25 window means | 0.21272 … 0.865847 | identical (1e-6) |
| α-syn / 20-AA / MKWVTFISLLLLFSSAYS | GRAVY (Biopython) | −0.4028571 / −0.49 / 1.2888889 | identical |

- **Edge cases:** empty → zeroed result; unknown residues skipped (poly-X → 0); case-insensitive.
  **Defect (robustness):** `windowSize = 0` silently acted as a 1-residue window and negative values
  crashed with an opaque range exception inside `string[start..end]` (e.g. `PredictDisorder("PPPPPEEEEE", -1)`).
  Fixed: `ArgumentOutOfRangeException` for `windowSize < 1` in `PredictDisorder` and `PredictDisorderRegions`.
- **Duplication:** `CalculateHydropathy` re-implemented the KD table + GRAVY mean already canonical in
  `SequenceStatistics.CalculateHydrophobicity` (same project, same semantics). Now delegates; private
  `Hydropathy` table removed. Behaviour-preserving (verified by existing C4 tests + new Biopython tests).
- **MCP:** `AnalysisTools.PredictDisorder` (Seqeron.Mcp.Analysis) delegates to `DisorderPredictor.PredictDisorder`.
- **Tests added** (`DisorderPredictor_DisorderPrediction_Tests.cs`): `PredictDisorder_AlphaSynuclein_MatchesIndependentTopIdpReference`,
  `PredictDisorder_Window7_OrderDisorderTransition_MatchesReference`, `PredictDisorder_NonPositiveWindow_Throws` (×3),
  `CalculateHydropathy_MatchesBiopythonGravy`, `CalculateHydropathy_NonStandardResiduesSkipped_DelegatesToCanonicalGravy`.

## Verdict & follow-ups
- Stage A PASS-WITH-NOTES, Stage B PASS-WITH-NOTES, **CLEAN**. The TOP-IDP method is the published
  algorithm itself (single-scale predictor), not a simplification of it; the class-level XML doc already
  states its accuracy relative to IUPred2A/flDPnn.
- Out of unit scope (for the B15 duplication sweep): `ProteinSequence.Gravy()` (Core, B02) and
  `ProteinMotifFinder` KD table (B14) duplicate the KD scale; private `CalculateShannonEntropy` (SEG,
  DISORDER-LC-001) is a candidate for the canonical entropy helper.
