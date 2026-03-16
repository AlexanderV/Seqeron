# Test Specification: DISORDER-PRED-001

**Test Unit ID:** DISORDER-PRED-001
**Area:** ProteinPred
**Algorithm:** Disorder Prediction
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-03-17

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | URL | Accessed |
|---|--------|---------------|-----|----------|
| 1 | Campen et al. (2008) TOP-IDP Scale | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC2676888/ | 2026-02-11 |
| 2 | Wikipedia — Intrinsically disordered proteins | 3 | https://en.wikipedia.org/wiki/Intrinsically_disordered_proteins | 2026-02-10 |
| 3 | Wikipedia — Hydrophilicity plot (Kyte-Doolittle) | 3 | https://en.wikipedia.org/wiki/Hydrophilicity_plot | 2026-02-10 |
| 4 | Wikipedia — Amino acid | 3 | https://en.wikipedia.org/wiki/Amino_acid | 2026-02-10 |
| 5 | Kyte & Doolittle (1982) | 1 | doi:10.1016/0022-2836(82)90515-0 | — |
| 6 | Dunker et al. (2001) | 1 | PMID 11381529 | — |
| 7 | Uversky et al. (2000) | 1 | — | — |

### 1.2 Key Evidence Points

1. Kyte-Doolittle hydropathy scale: 20 amino acid values from -4.5 (R) to 4.5 (I) — Kyte & Doolittle (1982), Wikipedia
2. Disorder-promoting amino acids: {A, R, G, Q, S, P, E, K} (8 AA) — Dunker et al. (2001)
3. Order-promoting amino acids: {W, C, F, I, Y, V, L, N} (8 AA) — Dunker et al. (2001)
4. Ambiguous amino acids: {D, H, M, T} (4 AA) — Dunker et al. (2001)
5. TOP-IDP scale: 20 amino acid disorder propensity values — Campen et al. (2008), Table 2, PMC2676888
6. IDPs distinguished by low mean hydropathy + high mean net charge — Uversky et al. (2000)
7. Charge at pH 7: R=+1, K=+1, H≈+0.1, D=-1, E=-1, others=0 — Wikipedia Amino acid
8. Proline and Glycine promote conformational flexibility and disorder — Wikipedia Amino acid

### 1.3 Documented Corner Cases

1. Sequences shorter than window size have boundary effects — Kyte & Doolittle (1982) window method
2. Unknown amino acids not in standard 20 — no authoritative guidance; implementation ignores them

### 1.4 Known Failure Modes / Pitfalls

1. Window boundary effects at sequence termini reduce prediction accuracy — inherent to sliding window methods

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `PredictDisorder(sequence, windowSize, threshold)` | DisorderPredictor | **Canonical** | Main prediction entry point |
| `CalculateDisorderScore(window)` | DisorderPredictor | **Internal** | Private; tested indirectly via PredictDisorder |
| `GetDisorderPropensity(char)` | DisorderPredictor | **Canonical** | Returns TOP-IDP propensity — Campen et al. (2008) |
| `IsDisorderPromoting(char)` | DisorderPredictor | **Canonical** | Dunker (2001) classification: true for {A,R,G,Q,S,P,E,K} |
| `DisorderPromotingAminoAcids` | DisorderPredictor | **Canonical** | 8 AA — Dunker et al. (2001) |
| `OrderPromotingAminoAcids` | DisorderPredictor | **Canonical** | 8 AA — Dunker et al. (2001) |
| `AmbiguousAminoAcids` | DisorderPredictor | **Canonical** | 4 AA {D,H,M,T} — Dunker et al. (2001) |
| `CalculateHydropathy(string)` | DisorderPredictor | **Canonical** | Mean Kyte-Doolittle hydropathy |

**Note:** All previously missing methods have been implemented. Scoring uses pure TOP-IDP averaging with published cutoff 0.542.

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | ResiduePredictions.Count == sequence.Length | Yes | Mathematical identity |
| INV-2 | All DisorderScore values ∈ [0, 1] | Yes | Normalized TOP-IDP averaging (inherently [0,1]) |
| INV-3 | OverallDisorderContent ∈ [0, 1] | Yes | Fraction definition |
| INV-4 | MeanDisorderScore ∈ [0, 1] | Yes | Average of normalized values ∈ [0, 1] |
| INV-5 | Case invariance: upper == lower results | Yes | Implementation ToUpperInvariant |
| INV-6 | DisorderPromoting ∪ OrderPromoting ∪ Ambiguous = all 20 AA | Yes | Dunker et al. (2001): 8+8+4=20 |
| INV-7 | Empty sequence → empty result (no predictions, no regions) | Yes | Standard edge case |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | PredictDisorder_EmptySequence_ReturnsEmptyResult | Empty string → zeroed result | Sequence="", no predictions, no regions, content=0, mean=0 | Standard edge case |
| M2 | PredictDisorder_ReturnsCorrectLength | Prediction count matches input length | ResiduePredictions.Count == 20 for "ACDEFGHIKLMNPQRSTVWY" | INV-1 |
| M3 | PredictDisorder_AllScoresInRange | All scores in [0,1] | Every ResiduePrediction.DisorderScore ∈ [0, 1] | INV-2, normalized TOP-IDP averaging |
| M4 | PredictDisorder_HydrophobicSequence_LowDisorder | Poly-Ile (hydropathy 4.5, TOP-IDP -0.486) → all ordered | OverallDisorderContent == 0.0 | Uversky (2000), Campen (2008): normalized 0.2127 < 0.542 |
| M5 | PredictDisorder_ChargedPolarSequence_HighDisorder | Poly-Glu (TOP-IDP 0.736, charge -1) → all disordered | OverallDisorderContent == 1.0 | Uversky (2000), Campen (2008): normalized 0.8660 ≥ 0.542 |
| M6 | PredictDisorder_ProlineRich_HighDisorder | Poly-Pro (highest TOP-IDP 0.987) → all disordered | OverallDisorderContent == 1.0, DisorderedRegions not empty | Campen et al. (2008): normalized 1.0 ≥ 0.542 |
| M7 | PredictDisorder_CaseInsensitive | Upper vs lower case produce same results | MeanDisorderScore equal within tolerance | INV-5 |
| M8 | GetDisorderPropensity_AllTwentyAminoAcids_MatchScale | All 20 values match TOP-IDP scale | Exact values verified | Campen et al. (2008) Table 2 |
| M8b | PredictDisorder_HomopolymericSequences_MatchExactNormalizedTopIdpScores | Exact normalized TOP-IDP score for poly-W/P/E/I | W→0.0, P→1.0, E→0.8660, I→0.2127 | Campen et al. (2008): S=(TOP-IDP+0.884)/1.871 |
| M9 | IsDisorderPromoting_DisorderPromotingResidues_ReturnsTrue | A, R, Q, E, G, K, P, S → true | True for each | Dunker et al. (2001) |
| M10 | IsDisorderPromoting_OrderPromotingResidues_ReturnsFalse | W, C, F, I, Y, V, L, N → false | False for each | Dunker et al. (2001) |
| M10b | IsDisorderPromoting_AmbiguousResidues_ReturnsFalse | D, H, M, T → false (ambiguous ≠ disorder-promoting) | False for each | Dunker et al. (2001) |
| M11 | DisorderPromotingAminoAcids_ContainsExpected | Contains 8 Dunker disorder-promoting AA | Contains A,E,G,K,P,Q,R,S | Dunker et al. (2001) |
| M12 | OrderPromotingAminoAcids_ContainsExpected | Contains 8 Dunker order-promoting AA | Contains C,F,I,L,N,V,W,Y | Dunker et al. (2001) |
| M13 | PredictDisorder_ResiduePredictionsHaveCorrectPositions | Position[i] == i and Residue[i] == sequence[i] | Verified for all positions | Mathematical identity |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | PredictDisorder_SingleResidue_Handles | Single Pro input → exact score | 1 prediction, score == 1.0 (normalized TOP-IDP for P) | Boundary condition, Campen (2008) |
| S2 | PredictDisorder_ShortSequence_Handles | Sequence shorter than window | Correct count, no crash | Boundary condition |
| S3 | PredictDisorder_UnknownResidue_Handles | "XXXXX" → all scores 0.0 | 5 predictions, all scores == 0.0 | Unknown AA tolerance |
| S4 | PredictDisorder_MeanDisorderScore_IsAverage | Mean equals average of all residue scores | Computed manually | INV-4 |
| S5 | PredictDisorder_OverallDisorderContent_IsFraction | Content = disordered count / length | Verified against predictions | INV-3 |
| S6 | GetDisorderPropensity_UnknownResidue_ReturnsZero | Unknown char returns 0.0 | 0.0 | Implementation GetValueOrDefault |
| S7 | PredictDisorder_LowercaseInput_HandledCorrectly | "ppppeeee" same as "PPPPEEEE" | Same scores | Case handling |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | PredictDisorder_MixedOrderedDisordered_FindsTransition | Ordered flanks + disordered middle | MeanDisorderScore > 0, at least one region detected | Functional verification |
| C2 | PredictDisorder_CustomWindowSize_Respected | Window size parameter works | Same count, different per-residue scores on mixed P/I | Parameter verification |
| C3 | AmbiguousAminoAcids_ContainsExpected | Contains {D, H, M, T} | 4 AA — Dunker (2001) | Classification completeness |
| C4 | CalculateHydropathy_ReturnsCorrectValues | Mean Kyte-Doolittle hydropathy | Verified against known values | Kyte & Doolittle (1982) |
| C5 | ClassificationSets_AreDisjointAndCoverAll20 | 8+8+4=20, pairwise disjoint | Verified | Dunker et al. (2001) |

---

## 5. Audit of Existing Tests

### 5.1 Test Files

- **Canonical:** `tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderPrediction_Tests.cs` (49 tests)
- **Future scope:** `tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictorTests.cs` (MoRF/LowComplexity for DISORDER-REGION-001)

### 5.2 Coverage Classification

All test cases are implemented in the canonical file `DisorderPredictor_DisorderPrediction_Tests.cs`.

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1: Empty sequence | ✅ Covered | `PredictDisorder_EmptySequence_ReturnsEmptyResult` |
| M2: Correct length | ✅ Covered | `PredictDisorder_AllTwentyAminoAcids_ReturnsCorrectLength` |
| M3: Scores in range | ✅ Covered | `PredictDisorder_AllScoresInZeroOneRange` — diverse sequences |
| M4: Hydrophobic → low disorder | ✅ Covered | Poly-Ile: OverallDisorderContent == 0.0 (all scores 0.2127 < 0.542) |
| M5: Charged → high disorder | ✅ Covered | Poly-Glu: OverallDisorderContent == 1.0 (all scores 0.8660 ≥ 0.542) |
| M6: Proline-rich → all disordered | ✅ Covered | OverallDisorderContent == 1.0 + DisorderedRegions not empty |
| M7: Case insensitive | ✅ Covered | `PredictDisorder_CaseInsensitive_SameResults` |
| M8: All 20 propensity values | ✅ Covered | `GetDisorderPropensity_AllTwentyAminoAcids_MatchScale` — all 20 exact values |
| M8b: Exact normalized scores | ✅ Covered | Poly-W→0.0, Poly-P→1.0, Poly-E→0.8660, Poly-I→0.2127 |
| M9: Disorder-promoting classification | ✅ Covered | All 8 AA tested: A,R,Q,E,G,K,P,S |
| M10: Order-promoting classification | ✅ Covered | All 8 AA tested: C,F,I,L,N,V,W,Y |
| M10b: Ambiguous classification | ✅ Covered | All 4 AA tested: D,H,M,T |
| M11: DisorderPromotingAminoAcids | ✅ Covered | 8 AA verified + count check |
| M12: OrderPromotingAminoAcids | ✅ Covered | 8 AA verified + count check |
| M13: Position/Residue correctness | ✅ Covered | `PredictDisorder_ResiduePredictionsHaveCorrectPositionsAndResidues` |
| S1: Single residue | ✅ Covered | Exact score 1.0 for single Pro (normalized TOP-IDP) |
| S2: Short sequence | ✅ Covered | `PredictDisorder_ShortSequence_HandledCorrectly` |
| S3: Unknown residue | ✅ Covered | All scores == 0.0 for "XXXXX" (unknown contributes nothing) |
| S4: Mean is average of scores | ✅ Covered | `PredictDisorder_MeanDisorderScore_IsAverageOfResidueScores` |
| S5: Content is fraction | ✅ Covered | `PredictDisorder_OverallDisorderContent_IsFraction` |
| S6: Unknown propensity zero | ✅ Covered | `GetDisorderPropensity_UnknownResidue_ReturnsZero` |
| S7: Lowercase input | ✅ Covered | `GetDisorderPropensity_LowercaseInput_SameAsUppercase` |
| C1-C5: COULD tests | ✅ Covered | All 5 COULD tests implemented |

### 5.3 Canonical Test Files

| File | Role | Test Count |
|------|------|------------|
| `DisorderPredictor_DisorderPrediction_Tests.cs` | DISORDER-PRED-001 canonical | 49 |
| `DisorderPredictorTests.cs` | DISORDER-REGION-001 (future) | ~12 |

---

## 6. Assumption Register

**Total assumptions:** 0

All parameters sourced from published peer-reviewed research. No implementation-specific assumptions remain.

---

## 7. Notes

Three classification sets (disorder/order/ambiguous) cover all 20 AA exactly — verified in C5.
