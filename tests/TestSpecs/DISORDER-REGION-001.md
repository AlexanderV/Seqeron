# Test Specification: DISORDER-REGION-001

**Test Unit ID:** DISORDER-REGION-001
**Area:** ProteinPred
**Algorithm:** Disordered Region Detection
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-09-28

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Campen et al. (2008) TOP-IDP Scale | 1 | https://doi.org/10.2174/092986608785849164 | 2026-02-12 |
| 2 | Dunker et al. (2001) Intrinsically disordered protein | 1 | https://doi.org/10.1016/s1093-3263(00)00138-8 | 2026-02-12 |
| 3 | van der Lee et al. (2014) Classification of IDRs | 1 | https://doi.org/10.1021/cr400525m | 2026-02-12 |
| 4 | Ward et al. (2004) Disorder prediction | 1 | https://doi.org/10.1016/j.jmb.2004.02.002 | 2026-02-12 |
| 5 | Wikipedia — Intrinsically disordered proteins | 4 | https://en.wikipedia.org/wiki/Intrinsically_disordered_proteins | 2026-02-12 |

### 1.2 Key Evidence Points

1. IDRs are contiguous segments where per-residue disorder scores exceed a threshold — Campen et al. (2008)
2. TOP-IDP prediction cutoff = 0.542, window = 21 residues — Campen et al. (2008)
3. IDRs can be classified by amino acid composition biases (proline-rich, acidic, basic, Ser/Thr-rich) — van der Lee et al. (2014)
4. Long IDRs (>30 residues) are functionally significant — Ward et al. (2004), Wikipedia citing Ward
5. The implementation scores residues using normalized TOP-IDP values and builds regions from contiguous runs exceeding the threshold

### 1.3 Documented Corner Cases

1. Empty predictions list → no regions
2. All residues ordered → no regions
3. All residues disordered → one region spanning entire sequence
4. Region at end of sequence (trailing region must be captured)
5. Short runs below minLength must be excluded
6. Region exactly at minLength boundary

### 1.4 Known Failure Modes / Pitfalls

1. Off-by-one in trailing region detection — if the region extends to the last residue, the "else" branch is never hit; end-of-loop handling needed
2. Window boundary effects blur order/disorder transitions — Campen et al. (2008)

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `IdentifyDisorderedRegions(predictions, threshold, minLen)` | DisorderPredictor | Canonical (private) | Tested indirectly via `PredictDisorder()` |
| `ClassifyDisorderedRegion(region)` | DisorderPredictor | Canonical (private) | Tested indirectly via `PredictDisorder()` |
| `PredictDisorder(sequence, windowSize, threshold, minRegionLength)` | DisorderPredictor | Public API | Entry point for testing both canonical methods |
| `ClassifyRegionFlavorMobiDbLite(regionSequence)` | DisorderPredictor | Canonical (public, opt-in) | Sourced MobiDB-lite 3.0 disorder-flavor label (Necci et al. 2020); does not affect boundaries |
| `PredictFlavorSubregionsMobiDbLite(sequence, regions?, lcMask?)` | DisorderPredictor | Canonical (public, opt-in) | Verbatim port of MobiDB-lite v3 `get_region_features` (7-residue mirrored window, math_morphology rmax 5, merge, runs ≥ 10 within IDRs) |
| `PredictDisorderRegions(sequence, …)` | DisorderPredictor | Public API | Same boundaries as `PredictDisorder`, `Confidence = NaN`, never blocked by the limitation policy |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | All regions have Start ≥ 0 and End < sequence.Length | Yes | Algorithm definition |
| INV-2 | Region.End ≥ Region.Start for every region | Yes | Algorithm definition |
| INV-3 | Region length (End - Start + 1) ≥ minRegionLength | Yes | Algorithm definition |
| INV-4 | Regions are non-overlapping and sorted by Start | Yes | Single-pass scan |
| INV-5 | MeanScore is in [0, 1] (normalized TOP-IDP range) | Yes | Campen et al. (2008) |
| INV-6 | Confidence is in [0, 1] | Yes | Formula definition |
| INV-7 | RegionType is one of: "Proline-rich", "Acidic", "Basic", "Ser/Thr-rich", "Long IDR", "Standard IDR" | Yes | Classification definition |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | AllOrderedProducesNoRegions | 30×W (lowest TOP-IDP, -0.884) → no disordered regions | 0 regions | Campen (2008) Table 2 |
| M2 | AllDisorderedProducesOneRegion | 30×P (highest TOP-IDP, 0.987) → one region spanning entire sequence | 1 region, Start=0, End=29 | Campen (2008) Table 2 |
| M3 | RegionBoundariesCorrect | Verify Start and End are correct for a known disordered sequence | Start ≥ 0, End < len, End ≥ Start | Algorithm definition |
| M5 | MinLengthFiltering | Disordered region shorter than minLength is excluded | 0 regions | Algorithm definition |
| M6 | TrailingRegionCaptured | Disorder at end of sequence → region includes last residue | Region.End = len - 1 | Algorithm definition |
| M7 | ProlineRichClassification | 30×P → classified as "Proline-rich" | RegionType = "Proline-rich" | Name: van der Lee (2014); threshold/algorithm: internal D1, D2 |
| M8 | AcidicClassification | 30×E → classified as "Acidic" | RegionType = "Acidic" | Name: van der Lee (2014); threshold/groups: internal D1, D4 |
| M9 | BasicClassification | K/R-rich sequence → classified as "Basic" | RegionType = "Basic" | Name: van der Lee (2014); threshold/groups: internal D1, D4 |
| M10 | SerThrRichClassification | 30×S → classified as "Ser/Thr-rich" | RegionType = "Ser/Thr-rich" | Name: van der Lee (2014); threshold/groups: internal D1, D4 |
| M11 | LongIdrClassification | Long disorder-promoting sequence with no dominant AA → "Long IDR" | RegionType = "Long IDR" | Ward (2004) |
| M12 | StandardIdrClassification | Short disorder-promoting sequence with no dominant AA → "Standard IDR" | RegionType = "Standard IDR" | Fallback |
| M13 | ConfidenceInRange | All region confidences in [0, 1] | Confidence ∈ [0, 1] | Cutoff: Campen (2008); formula: internal D3 |
| M14 | RegionsNonOverlapping | Multiple regions do not overlap | No start/end intersection | Algorithm definition |
| F1 | FlavorPolyampholyte | `RKDERKDE`: FCR=1.0>0.35, NCPR=0≤0.35 | `Polyampholyte` | Necci (2020); Das & Pappu (2013); `states.py:get_disorder_class` |
| F2 | FlavorPositivePolyelectrolyte | `RKRKRKRKRR`: f₊=1.0>0.35 | `PositivePolyelectrolyte` | Necci (2020); `get_disorder_class` |
| F3 | FlavorNegativePolyelectrolyte | `DEDEDEDEDD`: f₋=1.0>0.35 | `NegativePolyelectrolyte` | Necci (2020); `get_disorder_class` |
| F4 | FlavorChargeBeatsComposition | `RKRKPPPPPP`: FCR=0.4>0.35, f₊>0.35 (charge tested first) | `PositivePolyelectrolyte` | Necci (2020); `consensus.py` priority |
| F5 | FlavorFcrExactlyThreshold | f₊=0.35 (FCR=0.35, not >0.35) → no comp. | `WeaklyCharged` | Necci (2020); strict `>` in `get_disorder_class` |
| F6 | FlavorCysteineRich | `CCCCAAAAAA`: C=0.4≥0.32 | `CysteineRich` | Necci (2020); `is_enriched(threshold=0.32)` |
| F7 | FlavorProlineRich | `PPPPAAAAAA`: P=0.4≥0.32 | `ProlineRich` | Necci (2020); `is_enriched` |
| F8 | FlavorGlycineRich | `GGGGAAAAAA`: G=0.4≥0.32 | `GlycineRich` | Necci (2020); `is_enriched` |
| F9 | FlavorPolar | `SSTTNNQQAA`: {S,T,N,Q}=0.8≥0.32 | `Polar` | Necci (2020); `is_enriched(['S','T','N','Q'])` |
| F10 | FlavorCompositionPriority | `CCCCPPPPAA`: C and P both 0.4, C first | `CysteineRich` | Necci (2020); `consensus.py` C→P→G→polar |
| F11 | FlavorThresholdInclusive | 8/25 C = 0.32 exactly (`≥`) | `CysteineRich` | Necci (2020); `s >= threshold` |
| F12 | FlavorJustBelowThreshold | 7/25 C = 0.28 < 0.32 | `WeaklyCharged` | Necci (2020); `s >= threshold` |
| F13 | FlavorNoEnrichmentFallback | hydrophobic stretch, FCR=0, none enriched | `WeaklyCharged` | Necci (2020) |
| F15 | FlavorNullOrEmptyThrows | null / "" region sequence | `ArgumentException` | Input validation |
| F16 | FlavorBoundariesUnchanged | 30×P region stays `[0,29]`; flavor `ProlineRich` | Start=0, End=29 | Boundaries from TOP-IDP (Campen 2008), unaffected |
| FS1 | FlavorSubregions two blocks | A10+(RK)7+A11+P13+A10, IDR [0,57], no LC | PPE (9,24), ProlineRich (34,48) | MobiDB-lite v3 code run verbatim |
| FS2 | FlavorSubregions α-synuclein | TOP-IDP IDRs (10–43,47–66,94–139), no LC | PA (94,104), NPE (111,139); region-level label of 94–139 = PA | MobiDB-lite v3 code run verbatim |
| FS3 | LC beats polar | W10+(SQ)10+W10, LC mask 10–29 | LowComplexity (10,29); without mask Polar (9,30) | `consensus.py` priority G → LC → polar |
| FS4 | Morphology closes short gaps | E12+A3+E12 / E12+A8+E12 | NPE (0,26) / NPE (0,12),(19,31) | `states.py:math_morphology(rmax=5)` |
| FS5 | Min length 10 within IDR | E30, IDR (5,13) / (5,14) | none / NPE (5,14) | `feature_len_thr=10` |
| FS6 | Windowed split | G12+D12 | GlycineRich (0,10), NPE (11,23); region-level NPE | MobiDB-lite v3 code run verbatim |
| FS7 | Short / lowercase | "RK"; 12×e | none; NPE (0,11) | `tokenize` n = L−1 for short input |
| FS8 | Invalid input | "", null, mask length ≠ L, region out of range / End<Start | ArgumentException / ArgumentOutOfRangeException | Input validation |
| R1 | PredictDisorderRegions α-synuclein | P37840 | (10,43,Long IDR),(47,66,Standard IDR),(94,139,Acidic); MeanScore 0.5851639426, 0.5699485887, 0.6185213383; Confidence NaN | Independent Python TOP-IDP recompute (Campen 2008) |
