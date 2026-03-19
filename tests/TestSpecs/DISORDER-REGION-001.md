# Test Specification: DISORDER-REGION-001

**Test Unit ID:** DISORDER-REGION-001
**Area:** ProteinPred
**Algorithm:** Disordered Region Detection
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-02-12

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

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S2 | RegionAtStart | P(20)+W(30): region at start with exact End from window boundary | Count=1, Start=0, End=18 | Window boundary: 12P/21 at pos 18 |
| S3 | ExactMinLength | Region of exactly minLength residues → included | 1 region | Boundary case |
| S4 | JustBelowMinLength | Region of (minLength - 1) residues → excluded | 0 regions | Boundary case |
| S5 | MixedSequenceRegionDetection | Ordered-disordered-ordered → identifies central region | Start=16, End=33 | Exact window boundaries |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | ClassificationPriority | Sequence with both P>0.25 and E/D>0.25 → "Proline-rich" wins | RegionType = "Proline-rich" | Implementation priority |
| C2 | EmptySequence | PredictDisorder("") → no regions | 0 regions | Trivial |
| C3 | AcidicOverBasicPriority | Sequence with both E/D>0.25 and K/R>0.25 → "Acidic" wins | RegionType = "Acidic" | Priority chain verification |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs` — Canonical file: 24 tests for DISORDER-REGION-001.
- `tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictorTests.cs` — MoRF + LowComplexity tests only. Out of scope.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1: All ordered → no regions | ✅ Covered | Exact: Count=0 for 30×W |
| M2: All disordered → one region | ✅ Covered | Exact: Count=1, Start=0, End=29 for 30×P |
| M3: Region boundaries | ✅ Covered | Exact: Count=1, Start=0, End=29 for 30×E |
| M5: MinLength filtering | ✅ Covered | 30×P with minLen=31 → 0 regions |
| M6: Trailing region | ✅ Covered | Exact: Start=11, End=29 for W10+P20 (12P/21 → 0.571 > 0.542) |
| M7: Proline-rich classification | ✅ Covered | 30×P → "Proline-rich" |
| M8: Acidic classification | ✅ Covered | Exact: Count=1, 30×E → "Acidic" |
| M9: Basic classification | ✅ Covered | Exact: Count=1, 30×K → "Basic" |
| M10: Ser/Thr-rich | ✅ Covered | Exact: Count=1, 30×S → "Ser/Thr-rich" |
| M11: Long IDR | ✅ Covered | Exact: Count=1, Start=0, End=39, (EKQSP)×8 → "Long IDR" |
| M12: Standard IDR | ✅ Covered | Exact: Count=1, Start=0, End=19, (EKQSP)×4 → "Standard IDR" |
| M13: Confidence in range | ✅ Covered | [0,1] range for P(30), P(5), E(30), S(30) |
| M14: Non-overlapping | ✅ Covered | Exact: Count=2, sorted, non-overlapping |
| S2: Region at start | ✅ Covered | Exact: Count=1, Start=0, End=18 for P20+W30 (12P/21 at pos 18) |
| S3: Exact min length | ✅ Covered | 5×P with minLen=5 → included |
| S4: Below min length | ✅ Covered | 4×P with minLen=5 → excluded |
| S5: Central region | ✅ Covered | Exact: Start=16, End=33 for W15+P20+W15 |
| C1: Pro > Acidic priority | ✅ Covered | P15+E15 → "Proline-rich" |
| C2: Empty sequence | ✅ Covered | "" → 0 regions |
| C3: Acidic > Basic priority | ✅ Covered | E15+K15 → "Acidic" |
| INV-1/2: Bounds invariant | ✅ Covered | Start≥0, End<Length, End≥Start across 6 diverse inputs |
| INV-3: All regions ≥ minLen | ✅ Covered | All regions in multi-region sequence ≥ minLen |
| INV-5: MeanScore exact values | ✅ Covered | Exact: P=1.0, E≈0.866, K≈0.786, S≈0.655 |
| Confidence exact values | ✅ Covered | Exact: P=1.0, S≈0.246; P > S ordering |

### 5.3 Removed Tests (Duplicates)

| Test | Reason | Subsumed by |
|------|--------|-------------|
| M4: MeanScoreIsAverage (P→1.0) | Same assertion as INV-5 first line; Count=1 already in M2 | INV-5 |
| INV-7: ValidLabels | All 6 labels explicitly tested by M7–M12 | M7–M12 |
| Confidence-High (P→1.0) | Same P=1.0 assertion already in Confidence-Lower | Confidence-Lower |

### 5.4 Final State

| File | Role | Test Count |
|------|------|------------|
| `DisorderPredictor_DisorderedRegion_Tests.cs` | DISORDER-REGION-001 canonical | 24 |
| `DisorderPredictorTests.cs` | Future Test Units (MoRF, LowComplexity) | 5 |
| `DisorderPredictor_DisorderPrediction_Tests.cs` | DISORDER-PRED-001 canonical | ~50 (unchanged) |

---

## 6. Evidence Traceability

| Parameter | Value | Source | Status |
|-----------|-------|--------|--------|
| TOP-IDP scale values | 20 AA propensities | Campen et al. (2008) Table 2 | ✅ Verified |
| TOP-IDP cutoff | 0.542 | Campen et al. (2008) maximum-likelihood | ✅ Verified |
| Window size | 21 residues | Campen et al. (2008) web server | ✅ Verified |
| Disorder/Order AA sets | {A,R,G,Q,S,P,E,K} / {W,C,F,I,Y,V,L,N} | Dunker et al. (2001) | ✅ Verified |
| Long IDR threshold | >30 residues | Ward et al. (2004); van der Lee et al. (2014) | ✅ Verified |
| IDR subtype names | Proline-rich, Acidic, Basic, Ser/Thr-rich | van der Lee et al. (2014) — recognized subtypes | ✅ Verified |
| Classification threshold | 0.25 (= 5× random 1/20) | **Internal heuristic** — no published source | ⚠️ Design decision |
| Classification priority | Pro > Acidic > Basic > S/T > Long > Standard | **Internal heuristic** — no published source | ⚠️ Design decision |
| Confidence formula | (meanScore − 0.542) / (1.0 − 0.542) | **Internal heuristic** — not from Campen (2008) | ⚠️ Design decision |
| AA classification groups | E+D, K+R, S+T | Standard biochemistry — no IDR-specific source | ⚠️ Design decision |

---

## 7. Deviations and Assumptions

| ID | Item | Status | Detail |
|----|------|--------|--------|
| D1 | Classification enrichment threshold 0.25 | ⚠️ Internal | Defined as 5× random single-AA frequency (1/20). No published source. Previously falsely attributed to Das & Pappu (2013) f+/f− boundary; that paper's 0.25 is NCPR (net charge per residue) for globule/coil conformational state, an unrelated concept. |
| D2 | Classification priority order | ⚠️ Internal | Pro > Acidic > Basic > S/T > Long > Standard. No paper defines this ordering. Van der Lee (2014) lists IDR subtypes but provides no algorithmic classification scheme with priority. |
| D3 | Confidence formula | ⚠️ Internal | (meanScore − 0.542) / (1.0 − 0.542), clamped [0,1]. Campen (2008) defines prediction equation I = −(⟨Top-IDP⟩ − 0.542) but no confidence metric. This formula is a linear rescaling above the cutoff. |
| D4 | AA classification groups | ⚠️ Conventional | {E,D}=Acidic, {K,R}=Basic, {S,T}=Ser/Thr-rich. Standard biochemical side-chain property groupings. No IDR-specific paper mandates these exact group compositions for classification. |
