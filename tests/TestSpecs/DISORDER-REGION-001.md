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
| M4 | MeanScoreIsAverageOfResidueScores | For homo-polymeric disordered region, MeanScore = normalized TOP-IDP value | MeanScore ≈ expected | Campen (2008) |
| M5 | MinLengthFiltering | Disordered region shorter than minLength is excluded | 0 regions | Algorithm definition |
| M6 | TrailingRegionCaptured | Disorder at end of sequence → region includes last residue | Region.End = len - 1 | Algorithm definition |
| M7 | ProlineRichClassification | 30×P → classified as "Proline-rich" | RegionType = "Proline-rich" | van der Lee (2014) |
| M8 | AcidicClassification | 30×E → classified as "Acidic" | RegionType = "Acidic" | van der Lee (2014) |
| M9 | BasicClassification | K/R-rich sequence → classified as "Basic" | RegionType = "Basic" | van der Lee (2014) |
| M10 | SerThrRichClassification | 30×S → classified as "Ser/Thr-rich" | RegionType = "Ser/Thr-rich" | van der Lee (2014) |
| M11 | LongIdrClassification | Long disorder-promoting sequence with no dominant AA → "Long IDR" | RegionType = "Long IDR" | Ward (2004) |
| M12 | StandardIdrClassification | Short disorder-promoting sequence with no dominant AA → "Standard IDR" | RegionType = "Standard IDR" | Fallback |
| M13 | ConfidenceInRange | All region confidences in [0, 1] | Confidence ∈ [0, 1] | Invariant |
| M14 | RegionsNonOverlapping | Multiple regions do not overlap | No start/end intersection | Algorithm definition |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S2 | RegionAtStart | Disorder at very beginning then ordered → region starts at 0 | Region.Start = 0 | Boundary case |
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

- `tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictorTests.cs` — Contains ~12 tests originally kept for DISORDER-REGION-001 scope per DISORDER-PRED-001 consolidation plan.
- `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinderTests.cs` — Contains `PredictDisorderedRegions` tests, but these test `ProteinMotifFinder` not `DisorderPredictor`. Out of scope.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1: All ordered → no regions | ✅ Covered | Exact: Count=0 for 30×W |
| M2: All disordered → one region | ✅ Covered | Exact: Count=1, Start=0, End=29 for 30×P |
| M3: Region boundaries | ✅ Covered | Exact: Count=1, Start=0, End=29 for 30×E |
| M4: MeanScore correctness | ✅ Covered | Exact: MeanScore=1.0 for 30×P |
| M5: MinLength filtering | ✅ Covered | 30×P with minLen=31 → 0 regions |
| M6: Trailing region | ✅ Covered | Exact: Start=11, End=29 for W10+P20 (window boundary 12/21) |
| M7: Proline-rich classification | ✅ Covered | 30×P → "Proline-rich" |
| M8: Acidic classification | ✅ Covered | Exact: Count=1, 30×E → "Acidic" |
| M9: Basic classification | ✅ Covered | Exact: Count=1, 30×K → "Basic" |
| M10: Ser/Thr-rich | ✅ Covered | Exact: Count=1, 30×S → "Ser/Thr-rich" |
| M11: Long IDR | ✅ Covered | Exact: Count=1, Start=0, End=39, (EKQSP)×8 → "Long IDR" |
| M12: Standard IDR | ✅ Covered | Exact: Count=1, Start=0, End=19, (EKQSP)×4 → "Standard IDR" |
| M13: Confidence in range | ✅ Covered | [0,1] range for P, 5×P, E, S |
| M14: Non-overlapping | ✅ Covered | Exact: Count=2, sorted, non-overlapping |
| S2: Region at start | ✅ Covered | P20+W30 → Start=0 |
| S3: Exact min length | ✅ Covered | 5×P with minLen=5 → included |
| S4: Below min length | ✅ Covered | 4×P with minLen=5 → excluded |
| S5: Central region | ✅ Covered | Exact: Start=16, End=33 for W15+P20+W15 |
| C1: Pro > Acidic priority | ✅ Covered | P15+E15 → "Proline-rich" |
| C2: Empty sequence | ✅ Covered | "" → 0 regions |
| C3: Acidic > Basic priority | ✅ Covered | E15+K15 → "Acidic" |
| INV-3: All regions ≥ minLen | ✅ Covered | All regions in multi-region sequence ≥ minLen |
| INV-5: MeanScore exact values | ✅ Covered | Exact: E≈0.866, K≈0.786, S≈0.655 |
| INV-7: Valid labels | ✅ Covered | 6 sequences → all labels in valid set |
| Confidence-High | ✅ Covered | 30×P → Confidence=1.0 |
| Confidence-Lower | ✅ Covered | Exact: S confidence≈0.246; P > S ordering |
| MoRF tests | 🔁 Out of scope | Separate Test Unit (not in DISORDER-REGION-001) |
| LowComplexity tests | 🔁 Out of scope | Separate Test Unit (not in DISORDER-REGION-001) |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs`
  - All DISORDER-REGION-001 tests: region detection, classification, boundaries, edge cases
- **Remove from DisorderPredictorTests.cs:** Region detection tests (#region Region Detection Tests), Region classification tests (#region Region Classification Tests) — these move to canonical file
- **Keep in DisorderPredictorTests.cs:** MoRF tests, LowComplexity tests — these are out of scope for this Test Unit, will be handled in their own Test Units later
- **Remove weak/duplicate tests:** All weak tests replaced by evidence-based versions in canonical file

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `DisorderPredictor_DisorderedRegion_Tests.cs` | DISORDER-REGION-001 canonical | 26 |
| `DisorderPredictorTests.cs` | Future Test Units (MoRF, LowComplexity) | ~5 (kept) |
| `DisorderPredictor_DisorderPrediction_Tests.cs` | DISORDER-PRED-001 canonical | ~50 (unchanged) |

---

## 6. Evidence Traceability

All parameters are traceable to peer-reviewed sources:

| Parameter | Value | Source |
|-----------|-------|--------|
| Classification threshold | 0.25 | ≥5× random single-AA frequency (1/20); for charged pairs matches f+/f− boundary in Das & Pappu (2013) diagram-of-states |
| Classification priority | Pro > Acidic > Basic > S/T > Long > Standard | Most-specific single-AA bias first: P has highest TOP-IDP propensity (Campen 2008) → charge classes (Das & Pappu 2013) → S/T repeats (van der Lee 2014 Table 1) → length (Ward 2004) |
| Long IDR threshold | >30 residues | Ward et al. (2004); van der Lee et al. (2014): "substantial disordered segments of >30 amino acids" |
| Confidence formula | (meanScore − 0.542) / (1.0 − 0.542) | Normalized distance from TOP-IDP decision boundary (Campen 2008) |

---

## 7. Open Questions / Decisions

None. All behavior is traceable to evidence.
