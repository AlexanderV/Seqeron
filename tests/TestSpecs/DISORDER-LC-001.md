# Test Specification: DISORDER-LC-001

**Test Unit ID:** DISORDER-LC-001
**Area:** ProteinPred
**Algorithm:** Low Complexity Region Detection (SEG)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-02-12

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Wootton & Federhen (1993) Non-globular domains in protein sequences | 1 | https://doi.org/10.1016/0097-8485(93)85006-X | 2026-02-12 |
| 2 | Wootton & Federhen (1996) Analysis of compositionally biased regions | 1 | https://doi.org/10.1016/s0076-6879(96)66035-2 | 2026-02-12 |
| 3 | Shannon (1948) A Mathematical Theory of Communication | 1 | Bell System Technical Journal 27:379-423 | 2026-02-12 |

### 1.2 Key Evidence Points

1. SEG identifies low complexity regions via Shannon entropy — Wootton & Federhen (1993, 1996)
2. Two-pass approach: K1 trigger scan (window=12, threshold=2.2 bits), K2 extension (threshold=2.5 bits) — Wootton & Federhen (1993, 1996)
3. Shannon entropy: H = −Σ p_i · log₂(p_i) in bits, not normalized — Shannon (1948)
4. Overlapping triggered regions are merged — Wootton & Federhen (1993)
5. Default parameters: triggerWindow=12, K1=2.2, K2=2.5 — standard SEG defaults

### 1.3 Documented Corner Cases

1. Sequence shorter than triggerWindow → no regions
2. Homopolymer (entropy = 0) → entire sequence is one region
3. Maximum entropy (all 20 AAs in window) → no trigger
4. Region at end of sequence must be captured
5. K2 extension can pull additional residues beyond triggered span

### 1.4 Known Failure Modes / Pitfalls

1. Floating-point comparison at exact threshold boundaries
2. K2 extension + merge can produce regions whose total entropy slightly exceeds K2
3. Short sequences near triggerWindow length may have edge effects

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `PredictLowComplexityRegions(sequence, triggerWindow, triggerThreshold, extensionThreshold, minLength)` | DisorderPredictor | Public API | SEG algorithm entry point |
| `CalculateShannonEntropy(sequence, start, length)` | DisorderPredictor | Private helper | Raw Shannon entropy in bits |
| `MergeSegments(segments)` | DisorderPredictor | Private helper | Overlapping/adjacent merge |
| `ClassifyLowComplexityType(sequence, start, length)` | DisorderPredictor | Private helper | Dominant AA classification |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | All regions have Start ≥ 0, End < sequence.Length, Start ≤ End | Yes | Algorithm definition |
| INV-2 | Regions are non-overlapping and sorted by Start | Yes | Merge step guarantees |
| INV-3 | Region Type is non-null and non-empty | Yes | ClassifyLowComplexityType always returns a string |

---

## 4. Test Cases

### 4.1 MUST Tests (Required)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| S1 | Homopolymer_DetectedAsLC | 26×Q → entropy = 0 bits in all windows | 1 region (0, 25), "Q-rich" | Wootton (1993): H=0 < K1=2.2 |
| S2 | ComplexSequence_NoLC | All 20 AAs × 2 = 40 AA → H ≥ 3.58 bits | 0 regions | Wootton (1993): H >> K1 |
| S3 | DipeptideBlock_DetectedAsLC | 12A + 12L → max H = log₂(2) = 1.0 | 1 region (0, 23) | Wootton (1993): H < K1 |
| M1 | FourTypesTriggersK1 | AAABBBCCCDDD → H = log₂(4) = 2.0 ≤ 2.2 | 1 region | Wootton (1993) |
| M2 | TwelveDistinctNoTrigger | 12 distinct AAs → H = log₂(12) = 3.58 | 0 regions | Wootton (1993) |
| M3 | TwoSeparatedRegions | Poly-Q + long separator + Poly-A | Q-rich and A-rich both found | Wootton (1993): merge step |
| M6 | TypeClassification | 20×A → "A-rich" for single dominant AA | Type = "A-rich" | Internal classification |

### 4.2 SHOULD Tests (Edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | ShorterThanWindow | 6×Q with triggerWindow=12 | 0 regions | Minimum length guard |
| C2 | ExactWindowLength | 12×Q → exactly one window | 1 region (0, 11) | Boundary |
| C3 | EmptySequence | "" | 0 regions | Trivial |
| C4 | MinLengthFilters | 12×Q with minLength=15 | 0 regions | Post-filter |
| M4 | CustomThreshold | K1=0.5 on AAABBBCCCDDD (H=2.0 > 0.5) | 0 regions | Parameter override |
| M5 | CaseInsensitive | lower vs upper → same results | Identical regions | API contract |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_LowComplexity_Tests.cs` — 15 tests for DISORDER-LC-001.
- Previous tests in `DisorderPredictorTests.cs` removed (used old API with wrong parameters).

### 5.2 Coverage Classification

| Test Case ID | Status | Notes |
|--------------|--------|-------|
| S1: Homopolymer | ✅ Covered | Exact: Count=1, (0,25), "Q-rich" |
| S2: Complex | ✅ Covered | Exact: Is.Empty |
| S3: Dipeptide | ✅ Covered | Exact: Count=1, (0,23) |
| C1: Short | ✅ Covered | Is.Empty |
| C2: Exact window | ✅ Covered | Exact: Count=1, (0,11), "Q-rich" |
| C3: Empty | ✅ Covered | Is.Empty |
| C4: MinLength | ✅ Covered | Is.Empty with minLength=15 |
| M1: Four types | ✅ Covered | Count=1 |
| M2: Twelve distinct | ✅ Covered | Is.Empty |
| M3: Two separated | ✅ Covered | Q-rich + A-rich both found |
| M4: Custom threshold | ✅ Covered | Is.Empty with K1=0.5 |
| M5: Case insensitive | ✅ Covered | Same count/positions |
| M6: Type classification | ✅ Covered | "A-rich" |
| INV-1: Valid boundaries | ✅ Covered | 3 diverse inputs |
| INV-2: No overlapping | ✅ Covered | 2 inputs |

### 5.3 Final State

| File | Role | Test Count |
|------|------|------------|
| `DisorderPredictor_LowComplexity_Tests.cs` | DISORDER-LC-001 canonical | 15 |

---

## 6. Evidence Traceability

| Parameter | Value | Source | Status |
|-----------|-------|--------|--------|
| SEG trigger window | 12 residues | Wootton & Federhen (1993, 1996) | ✅ Standard default |
| K1 trigger threshold | 2.2 bits | Wootton & Federhen (1993, 1996) | ✅ Standard default |
| K2 extension threshold | 2.5 bits | Wootton & Federhen (1993, 1996) | ✅ Standard default |
| Entropy formula | H = −Σ p_i · log₂(p_i) | Shannon (1948) | ✅ Standard |
| Merge step | overlapping/adjacent segments unified | Wootton & Federhen (1993) | ✅ Standard |
| Type classification | dominant AA fraction > 0.5 → "X-rich" | **Internal heuristic** | ⚠️ Design decision |

---

## 7. Deviations and Assumptions

| ID | Item | Status | Detail |
|----|------|--------|--------|
| D1 | Type classification heuristic | ⚠️ Internal | Dominant AA > 50% = "X-rich", otherwise top two AAs. SEG itself does not classify regions — it only identifies them. |
| D2 | minLength default = 1 | ⚠️ Internal | Original SEG has no minimum length filter. Parameter added for practical use. |
