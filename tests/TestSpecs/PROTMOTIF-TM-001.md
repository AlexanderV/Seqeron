# Test Specification: PROTMOTIF-TM-001

**Test Unit ID:** PROTMOTIF-TM-001
**Area:** ProteinMotif
**Algorithm:** Transmembrane Helix Prediction (Kyte-Doolittle hydropathy sliding window)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Kyte & Doolittle (1982), J Mol Biol 157:105-132 | 1 | https://doi.org/10.1016/0022-2836(82)90515-0 | 2026-06-14 |
| 2 | Davidson College — Kyte-Doolittle background | 4 | https://gcat.davidson.edu/DGPB/kd/kyte-doolittle-background.htm | 2026-06-14 |
| 3 | QIAGEN CLC — Hydrophobicity scales | 3 | https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Hydrophobicity_scales.html | 2026-06-14 |
| 4 | Biopython Bio.SeqUtils.ProtParam | 3 | https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/ProtParam.py | 2026-06-14 |

### 1.2 Key Evidence Points

1. TM detection uses a sliding window of **19** residues and a threshold of **1.6** on the window mean — Davidson background page (citing Kyte & Doolittle 1982).
2. The profile point is the **arithmetic mean** of the window's per-residue hydropathy scores, sliding one residue at a time — Davidson page; Biopython `protein_scale` with edge weight 1.0.
3. Kyte-Doolittle scale values (one-letter): I 4.5, V 4.2, L 3.8, F 2.8, C 2.5, M 1.9, A 1.8, G −0.4, T −0.7, S −0.8, W −0.9, Y −1.3, P −1.6, H −3.2, E/Q/D/N −3.5, K −3.9, R −4.5 — QIAGEN and Davidson scale tables (matching exactly).
4. A transmembrane α-helix needs ≈18–21 residues to span the ≈30 Å bilayer — textbook/biology source; justifies the 19-residue window and the minimum-span filter equal to the window width.

### 1.3 Documented Corner Cases

- Window longer than sequence → no profile point computable → no segments (Davidson page).
- Non-standard residues (X, B, Z, *) have no scale value → excluded from the window mean (Biopython scale coverage).

### 1.4 Known Failure Modes / Pitfalls

1. Mapping the threshold-crossing window **run** back to residue boundaries is not prescribed by the source; the reported `End` uses a stated implementation convention (`lastProfileIndex + windowSize`) — Evidence §Assumptions.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `PredictTransmembraneHelices(string, int, double)` | ProteinMotifFinder | Canonical | Kyte-Doolittle window/threshold scan |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every reported segment `Score` (peak window mean) ≥ `threshold`. | Yes | Detection rule (Davidson page) |
| INV-2 | Every reported segment has `0 ≤ Start ≤ End ≤ length−1`. | Yes | Output-coordinate contract |
| INV-3 | A uniform sequence of one residue (length ≥ window) yields a profile equal to that residue's scale value at every point. | Yes | Window mean = mean of identical values (QIAGEN scale) |
| INV-4 | Null / empty / shorter-than-window input yields no segments. | Yes | Window undefined below `windowSize` |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Single hydrophobic stretch | `D`×10 + `L`×20 + `D`×10, window 19, threshold 1.6 | Exactly one segment (Start=5, End=34), Score=3.8 | Davidson rule + KD scale (QIAGEN); End = last residue covered by a passing window |
| M2 | All-hydrophilic sequence | `D`×40 (D=−3.5 < 1.6) | No segments | KD scale value D=−3.5 |
| M3 | Exactly one window of poly-Leu | `L`×19 | One segment (Start=0, End=18, Score=3.8) | Window mean = 3.8 ≥ 1.6 |
| M4 | Scale value reproduction | 19-residue uniform window for each of I, V, R reproduces its KD value as the peak score | Score = 4.5 (I), 4.2 (V); R → no segment (−4.5 < 1.6) | QIAGEN / Davidson scale tables |
| M5 | Null input | `null` | Empty result | Window undefined |
| M6 | Empty input | `""` | Empty result | Window undefined |
| M7 | Shorter than window | length 18 < 19 | Empty result | Window undefined |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Non-standard residue excluded | window with one `X` among 18 `L` (window 19) → mean over the 18 L = 3.8 | Segment with Score=3.8 | Biopython scale coverage |
| S2 | Custom higher threshold suppresses weak segment | `D`×10 + `A`×20 + `D`×10 (A=1.8), threshold 2.0 | No segment (peak mean < 2.0) | Threshold is a stated parameter |
| S3 | Lowercase input accepted | lowercase of M1 sequence | Same segment as M1 | Case-insensitive contract |
| S4 | Non-positive window | windowSize = 0 | Empty result | Guarded invalid input |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Invariant (property) | For a multi-residue sequence with segments, every Score ≥ threshold and Start ≤ End | All segments satisfy INV-1/INV-2 | Property-based |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Searched `tests/Seqeron/Seqeron.Genomics.Tests/` for `PredictTransmembraneHelices`. No dedicated test file exists. `ProteinMotifFinderTests.cs` and `ProteinMotifFinder_MotifSearch_Tests.cs` do not exercise this method. This is a new unit; all planned cases start ❌ Missing.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new unit |
| M2 | ❌ Missing | new unit |
| M3 | ❌ Missing | new unit |
| M4 | ❌ Missing | new unit |
| M5 | ❌ Missing | new unit |
| M6 | ❌ Missing | new unit |
| M7 | ❌ Missing | new unit |
| S1 | ❌ Missing | new unit |
| S2 | ❌ Missing | new unit |
| S3 | ❌ Missing | new unit |
| S4 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_PredictTransmembraneHelices_Tests.cs` — all cases above.
- **Remove:** none (no prior tests for this method).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `ProteinMotifFinder_PredictTransmembraneHelices_Tests.cs` | Canonical | 12 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | S1 | ❌ Missing | Implemented | ✅ Done |
| 9 | S2 | ❌ Missing | Implemented | ✅ Done |
| 10 | S3 | ❌ Missing | Implemented | ✅ Done |
| 11 | S4 | ❌ Missing | Implemented | ✅ Done |
| 12 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 12
**✅ Done:** 12 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | exact (Start,End,Score) asserted |
| M2 | ✅ Covered | empty asserted |
| M3 | ✅ Covered | exact (0,18,3.8) asserted |
| M4 | ✅ Covered | I/V scores + R no-segment asserted |
| M5 | ✅ Covered | null → empty |
| M6 | ✅ Covered | empty → empty |
| M7 | ✅ Covered | len 18 → empty |
| S1 | ✅ Covered | X excluded, score 3.8 |
| S2 | ✅ Covered | threshold 2.0 → no segment |
| S3 | ✅ Covered | lowercase == uppercase |
| S4 | ✅ Covered | window 0 → empty |
| C1 | ✅ Covered | property over a real sequence |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Segment `End` reported as `lastProfileIndex + windowSize − 1` (clamped) = last residue covered by an above-threshold window — output coordinate convention, not a detection-rule change | M1, M3 expected `End` values |

---

## 7. Open Questions / Decisions

1. None. The detection rule (window 19, threshold 1.6, arithmetic-mean profile) and scale values are fully source-backed; the only non-source-prescribed item (segment-end coordinate mapping) is documented as an output convention and does not affect which residues qualify as transmembrane.
