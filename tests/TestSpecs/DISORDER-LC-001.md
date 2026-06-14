# Test Specification: DISORDER-LC-001

**Test Unit ID:** DISORDER-LC-001
**Area:** ProteinPred
**Algorithm:** Low-Complexity Region Detection in Protein Sequences (SEG; Wootton & Federhen)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Wootton & Federhen (1993), Comput. Chem. 17(2):149–163 | 1 | https://doi.org/10.1016/0097-8485(93)85006-X | 2026-06-14 |
| 2 | Wootton & Federhen (1996), Methods Enzymol 266:554–571 | 1 | https://doi.org/10.1016/S0076-6879(96)66035-2 | 2026-06-14 |
| 3 | NCBI BLAST `blast_seg.c` (reference impl) | 3 | https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/blast__seg_8c.html | 2026-06-14 |
| 4 | SEG program help (GCG mirror) + `ncbi-seg` manpage | 3 | https://bip.weizmann.ac.il/education/materials/gcg/seg.html | 2026-06-14 |

### 1.2 Key Evidence Points

1. Default parameters: trigger window W = 12, trigger complexity K1 = 2.2 bits, extension complexity K2 = 2.5 bits — Sources 3 (`kSegWindow=12`, `kSegLocut=2.2`, `kSegHicut=2.5`), 4 (`-WINdow=12`, `-LOWcut=2.2`, `-HIGhcut=2.5`).
2. Complexity = Shannon entropy of the window residue composition, in bits/residue: H = −Σ pᵢ·log₂(pᵢ) — Source 3 (`s_Entropy`, normalized by window length, base-2 via `NCBIMATH_LN2`); Source 1/2 ("Shannon's Entropy").
3. Maximum complexity for a 20-letter amino-acid alphabet = log₂(20) ≈ 4.322 bits/residue — Source 4.
4. Two-stage scan: stage 1 marks windows with complexity ≤ K1 (trigger); stage 2 extends triggered segments while complexity ≤ K2 — Source 4.

### 1.3 Documented Corner Cases

- Sequence shorter than W: no full trigger window → no segments (Source 4, window = minimum first-stage segment size).
- W-distinct window: H = log₂(W) ≈ 3.585 for W=12 > K2 → never low-complexity (Source 4 max-complexity statement).

### 1.4 Known Failure Modes / Pitfalls

1. The entropy must normalize by the window length; only standard residues contribute to the composition (Source 3, `s_Entropy` over the composition state vector).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `PredictLowComplexityRegions(string, int, double, double, int)` | DisorderPredictor | **Canonical** | SEG trigger+extension; deep evidence-based testing |
| `CalculateShannonEntropy` (private) | DisorderPredictor | **Internal** | Tested indirectly via segment boundaries |
| `ClassifyLowComplexityType` (private) | DisorderPredictor | **Internal** | Convenience label; tested via region `Type` |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every returned segment has 0 ≤ Start ≤ End < sequence length | Yes | Coordinate contract (0-based inclusive) |
| INV-2 | Returned segments are non-overlapping and non-adjacent (merged) | Yes | Source 4 (overlapping segments merged into contigs) |
| INV-3 | A window of L identical residues has complexity 0 ≤ K1 (always triggers) | Yes | Source 3 (`s_Entropy`); H=0 for single residue |
| INV-4 | A window of W distinct residues has complexity log₂(W) > K2 (never low-complexity) | Yes | Source 4 (max complexity log₂20) |
| INV-5 | Every segment length ≥ `minLength` | Yes | minimum-length post-filter |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Homopolymer ≥ W | 26×Q; every window H=0 ≤ K1 | exactly 1 segment (0, 25), Type "Q-rich" | Source 3 (`s_Entropy`=0), Source 4 |
| M2 | Max-complexity sequence | all-20 AA ×2 (40 AA); every 12-window has ≥12 distinct, H ≥ log₂12 ≈ 3.585 > K2 | empty | Source 4 max-complexity; H=3.585 |
| M3 | Four-types window triggers at default | `AAABBBCCCDDD`, H=2.0 ≤ K1=2.2 | exactly 1 segment (0,11) | hand-derived H=2.0 vs K1 |
| M4 | Four-types window does NOT trigger at strict K1 | same input, K1=0.5; H=2.0 > 0.5 | empty | hand-derived H=2.0 vs strict cutoff |
| M5 | Two-residue block | 12×A + 12×L; each window H ≤ 1.0 ≤ K1 | 1 merged segment (0, 23) | hand-derived H=1.0 |
| M6 | Two separated runs | 20×Q + 60-AA high-complexity spacer + 20×A | exactly 2 segments (0,34) Q-rich and (67,99) A-rich, separated by a high-complexity gap | stage1+2 semantics; spacer H > K2; boundaries hand-derived from trigger spans (0,24)/(75,99) extended while seg H ≤ K2 |
| M7 | Dominant-residue label > 50 % | 20×A | Type "A-rich" | dominant-fraction label rule |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | minLength filter | 12×Q with minLength=15 | empty (length 12 < 15) | post-filter |
| S2 | Sequence shorter than W | 6×Q | empty | no full trigger window |
| S3 | Exactly window length | 12×Q | 1 segment (0, 11), "Q-rich" | single window boundary |
| S4 | Case-insensitivity | lowercase vs uppercase poly-Q | identical segments | input upper-cased |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Empty string | "" | empty | length 0 < W |
| C2 | Invariant: boundaries | property over several inputs | all segments satisfy INV-1 | property-based |
| C3 | Invariant: non-overlap | property over several inputs | segments satisfy INV-2 | property-based |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_LowComplexity_Tests.cs` — pre-existing fixture for this method (S1–S3, C1–C4, M1–M6, INV1–INV2). Not previously backed by a template Evidence/TestSpec pair; several assertions are permissive: `M3_TwoSeparatedLCRegions` used `Contains` substring checks rather than exact boundaries/count; `M1`/several others used `.Count.EqualTo(1)` without exact-entropy evidence notes; `M5` (case) compared counts only.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 Homopolymer (existing S1) | ⚠ Weak | Right values but no exact-entropy evidence note; rewrite under canonical IDs |
| M2 Max-complexity (existing S2/M2) | ⚠ Weak | Rewrite with H=3.585 message and a 40-AA worked input |
| M3 Four-types triggers (existing M1) | ⚠ Weak | `Count.EqualTo(1)` only; add exact H=2.0 message + boundary (0,11) |
| M4 Strict-K1 no trigger (existing M4) | ⚠ Weak | OK logic; rewrite with H=2.0 vs 0.5 message |
| M5 Two-residue block (existing S3) | ⚠ Weak | Lacked exact H message |
| M6 Two separated runs (existing M3) | ⚠ Weak | `Contains` substring; rewrite with exact count=2 + boundary checks |
| M7 Dominant label (existing M6) | ⚠ Weak | OK; rewrite under canonical ID |
| S1 minLength (existing C4) | ✅ Covered | Keep semantics, restate |
| S2 short<W (existing C1) | ✅ Covered | Keep |
| S3 exactly W (existing C2) | ✅ Covered | Keep |
| S4 case-insensitive (existing M5) | ⚠ Weak | Compared counts only; tighten to identical boundaries |
| C1 empty (existing C3) | ✅ Covered | Keep |
| C2 INV boundaries (existing INV1) | ⚠ Weak | Permissive; keep as property but add message detail |
| C3 INV non-overlap (existing INV2) | ✅ Covered | Keep |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_LowComplexity_Tests.cs` — rewrite as the single canonical fixture for DISORDER-LC-001 with evidence-grounded IDs (M/S/C/INV), exact entropy values, and assertion messages.
- **Remove:** none (the broad `DisorderPredictorTests.cs` does not test this method; no duplicate file).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `DisorderPredictor_LowComplexity_Tests.cs` | Canonical fixture for DISORDER-LC-001 | 14 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ⚠ Weak | Rewrote with H=0 evidence + exact boundary | ✅ Done |
| 2 | M2 | ⚠ Weak | Rewrote with H=3.585 message, 40-AA input | ✅ Done |
| 3 | M3 | ⚠ Weak | Rewrote with H=2.0 + exact boundary (0,11) | ✅ Done |
| 4 | M4 | ⚠ Weak | Rewrote with H=2.0 vs 0.5 message | ✅ Done |
| 5 | M5 | ⚠ Weak | Rewrote with H=1.0 + boundary (0,23) | ✅ Done |
| 6 | M6 | ⚠ Weak | Rewrote with exact count=2 + per-region boundary | ✅ Done |
| 7 | M7 | ⚠ Weak | Rewrote dominant label | ✅ Done |
| 8 | S4 | ⚠ Weak | Tightened to identical boundaries | ✅ Done |
| 9 | C2 | ⚠ Weak | Kept property test with detailed messages | ✅ Done |

**Total items:** 9
**✅ Done:** 9 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | Exact boundary + Type, evidence H=0 |
| M2 | ✅ | Empty result, H=3.585 > K2 |
| M3 | ✅ | 1 segment (0,11), H=2.0 ≤ K1 |
| M4 | ✅ | Empty at strict K1=0.5 |
| M5 | ✅ | Merged (0,23), H=1.0 |
| M6 | ✅ | Two segments, exact boundaries |
| M7 | ✅ | Type "A-rich" |
| S1 | ✅ | minLength filter |
| S2 | ✅ | short<W empty |
| S3 | ✅ | exactly W → (0,11) |
| S4 | ✅ | identical boundaries lower/upper |
| C1 | ✅ | empty string |
| C2 | ✅ | INV-1 boundaries property |
| C3 | ✅ | INV-2 non-overlap property |

**In-scope cases:** 14 | **✅:** 14

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Region-type label (`"X-rich"`/`"X/Y-rich"`) is a presentation extension (> 50 % dominant-residue rule) | M1, M7 (Type assertions) |
| 2 | Greedy single-residue extension equals merged-window contig growth for the test inputs | M5, M6 boundaries |

---

## 7. Open Questions / Decisions

1. None — the core SEG parameters and complexity measure are fully source-backed; the only non-source items (textual label, extension granularity) are documented assumptions and do not affect segment boundaries for the tested inputs.
