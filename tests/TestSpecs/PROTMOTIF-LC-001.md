# Test Specification: PROTMOTIF-LC-001

**Test Unit ID:** PROTMOTIF-LC-001
**Area:** ProteinMotif
**Algorithm:** Low-Complexity Region Detection (SEG, Wootton & Federhen 1993)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Wootton & Federhen (1993), Computers & Chemistry 17(2):149–163 | 1 | https://doi.org/10.1016/0097-8485(93)85006-X | 2026-06-14 |
| 2 | NCBI `ncbi-seg` man page | 2 | https://manpages.ubuntu.com/manpages/focal/man1/ncbi-seg.1.html | 2026-06-14 |
| 3 | NCBI C++ Toolkit `blast_seg.c` (s_Entropy, kSegWindow/Locut/Hicut) | 3 | https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/blast__seg_8c.html | 2026-06-14 |
| 4 | SeqComplex `SeqComplex.pm` (subs `ce`, `cwf`, `log_k`) | 3 | https://github.com/caballero/SeqComplex | 2026-06-14 |
| 5 | Mier et al., Bioinformatics 22(24):2980 (Shannon entropy −Σ pᵢ log pᵢ) | 1 | https://academic.oup.com/bioinformatics/article/22/24/2980/208627 | 2026-06-14 |
| 6 | Pei & Grishin, Bioinformatics 21(2):160 (SEG two-pass + defaults) | 1 | https://academic.oup.com/bioinformatics/article/21/2/160/187330 | 2026-06-14 |

### 1.2 Key Evidence Points

1. SEG measures local complexity over a sliding window; complexity is in **bits per residue**, range 0 to log₂(20) = 4.322 for amino acids — Source 2.
2. Default parameters: **W = 12, K1 (trigger) = 2.2 bits, K2 (extension) = 2.5 bits** — Sources 2, 3, 6.
3. The operational complexity is **Shannon entropy K = −Σ pᵢ·log₂ pᵢ** over residue counts in the window (NCBI `s_Entropy`; SeqComplex `ce`: `ce -= r * log2(r)`) — Sources 3, 4, 5.
4. SEG is a **two-pass algorithm**: trigger windows with complexity ≤ K1 mark raw low-complexity segments, extended while complexity ≤ K2 — Sources 3, 6.
5. SEG "reflects residue composition … with no regard of patterns or periodicity" — Source 6.

### 1.3 Documented Corner Cases

- Sequence shorter than window W: no complete trigger window exists → no segment can be triggered (Source 2).
- Homopolymer window: single symbol → entropy 0 (lowest possible complexity) (Sources 4, 5).
- K2 must exceed K1 to be effective in extension (Source 2).

### 1.4 Known Failure Modes / Pitfalls

1. Confusing the bits/residue entropy form with the normalized [0,1] WF variant (universalmotif) — units differ; the official spec uses bits/residue — Sources 2, 5.
2. Using an ad-hoc "dominant single-residue frequency ≥ threshold" rule instead of an entropy/complexity measure — not traceable to Wootton & Federhen; corrected in this unit.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindLowComplexityRegions(string, int windowSize, double triggerComplexity, double extensionComplexity)` | ProteinMotifFinder | **Canonical** | SEG bits/residue complexity, two-pass trigger/extension |
| `CalculateSegComplexity(ReadOnlySpan<char>)` | ProteinMotifFinder | **Internal** | per-window Shannon entropy in bits/residue; tested indirectly via region complexity |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Per-window complexity K ∈ [0, log₂(20)] for amino-acid input | Yes | Source 2 (range 0–4.322) |
| INV-2 | A homopolymer window has complexity exactly 0 | Yes | Sources 4, 5 (single symbol entropy = 0) |
| INV-3 | Reported region complexity ≤ extensionComplexity (K2); the triggering minimum ≤ triggerComplexity (K1) | Yes | Sources 2, 3 (two-pass cutoffs) |
| INV-4 | Region boundaries are 0-based inclusive and lie within the sequence | Yes | windowed scan over sequence |
| INV-5 | Result is deterministic and order-independent (left-to-right by Start) | Yes | deterministic scan |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Homopolymer complexity | `CalculateSegComplexity` on 12×'A' | 0.0 bits | INV-2; Sources 4,5 |
| M2 | Biased window complexity | 11×'A'+1×'B', L=12 | 0.413817 bits (±1e-6) | K=−Σpᵢlog₂pᵢ; Source 5 |
| M3 | Equal-split window | 6×'A'+6×'B', L=12 | 1.0 bits | formula; Source 5 |
| M4 | Asymmetric window | 10×'A'+2×'B', L=12 | 0.650022 bits | formula; Source 5 |
| M5 | Max-diversity window | 12 distinct residues | log₂(12)=3.584963 bits | formula; Source 5 |
| M6 | Poly-Q tract detected | diverse flank + 20×'Q' + diverse flank | exactly one region spanning the Q tract; complexity ≤ K2 | Sources 2,3,6 |
| M7 | Diverse protein, no region | 20 distinct residues only | empty result (all windows > K2) | M5 + K2=2.5; Sources 2,6 |
| M8 | Default parameters | call with defaults | W=12, K1=2.2, K2=2.5 govern detection (poly-A 12-mer in 12-mer-flanked seq detected) | Sources 2,3,6 |
| M9 | Two separated tracts | poly-G tract + diverse spacer + poly-S tract | two distinct regions, correct dominant composition order | Source 6 (independent segments) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Sequence shorter than window | length < W | empty result | no complete window (Source 2) |
| S2 | Region within bounds | poly-tract region Start≥0, End<len, Start≤End | bounds valid | INV-4 |
| S3 | Lowercase input | mixed-case poly tract | detected (case-insensitive) | normalization to upper |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Complexity bound | any region complexity in [0, log₂20] | within bounds | INV-1 |
| C2 | Determinism | two runs equal | identical regions | INV-5 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinderTests.cs` (lines 103–136) — three tests for the *old* `FindLowComplexityRegions` using the invented "dominant-AA frequency ≥ threshold" rule and `DominantAa`/`Frequency` fields. These encode a non-conforming behavior and must be removed (superseded by the SEG implementation and the new canonical file).
- No `ProteinMotifFinder_FindLowComplexityRegions_Tests.cs` existed prior to this unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| Old `FindLowComplexityRegions_PolyAlanine_Finds` | 🔁 Duplicate / non-conforming | tests removed invented `Frequency` field; superseded by M1/M6 |
| Old `FindLowComplexityRegions_Diverse_ReturnsEmpty` | ⚠ Weak | existence-only, invented rule; superseded by M7 |
| Old `FindLowComplexityRegions_MultipleRegions_FindsAll` | ⚠ Weak | invented `DominantAa`; superseded by M9 |
| M1–M9, S1–S3, C1–C2 | ❌ Missing | new canonical file |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_FindLowComplexityRegions_Tests.cs` — all SEG evidence-based tests.
- **Remove:** the `Low Complexity Tests` region in `ProteinMotifFinderTests.cs` (non-conforming, refers to removed `DominantAa`/`Frequency`).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `ProteinMotifFinder_FindLowComplexityRegions_Tests.cs` | canonical | 14 |
| `ProteinMotifFinderTests.cs` (Low Complexity region) | removed | 0 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented homopolymer complexity test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented 11A/1B complexity test | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented 6A/6B complexity test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented 10A/2B complexity test | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented 12-distinct complexity test | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented poly-Q detection test | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented diverse no-region test | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented default-parameter test | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented two-tract test | ✅ Done |
| 10 | S1 | ❌ Missing | Implemented short-sequence test | ✅ Done |
| 11 | S2 | ❌ Missing | Implemented bounds test | ✅ Done |
| 12 | S3 | ❌ Missing | Implemented lowercase test | ✅ Done |
| 13 | C1 | ❌ Missing | Implemented complexity-bound test | ✅ Done |
| 14 | C2 | ❌ Missing | Implemented determinism test | ✅ Done |
| 15 | Old tests | 🔁/⚠ | Removed Low Complexity region from ProteinMotifFinderTests.cs | ✅ Done |

**Total items:** 15
**✅ Done:** 15 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | exact 0.0 |
| M2 | ✅ Covered | 0.413817 ±1e-6 |
| M3 | ✅ Covered | 1.0 |
| M4 | ✅ Covered | 0.650022 ±1e-6 |
| M5 | ✅ Covered | 3.584963 ±1e-6 |
| M6 | ✅ Covered | single poly-Q region |
| M7 | ✅ Covered | empty |
| M8 | ✅ Covered | defaults verified |
| M9 | ✅ Covered | two regions |
| S1 | ✅ Covered | empty |
| S2 | ✅ Covered | bounds |
| S3 | ✅ Covered | case-insensitive |
| C1 | ✅ Covered | bound |
| C2 | ✅ Covered | deterministic |
| Old tests | ✅ Covered | removed |

In-scope cases: 14 test cases (M1–M9, S1–S3, C1–C2), all ✅.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Complexity = Shannon entropy in bits/residue (the man-page "bits, max log₂20" form of eq. 3) | implementation, M1–M9 |
| 2 | Sequence shorter than window → empty result (spec does not define a value) | S1 |

---

## 7. Open Questions / Decisions

1. **Decision:** The pre-existing `FindLowComplexityRegions` used an invented "dominant single-AA frequency ≥ 0.4" rule with no source basis. Per the implementation-conformance policy this is a defect; the method was rewritten to the SEG bits/residue complexity measure, and the return type changed from `(Start, End, DominantAa, Frequency)` to `(Start, End, Complexity)`. The MCP wrapper and old tests were updated accordingly (in-scope correctness fix).
2. **Decision:** Among the two interconvertible forms of eq. (3), the Shannon-entropy bits/residue form is used because its units and range are stated verbatim in the official SEG man page; recorded as Assumption 1.
