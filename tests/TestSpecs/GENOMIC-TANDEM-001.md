# Test Specification: GENOMIC-TANDEM-001

**Test Unit ID:** GENOMIC-TANDEM-001
**Area:** Analysis
**Algorithm:** Tandem Repeat Detection (`GenomicAnalyzer.FindTandemRepeats`)
**Status:** ☑ Complete (consolidated into REP-TANDEM-001 — see §7)
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Benson, G. (1999) *Tandem Repeats Finder*, *Nucleic Acids Research* 27(2):573–580 | 1 | https://doi.org/10.1093/nar/27.2.573 | 2026-06-14 |
| 2 | Wikipedia "Tandem repeat" (cited primaries) | 4 | https://en.wikipedia.org/wiki/Tandem_repeat | 2026-06-14 |

### 1.2 Key Evidence Points

1. A tandem repeat is "two or more contiguous … copies of a pattern of nucleotides" (k ≥ 2; copies directly adjacent) — Benson 1999 [1].
2. Worked example: `ATTCG ATTCG ATTCG` — unit `ATTCG` (period 5) repeated 3 times, contiguous — Wikipedia [2].
3. Period = unit length for an exact repeat; copy number = number of copies reported per repeat — Benson 1999 [1].
4. Detection of tandem repeats in strings can use suffix trees/arrays; exact-match detection over a single text is the relevant case — Wikipedia [2].

### 1.3 Documented Corner Cases

- Fewer than two contiguous copies is not a tandem repeat (k ≥ 2 floor) — Benson 1999 [1].
- The same region can satisfy several period sizes; each interpretation meeting the threshold is a valid tandem repeat under the definition.

### 1.4 Known Failure Modes / Pitfalls

1. Benson's TRF detects *approximate* copies; the repository detector detects only *exact* contiguous copies (documented simplification). Output agrees with the formal definition over exact repeats — Benson 1999 [1].

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindTandemRepeats(sequence, minUnitLength, minRepetitions)` | GenomicAnalyzer | **Canonical** | Identical to the canonical method of REP-TANDEM-001. |

> **This method is already implemented and tested under REP-TANDEM-001.** GENOMIC-TANDEM-001 is a duplicate Registry entry for the same method of the same class (see §7).

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every reported repeat has Repetitions ≥ `minRepetitions` ≥ 2 (k ≥ 2) | Yes | Benson 1999 "two or more contiguous copies" [1] |
| INV-2 | `TotalLength = Unit.Length × Repetitions` | Yes | period × copy number = repeat length [1] |
| INV-3 | `Position + Unit.Length × Repetitions ≤ sequence.Length` (within bounds) | Yes | copies are contiguous and inside the sequence [1][2] |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Exact trinucleotide unit | `ATGATGATG`, minUnit 3, minReps 3 | `ATG`, pos 0, 3 copies | Benson/Wikipedia definition [1][2] |
| M2 | Worked example | `ATTCGATTCGATTCG`, minUnit 5, minReps 2 | `ATTCG`, pos 0, 3 copies, total length 15 | Wikipedia worked example [2] |
| M3 | minRepetitions floor | sequence with only 2 copies, minReps 3 | not reported | k ≥ 2 / threshold [1] |
| M4 | No repeat / empty | `ACGT` / `""` | empty result | definition (no contiguous copies) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Position is 0-based | tandem starting at index 0 | Position == 0 | contract |
| S2 | TotalLength invariant | any reported repeat | Unit.Length × Repetitions | INV-2 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Telomere motif | `TTAGGG` ×4 | unit `TTAGGG`, 4 copies | vertebrate telomere repeat |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- The canonical fixture for `GenomicAnalyzer.FindTandemRepeats` already exists at
  `tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzer_TandemRepeat_Tests.cs`, delivered and
  committed under **REP-TANDEM-001** (27 tests: 13 MUST, 5 SHOULD, 2 COULD, 3 property, 4 summary).
- That fixture covers M1 (`SimpleTrinucleotide`), the worked-example equivalents (M2 telomere
  `TTAGGG`, tetranucleotide, pentanucleotide forensic), the `minRepetitions`/`minUnitLength`
  thresholds (M3), no-repeat and empty inputs (M4), 0-based position (S1), and the TotalLength
  invariant (S2) plus property tests P1–P3.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 exact trinucleotide | ✅ Covered | `FindTandemRepeats_SimpleTrinucleotide_FindsRepeat` |
| M2 worked example | ✅ Covered | `..._TelomereRepeat_TTAGGG` / tetra-/pentanucleotide tests (same unit≥2 contiguous form) |
| M3 minRepetitions floor | ✅ Covered | `..._MinRepetitionsFilter_RespectsThreshold` |
| M4 no-repeat / empty | ✅ Covered | `..._NoRepeatsFound_ReturnsEmpty`, `..._EmptySequence_ReturnsEmpty` |
| S1 position 0-based | ✅ Covered | `..._PositionCorrect_ZeroBased` |
| S2 TotalLength invariant | ✅ Covered | `..._TotalLength_InvariantHolds` |
| C1 telomere motif | ✅ Covered | `..._TelomereRepeat_TTAGGG` |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzer_TandemRepeat_Tests.cs` — already the single canonical fixture for `FindTandemRepeats` (owned by REP-TANDEM-001).
- **Remove:** nothing. **Do NOT create** a second test file for the same method — that would violate the duplicate-elimination rule (one canonical file per method). GENOMIC-TANDEM-001 reuses the existing fixture.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `GenomicAnalyzer_TandemRepeat_Tests.cs` | Canonical fixture for `FindTandemRepeats` (shared with REP-TANDEM-001) | 27 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | (all M/S/C) | ✅ Covered | None — already covered by the canonical fixture; no new/duplicate tests created | ✅ Done |

**Total items:** 1
**✅ Done:** 1 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | Covered by canonical fixture |
| M2 | ✅ | Covered by canonical fixture |
| M3 | ✅ | Covered by canonical fixture |
| M4 | ✅ | Covered by canonical fixture |
| S1 | ✅ | Covered by canonical fixture |
| S2 | ✅ | Covered by canonical fixture |
| C1 | ✅ | Covered by canonical fixture |

All 7 in-scope cases ✅ (7 of 7). No ❌/⚠ remain.

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Detector reports exact contiguous copies only (not Benson's approximate copies); output matches the formal definition over exact repeats. | §1.4, algorithm doc §5.3 |

---

## 7. Open Questions / Decisions

1. **DECISION — duplicate Registry entry.** The Processing Registry contains two entries for the identical method on the same class: **REP-TANDEM-001** ("Tandem Repeat Detection", Repeats, ☑ Complete) and **GENOMIC-TANDEM-001** ("Tandem Repeat Detection", Analysis). Both list `GenomicAnalyzer.FindTandemRepeats(...)` as canonical with the same O(n²) complexity, and the method-index table maps `FindTandemRepeats` to REP-TANDEM-001 (line ~4630) and `FindTandemRepeats (GenomicAnalyzer)` to GENOMIC-TANDEM-001 (line ~4833) — the same method twice. REP-TANDEM-001 already ships the implementation, the canonical test fixture (`GenomicAnalyzer_TandemRepeat_Tests.cs`), Evidence and an algorithm doc (`docs/algorithms/Repeat_Analysis/Tandem_Repeat_Detection.md`). Per the prompt's duplicate-elimination rule ("one canonical test file per unit; no duplicate tests remain") and workflow-control rule ("note the conflict in the TestSpec and update the checklist entry"), GENOMIC-TANDEM-001 is **resolved by consolidation**: no new production code and no duplicate test file are created; this unit reuses the existing implementation and canonical fixture. Evidence (Benson 1999 + Wikipedia "Tandem repeat") was independently re-retrieved this session to confirm the behavior is correct and source-backed — the formal definition (two or more contiguous copies, k ≥ 2; worked example `ATTCG`×3) matches the implementation.
