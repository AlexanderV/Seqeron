# Test Specification: REP-TANDEM-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | REP-TANDEM-001 |
| **Area** | Repeats |
| **Title** | Tandem Repeat Detection |
| **Status** | ☑ Complete |
| **Created** | 2026-01-22 |
| **Last Updated** | 2026-01-22 |

---

## Methods Under Test

| Method | Class | Type | Test Priority |
|--------|-------|------|---------------|
| `FindTandemRepeats(seq, minUnit, maxUnit, minReps)` | GenomicAnalyzer | Canonical | Deep testing |
| `GetTandemRepeatSummary(seq, minRepeats)` | RepeatFinder | Summary/Delegate | Smoke testing |

---

## Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| [Wikipedia - Tandem repeat](https://en.wikipedia.org/wiki/Tandem_repeat) | Definition | Adjacent repeating patterns, 8% of human genome |
| [Wikipedia - Microsatellite](https://en.wikipedia.org/wiki/Microsatellite) | Classification | STR = 1-6bp, mutation via slippage, forensic/disease applications |
| Richard et al. (2008) | Review | Comprehensive genomics of DNA repeats |

---

## Test Categories

### MUST Tests (Required for DoD)

All MUST tests are justified by evidence or explicitly marked.

| ID | Test Name | Rationale | Evidence |
|----|-----------|-----------|----------|
| M1 | SimpleTrinucleotide_FindsRepeat | Core algorithm verification - detect ATG×3 | Wikipedia definition |
| M2 | DinucleotideRepeat_FindsRepeat | Common STR type (forensic marker) | Wikipedia - microsatellite forensics |
| M3 | MononucleotideRepeat_FindsRepeat | Homopolymer runs are valid tandems | Wikipedia - microsatellite types |
| M4 | TetranucleotideRepeat_FindsRepeat | Forensic STRs use tetra/penta repeats | Wikipedia - forensic fingerprinting |
| M5 | NoRepeatsFound_ReturnsEmpty | Edge case - sequence without tandems | Standard edge case |
| M6 | EmptySequence_ReturnsEmpty | Boundary - empty input | Standard boundary |
| M7 | MinRepetitionsFilter_RespectsThreshold | Parameter validation | Algorithm specification |
| M8 | MinUnitLengthFilter_RespectsThreshold | Parameter validation | Algorithm specification |
| M9 | PositionCorrect_ZeroBased | Verify position accuracy | Implementation contract |
| M10 | RepetitionCount_Accurate | Count verification | Core correctness |
| M11 | TotalLength_InvariantHolds | TotalLength = Unit.Length × Repetitions | Invariant |
| M12 | FullSequence_Reconstructable | FullSequence matches actual occurrence | Invariant |
| M13 | CAGExpansion_HuntingtonsPattern | Disease-relevant trinucleotide | Wikipedia - trinucleotide disorders |

### SHOULD Tests (Important but not blocking)

| ID | Test Name | Rationale | Evidence |
|----|-----------|-----------|----------|
| S1 | LongRepeat_HandlesCorrectly | Sequences with many repeats | Robustness |
| S2 | EntireSequenceIsRepeat_SingleResult | Edge case - whole sequence | Standard edge case |
| S3 | AdjacentDifferentRepeats_FindsBoth | Multiple distinct patterns | Common scenario |
| S4 | MaxUnitLength_Boundary | Parameter boundary | Algorithm limits |
| S5 | CaseSensitivity_UpperCase | Verify case handling | DnaSequence normalizes |

### COULD Tests (Nice to have)

| ID | Test Name | Rationale |
|----|-----------|-----------|
| C1 | PerformanceBaseline_MediumSequence | Track O(n²) performance |
| C2 | ForensicSTR_RealWorld | Real forensic marker patterns |

---

## Summary Tests (GetTandemRepeatSummary)

These tests verify the delegate method which wraps FindMicrosatellites.

| ID | Test Name | Rationale |
|----|-----------|-----------|
| D1 | MixedRepeats_CorrectAggregation | Summary statistics accurate |
| D2 | NoRepeats_ZeroValues | Edge case |
| D3 | LongestRepeat_Identified | Correct identification |
| D4 | MononucleotideCount_Correct | Category counting |

---

## Existing Test Audit

### GenomicAnalyzerTests.cs (2 tests)

| Test | Classification | Action |
|------|----------------|--------|
| `FindTandemRepeats_SimpleTandem_FindsIt` | Weak | Migrate and strengthen |
| `FindTandemRepeats_DiNucleotide_FindsIt` | Weak | Migrate and strengthen |

### RepeatFinderTests.cs (4 tests)

| Test | Classification | Action |
|------|----------------|--------|
| `GetTandemRepeatSummary_MixedRepeats_CorrectSummary` | Covered | Keep |
| `GetTandemRepeatSummary_NoRepeats_ZeroSummary` | Covered | Keep |
| `GetTandemRepeatSummary_MononucleotideCount_Correct` | Covered | Keep |
| `GetTandemRepeatSummary_LongestRepeat_Identified` | Covered | Keep |

### RepeatFinder_Microsatellite_Tests.cs (3 tests)

| Test | Classification | Action |
|------|----------------|--------|
| `GetTandemRepeatSummary_WithMixedRepeats_CorrectStatistics` | Duplicate | Remove |
| `GetTandemRepeatSummary_LongestRepeat_CorrectlyIdentified` | Duplicate | Remove |
| `GetTandemRepeatSummary_NoRepeats_ZeroValues` | Duplicate | Remove |

---

## Consolidation Plan

1. **Canonical test file:** `GenomicAnalyzer_TandemRepeat_Tests.cs`
   - Comprehensive tests for `FindTandemRepeats`
   - All MUST tests (M1-M13)
   - Selected SHOULD tests (S1-S5)

2. **Summary tests:** Keep in `RepeatFinderTests.cs`
   - `GetTandemRepeatSummary` tests (D1-D4)
   - Already properly located

3. **Remove duplicates from:** `RepeatFinder_Microsatellite_Tests.cs`
   - Remove 3 duplicate GetTandemRepeatSummary tests
   - This file focuses on microsatellite detection (REP-STR-001)

---

## Test Counts

| Category | Count |
|----------|-------|
| MUST | 13 |
| SHOULD | 5 |
| COULD | 2 |
| Summary (delegate) | 4 |
| **Total** | 24 |

### ASSUMPTION Count: 0

All tests are justified by:
- Wikipedia sources (Tandem repeat, Microsatellite)
- Standard edge case coverage
- Algorithm specification/invariants

---

## Open Questions / Decisions

None. Algorithm behavior is well-defined by sources and implementation.
