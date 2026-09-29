# Test Specification: REP-TANDEM-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | REP-TANDEM-001 |
| **Area** | Repeats |
| **Title** | Tandem Repeat Detection |
| **Status** | ☑ Complete |
| **Created** | 2026-01-22 |
| **Last Updated** | 2026-09-29 |

---

## Methods Under Test

| Method | Class | Type | Test Priority |
|--------|-------|------|---------------|
| `FindTandemRepeats(seq, minUnitLength, minRepetitions)` | GenomicAnalyzer | Canonical | Deep testing |
| `GetTandemRepeatSummary(seq, minRepeats)` | RepeatFinder | Summary/Delegate | Smoke testing |

---

## Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| [Wikipedia - Tandem repeat](https://en.wikipedia.org/wiki/Tandem_repeat) | Definition | Adjacent repeating patterns; 8% of human genome; >50 diseases; detection via suffix trees/arrays |
| [Wikipedia - Microsatellite](https://en.wikipedia.org/wiki/Microsatellite) | Classification | STR = 1–6 bp (up to 10 bp by some authors); mutation via slippage (~1 per 1,000 generations); forensic STRs are tetra-/pentanucleotide only |
| Richard et al. (2008) | Review | Comparative genomics of DNA repeats in eukaryotes, MMBR 72(4):686–727 |
| MISA `misa.pl` v1.0 (Thiel et al. 2003, TAG 106:411) — source opened (raw GitHub mirror, 2026-09-29) | Reference tool | `.statistics`: total SSRs; "Distribution to different repeat type classes" = one count per unit size 1–6; "Frequency of identified SSR motifs" counts raw motifs (a second table groups rotations + reverse complement) |
| Krait `src/statistics.py` (Du et al. 2018, Bioinformatics 34:681) — source opened (raw GitHub, 2026-09-29) | Reference tool | Type table Mono…Hexa; "Length (bp)" = `SUM(length)`; density = length / valid (ACGT) size |

---

## Test Categories

### MUST Tests (Required for DoD)

All MUST tests are justified by evidence.

| ID | Test Name | Rationale | Evidence |
|----|-----------|-----------|----------|
| M1 | SimpleTrinucleotide_FindsRepeat | Core algorithm verification — detect ATG×3 | Wikipedia (Tandem repeat): definition |
| M2 | DinucleotideRepeat_FindsRepeat | Most common microsatellite type in human genome (50,000–100,000 loci) | Wikipedia (Microsatellite): structures, locations |
| M3 | MononucleotideRepeat_FindsRepeat | Homopolymer runs are valid tandems (1 bp unit) | Wikipedia (Microsatellite): 1–6 bp classification |
| M4 | TetranucleotideRepeat_FindsRepeat | Forensic STRs are tetra-/pentanucleotide repeats | Wikipedia (Microsatellite): forensic fingerprinting |
| M5 | NoRepeatsFound_ReturnsEmpty | Edge case — sequence without tandems | Standard edge case |
| M6 | EmptySequence_ReturnsEmpty | Boundary — empty input | Standard boundary |
| M7 | MinRepetitionsFilter_RespectsThreshold | Parameter validation | Algorithm specification |
| M8 | MinUnitLengthFilter_RespectsThreshold | Parameter validation | Algorithm specification |
| M9 | PositionCorrect_ZeroBased | Verify position accuracy | Implementation contract |
| M10 | RepetitionCount_Accurate | Count verification | Core correctness |
| M11 | TotalLength_InvariantHolds | TotalLength = Unit.Length × Repetitions | Invariant |
| M12 | FullSequence_Reconstructable | FullSequence matches actual occurrence | Invariant |
| M13 | PentanucleotideRepeat_ForensicSTR | Pentanucleotide repeat — Penta E forensic locus (AAAGA) | Wikipedia (Microsatellite): forensic STRs are tetra-/pentanucleotide |

### SHOULD Tests (Important but not blocking)

| ID | Test Name | Rationale | Evidence |
|----|-----------|-----------|----------|
| S1 | LongRepeat_HandlesCorrectly | Sequences with many repeats | Robustness |
| S2 | EntireSequenceIsRepeat_CoversFullLength | Edge case — whole sequence is one repeat | Standard edge case |
| S3 | AdjacentDifferentRepeats_FindsBoth | Multiple distinct patterns | Common scenario |
| S4 | HexanucleotideRepeat_Boundary | Hexanucleotide (6 bp) upper boundary of microsatellite range | Wikipedia (Microsatellite): 1–6 bp classification |
| S5 | CaseSensitivity_UpperCase | Verify case handling | DnaSequence normalizes to uppercase |

### COULD Tests (Nice to have)

| ID | Test Name | Rationale | Status |
|----|-----------|-----------|--------|
| C1 | PerformanceBaseline_MediumSequence | O(n²) brute-force completes within 30 s for 2 kb sequence | ☑ Implemented |
| C2 | TelomereRepeat_TTAGGG | Vertebrate telomere repeat motif TTAGGG × 4 | ☑ Implemented |

---

## Property Tests (Invariants)

| ID | Test Name | Rationale |
|----|-----------|----------|
| P1 | AllResults_SatisfyMinRepetitions | Every result has Repetitions ≥ minRepetitions |
| P2 | AllResults_SatisfyMinUnitLength | Every result has Unit.Length ≥ minUnitLength |
| P3 | AllResults_WithinSequenceBounds | Position + TotalLength ≤ sequence length |

---

## Summary Tests (GetTandemRepeatSummary)

These tests verify the delegate method which wraps FindMicrosatellites.

| ID | Test Name | Rationale |
|----|-----------|-----------|
| D1 | MixedRepeats_CorrectAggregation | Summary statistics accurate |
| D2 | NoRepeats_ZeroValues | Edge case |
| D3 | LongestRepeat_Identified | Correct identification |
| D4 | MononucleotideCount_Correct | Category counting |
| D5 | PentaAndHexa_CountedAndSumToTotal | All six classes counted (MISA/Krait); counts sum to TotalRepeats; sum-of-lengths 71 vs union coverage 44/46 |
| D6 | OverlappingRuns_SumExceedsCoveredBases | `AAAAATATATAT`: 13 repeat bases, 100 % coverage |
| D7 | NoRepeats / EmptySequence → LongestRepeat and MostFrequentUnit null | Null contract (was a default Position-0 record) |
| D8 | HigherMinRepeats_PartialCoverage | minRepeats 4 → 14/32 = 43.75 % |
| D9 | InvalidArguments_Throw | null → ArgumentNullException; minRepeats < 2 → ArgumentOutOfRangeException |

Summary expected values come from an independent Python reference (brute-force maximal primitive runs
+ aggregation); the code agreed on 9000/9000 random sequences (2026-09-29). Per-class totals of the
same reference agreed exactly with running `misa.pl` for minRepeats ≥ 4; for minRepeats 2–3 MISA differs
only by the run-level conventions documented under REP-STR-001 (greedy consumption of non-primitive
regions, same-period overlapping runs).

---

## Test Counts

| Category | Count |
|----------|-------|
| MUST | 13 |
| SHOULD | 5 |
| COULD | 2 |
| Property (invariants) | 3 |
| Summary (delegate) | 9 |
| **Total** | 32 |

### Deviations and Assumptions

None. All test data and evidence verified against external sources (Wikipedia, accessed 2026-03-01).
