# Test Specification: REP-TANDEM-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | REP-TANDEM-001 |
| **Area** | Repeats |
| **Title** | Tandem Repeat Detection |
| **Status** | ☑ Complete |
| **Created** | 2026-01-22 |
| **Last Updated** | 2026-09-30 |

---

## Methods Under Test

| Method | Class | Type | Test Priority |
|--------|-------|------|---------------|
| `FindTandemRepeats(seq, minUnitLength, minRepetitions)` | GenomicAnalyzer | Canonical | Deep testing |
| `GetTandemRepeatSummary(seq, minRepeats)` | RepeatFinder | Summary/Delegate | Smoke testing |
| `GetTandemRepeatSummary(seq, IReadOnlyDictionary<int,int> minRepeatsByUnitLength)` | RepeatFinder | MISA per-unit-size thresholds | Deep (D10) |
| `GetTandemRepeatSummary(string, …)` (uniform / map / map + `MicrosatelliteScanMode`) | RepeatFinder | N/IUPAC-tolerant raw-string overloads (WP11) | Deep (D13) |
| `GetCanonicalMotifClass`, `GetCanonicalMotifFrequencies` | RepeatFinder | MISA repeat-type classes | Deep (D11) |
| `GetStandardMotif(motif, level)`, `GetStandardMotifFrequencies` | RepeatFinder | Krait standard motifs | Deep (D12) |

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
| D10 | `GetTandemRepeatSummary_MisaThresholds_StatSequence` (+ map validation: unit length outside 1–6 throws) | misa.pl default ini: 7 SSRs, class counts 2 / 4 / 1 |
| D11 | `GetCanonicalMotifClass_MatchesMisaStatisticsRowName` (10), `GetCanonicalMotifFrequencies_StatSequence_MatchesMisaClassifiedTable` | misa.pl `.statistics` "Frequency of classified repeat types (considering sequence complementary)": AC/CA/GT/TG → AC/GT; table A/T 2, AC/GT 4, ACAT/ATGT 1 |
| D12 | `GetStandardMotif_MatchesKraitAllLevels` (8 motifs × levels 0–4) | Krait `motif.py` `StandardMotif(level).standard()` (A < T < C < G order; level 2 = rotations + reverse complement, e.g. ACAT → ATAC; Krait GUI default level 3) |

| D13 | `RepeatFinder_TandemSummaryString_Tests` (6) | Raw-string overloads (WP11): N/IUPAC never form or extend an SSR (MISA `[acgt]`); `PercentageOfSequence` denominator = full length incl. N (misa.pl `length $seq`); `MisaRegex` counts = misa.pl `.statistics` (93-mer: 6 SSRs, 3/1/1/1; 87-mer: default 3 SSRs 1/2, uniform 3 copies 4 SSRs 1/3); ACGT-only input = DnaSequence overloads; maximal-runs class counts = sum over ACGT pieces; null/empty → empty summary |

2026-10-01 (WP11): string overloads with `MicrosatelliteScanMode.MisaRegex` vs real `perl misa.pl` `.statistics` on 6 000 random
sequences (5 684 containing N/IUPAC/lower case, 1 882 307 bp) × 6 `misa.ini` definitions: total size, total SSRs, SSR-containing
sequences, sequences with > 1 SSR and the six per-unit-size counts identical (aggregate, all 6 definitions); per sequence (one
misa.pl run each) for the default and the `1-3 2-2 3-2 4-2 5-2 6-2` definitions: 12 000 / 12 000 identical.

Summary expected values come from an independent Python reference (brute-force maximal primitive runs
+ aggregation); the code agreed on 9000/9000 random sequences (2026-09-29). Per-class totals of the
same reference agreed exactly with running `misa.pl` for minRepeats ≥ 4; for minRepeats 2–3 MISA differs
only by the run-level conventions documented under REP-STR-001 (greedy consumption of non-primitive
regions, same-period overlapping runs).

2026-09-30: MISA class names equal real misa.pl `.statistics` rows for all 5 356 primitive motifs of 1–6 bp
(one misa.pl run per motif); Krait standard motifs equal Krait's `StandardMotif.standard()` for all 5 460 motifs
of 1–6 bp at levels 0–4 (note: Krait caches in a class-level dict shared by all levels; the reference was run with
a fresh cache per level). The classified table built from misa.pl's own SSR list equals misa.pl's table in all
6 `misa.ini` configurations (6 048 sequences); built from this library's SSR list it differs only through the
documented SSR-list conventions. Tests: `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_MisaCompound_Tests.cs`.

---

## Test Counts

| Category | Count |
|----------|-------|
| MUST | 13 |
| SHOULD | 5 |
| COULD | 2 |
| Property (invariants) | 3 |
| Summary (delegate) | 13 |
| **Total** | 36 |

### Deviations and Assumptions

None. All test data and evidence verified against external sources (Wikipedia, accessed 2026-03-01).
