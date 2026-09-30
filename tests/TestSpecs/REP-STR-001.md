# Test Specification: REP-STR-001

**Test Unit:** Microsatellite Detection (STR) — perfect (default) + approximate / imperfect / interrupted (opt-in, TRF model)
**Canonical Class:** `RepeatFinder`
**Primary Method:** `FindMicrosatellites` (perfect) + `FindApproximateTandemRepeats` (approximate, Benson 1999)
**Status:** Complete
**Last Updated:** 2026-09-30

> §1–§8 below cover the perfect-STR detector. **§9 adds the opt-in approximate (TRF) detector**
> (`FindApproximateTandemRepeats`), added for the REP-STR-001 limitation fix; its evidence is in
> `docs/Evidence/REP-STR-001-Evidence.md`.

---

## 1. Scope

### Methods Under Test

| Method | Class | Type | Test Depth |
|--------|-------|------|------------|
| `FindMicrosatellites(DnaSequence, int, int, int)` | RepeatFinder | Canonical | Deep |
| `FindMicrosatellites(string, int, int, int)` | RepeatFinder | Overload | Deep |
| `FindMicrosatellites(DnaSequence, ..., CancellationToken, IProgress<double>)` | RepeatFinder | Cancellable + progress | Deep (C02/C03) |
| `FindMicrosatellites(string, ..., CancellationToken, IProgress<double>)` | RepeatFinder | Cancellable + progress | Deep (C02/C03) |
| `FindMicrosatellites(DnaSequence \| string, IReadOnlyDictionary<int,int>, CancellationToken, IProgress<double>)` | RepeatFinder | MISA per-unit-size thresholds | Deep (§11) |
| `MisaDefaultMinRepeats` / `MisaDefaultMaxInterruption` | RepeatFinder | MISA `misa.ini` defaults | Deep (§11) |
| `FindCompoundMicrosatellites(...)`, `AssembleCompoundMicrosatellites(string, IEnumerable<MicrosatelliteResult>, int)` | RepeatFinder | MISA compound SSRs (c / c*) | Deep (§11) |

### Supporting Methods
| Method | Class | Test Approach |
|--------|-------|---------------|
| `GetTandemRepeatSummary` | RepeatFinder | Integration via microsatellite |

---

## 2. Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| Wikipedia: Microsatellite | Encyclopedia | Definition: 1-6 bp motifs repeated 5-50 times; terminology (mono-/di-/tri-nucleotide); forensic use of tetra-/penta-nucleotide |
| Wikipedia: Trinucleotide repeat disorder | Encyclopedia | CAG repeat thresholds: HD normal 6-35, pathogenic 36-250; disease-specific repeat ranges |
| Richard GF et al. (2008) MMBR | Peer-reviewed | Comprehensive review of repeat dynamics |
| Tóth G et al. (2000) Genome Res | Peer-reviewed | Microsatellite distribution analysis |
| Thiel T et al. (2003) TAG 106:411 — MISA `misa.pl` source | Reference tool | Per-motif-size leftmost regex `([acgt]{p})\2{k-1,}`; reject non-primitive ("false type") motifs; ACGT-only |
| Du L et al. (2018) Bioinformatics 34:681 — Krait / pytrf 1.5.0 `str.c` | Reference tool | Run seeded at its start, `repeat = length / p` (complete copies), `N` skipped |
| Kolpakov & Kucherov (1999) FOCS | Peer-reviewed | Maximal repetition (run) definition: left/right-maximal, minimal period |

---

## 3. Test Categories

### 3.1 MUST Tests (Evidence-Based)

| ID | Test Case | Evidence | Rationale |
|----|-----------|----------|-----------|
| M01 | Mononucleotide repeat detection (A×n) | Wikipedia: "TATATATATA is a dinucleotide microsatellite" - confirms mono/di classification | Verify basic mono detection |
| M02 | Dinucleotide repeat detection (CA×n) | Wikipedia: AC/CA common in human genome | Common biological motif |
| M03 | Trinucleotide CAG repeat detection | Wikipedia: Huntington's disease CAG repeats | Critical medical relevance |
| M04 | Tetranucleotide GATA detection | Wikipedia: forensic markers use tetra-nucleotide | Forensic application |
| M05 | Empty sequence returns empty | Standard edge case | Defensive programming |
| M06 | minRepeats filter respected | Algorithm specification | Contract validation |
| M07 | RepeatType classification correct | Wikipedia terminology | Type verification |
| M08 | FullSequence property correctness | Invariant: unit × count | Mathematical correctness |
| M09 | Position accuracy | Invariant: position in range | Location correctness |
| M10 | Null DnaSequence throws | Defensive programming | Error handling |
| M11 | Invalid parameters throw (DnaSequence + string) | API contract | Error handling |
| M12 | Redundant unit filtering | Implementation: "ATAT" → "AT"×2 | Avoid duplicate reporting |
| M13 | Wikipedia TATA example: "TATATATATA" | Wikipedia: exact dinucleotide example (Structures section) | Source fidelity |
| M14 | Wikipedia GTC example: "GTCGTCGTCGTCGTC" | Wikipedia: exact trinucleotide example (Structures section) | Source fidelity |
| M15 | GGAT myoglobin repeat | Wikipedia: first microsatellite characterized (Weller et al. 1984) | Historical reference |
| M16 | String overload parameter validation parity | API contract: string overloads must validate identically to DnaSequence | Eliminates silent invalid-parameter bugs |

### 3.2 SHOULD Tests (Quality/Coverage)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| S01 | Multiple different repeats in sequence | Real-world: genomes contain many STRs |
| S02 | String overload parity with DnaSequence | API consistency |
| S03 | Hexanucleotide repeat detection | Complete unit length coverage |
| S04 | Case insensitivity | Robust input handling |
| S05 | Non-standard characters (N) | IUPAC ambiguity handling |
| S06 | Adjacent different repeat types | Complex pattern handling |
| S07 | TandemRepeatSummary accuracy | Summary statistics |

### 3.3 COULD Tests (Extended)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| C01 | Large sequence performance | Scalability |
| C02 | Cancellation mid-operation | Async operation support |
| C03 | Progress reporting (`IProgress<double>` on the cancellable overloads) | User feedback |

---

## 4. Test Data

### Biological Test Cases (Evidence-Based)

| Name | Sequence | Expected Result | Source |
|------|----------|-----------------|--------|
| Huntington CAG | `ATGCAGCAGCAGCAGCAGTGA` | GCA×5 at position 2 (run-start phase; = CAG×5 locus) | Wikipedia: HD has CAG repeats |
| Dinucleotide CA | `AAACACACACACAAA` | AC×5 at position 2 (run of 11 bp, 5 complete copies) | Wikipedia: common microsatellite |
| Mononucleotide A | `ACGTAAAAAACGT` | A×6 at position 4 | Basic mononucleotide |
| Tetranucleotide GATA | `AAGATAGATAGATAGATAAA` | GATA-family×4 | Wikipedia: forensic marker |
| EcoRI site as repeat | `GAATTCGAATTCGAATTC` | GAATTC×3 | Hexanucleotide example |
| Wikipedia TA example | `TATATATATA` | TA×5 at position 0 | Wikipedia: "TATATATATA is a dinucleotide microsatellite" |
| Wikipedia GTC example | `GTCGTCGTCGTCGTC` | GTC×5 at position 0 | Wikipedia: "GTCGTCGTCGTCGTC is a trinucleotide microsatellite" |
| GGAT myoglobin | `AAGGATGGATGGATGGATAA` | GGAT-family×4 | Wikipedia: first microsatellite (Weller et al. 1984) |

### Edge Case Data

| Name | Sequence | minRepeats | Expected |
|------|----------|------------|----------|
| Empty | "" | 3 | Empty |
| Too short | "AT" | 3 | Empty |
| Exactly minRepeats | "ATATAT" | 3 | AT×3 |
| Below threshold | "ATAT" | 3 | Empty |

---

## 5. Invariants to Assert

```csharp
// For each result r in FindMicrosatellites(seq, minUnit, maxUnit, minReps):
Assert.That(r.RepeatCount, Is.GreaterThanOrEqualTo(minReps));
Assert.That(r.RepeatUnit.Length, Is.InRange(minUnit, maxUnit));
Assert.That(r.TotalLength, Is.EqualTo(r.RepeatUnit.Length * r.RepeatCount));
Assert.That(r.FullSequence, Is.EqualTo(string.Concat(Enumerable.Repeat(r.RepeatUnit, r.RepeatCount))));
Assert.That(r.Position, Is.InRange(0, seq.Length - r.TotalLength));
Assert.That(seq.Substring(r.Position, r.TotalLength), Is.EqualTo(r.FullSequence));
```

---

## 6. Audit of Existing Tests

### 6.1 Discovery Summary

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_Microsatellite_Tests.cs`
- **Supporting file:** `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinderTests.cs` — TandemRepeatSummary tests (S07)
- **MISA / progress file (2026-09-30):** `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_MisaCompound_Tests.cs` — §11, C02, C03
- **Cross-reference:** `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Core/PerformanceExtensionsTests.cs` — cancellation smoke (separate test unit)

### 6.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| **MUST Tests** | | |
| M01 — Mononucleotide detection | ✅ Covered | Exact: A×6, pos=4, len=6, type=Mononucleotide |
| M02 — Dinucleotide CA detection | ✅ Covered | 2026-09: count=1, AC×5 pos=2 (the CA×5@3 rotation of the same run is no longer reported) |
| M03 — Trinucleotide CAG detection | ✅ Covered | 2026-09: count=1, GCA×5 pos=2 (run starts at the G; CAG×5@3 rotation not re-reported) |
| M04 — Tetranucleotide GATA detection | ✅ Covered | 2026-09: count=1, AGAT×4 pos=1 |
| M05 — Empty sequence returns empty | ✅ Covered | — |
| M06 — minRepeats filter respected | ✅ Covered | 5 tests: MinRepeatsInvariant (×3) + ExactlyMinRepeats + BelowMinRepeats |
| M07 — RepeatType classification | ✅ Covered | Strengthened: 6 TestCases, exact count=1, unit×5, pos=0 |
| M08 — FullSequence / TotalLength invariants | ✅ Covered | 2 invariant tests with loop over all results |
| M09 — Position accuracy | ✅ Covered | 2 tests: PositionInvariant + SequenceAtPosition |
| M10 — Null DnaSequence throws | ✅ Covered | ArgumentNullException |
| M11 — Invalid parameters throw (DnaSequence) | ✅ Covered | 3 tests: minUnit=0, max<min, minRepeats=1 |
| M12 — Redundant unit filtering | ✅ Covered | NEW: AT×4 not ATAT×2 from "ATATATAT" |
| M13 — Wikipedia TATA example | ✅ Covered | Exact: TA×5, pos=0, len=10 |
| M14 — Wikipedia GTC example | ✅ Covered | Exact: GTC×5, pos=0, len=15 |
| M15 — GGAT myoglobin repeat | ✅ Covered | Strengthened: count=1, GGAT×4, pos=2, len=16 |
| M16 — String overload validation parity | ✅ Covered | 3 tests: minUnit=0, max<min, minRepeats=1 |
| **SHOULD Tests** | | |
| S01 — Multiple different repeats | ✅ Covered | Strengthened: exact count=3, all results verified |
| S02 — String overload parity | ✅ Covered | Compares DnaSequence vs string results field-by-field |
| S03 — Hexanucleotide detection | ✅ Covered | GAATTC×3, exact values |
| S04 — Case insensitivity | ✅ Covered | lowercase "cagcagcagcag" → CAG×4 |
| S05 — Non-standard characters (N) | ✅ Covered | DnaSequence rejects N; 2026-09: string overload never reports a unit containing N (MISA `[acgt]`, pytrf skips N) |
| S06 — Adjacent different repeat types | ✅ Covered | NEW: A×5 pos=0 + CAG×3 pos=5 |
| S07 — TandemRepeatSummary accuracy | ✅ Covered | 4 tests in RepeatFinderTests.cs |
| **COULD Tests** | | |
| C01 — Large sequence performance | ❌ Missing | Benchmark-only (SuffixTree.Benchmarks), not unit-testable |
| C02 — Cancellation mid-operation | ✅ Covered | 2026-09-30: pre-cancelled token → `OperationCanceledException` on enumeration (DnaSequence, string, map overloads); token cancelled from the first progress callback → throws at the next check, no further report (`RepeatFinder_MisaCompound_Tests`) |
| C03 — Progress reporting | ✅ Covered | 2026-09-30 (was wrongly listed "not implemented": the cancellable overloads take `IProgress<double>`): values non-decreasing, in [0, 1), final report exactly 1.0; results identical to the non-cancellable overload; short input → only the final 1.0; map overload likewise |
| **Edge Cases** | | |
| SequenceTooShort | ✅ Covered | "AT" with minRepeats=3 → empty |
| EntireSequenceIsRepeat | ✅ Covered | CAG×10, exact: count=1, pos=0, len=30 |

### 6.3 Consolidation Plan

- **Canonical file:** `RepeatFinder_Microsatellite_Tests.cs` — all microsatellite detection tests
- **Remove:** 14 duplicate microsatellite tests from `RepeatFinderTests.cs` (Mono, Di, Tri, Tetra, NoRepeats, Multiple, MinRepeats, StringOverload, Empty, FullSequence, Null, InvalidMin, InvalidMax, InvalidRepeats) — all covered by stronger tests in canonical file
- **Keep:** `RepeatFinderTests.cs` retains TandemRepeatSummary (4 tests) + InvertedRepeats null (1 test)

### 6.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `RepeatFinder_Microsatellite_Tests.cs` | Canonical (REP-STR-001) | 34 |
| `RepeatFinderTests.cs` | TandemRepeatSummary + InvertedRepeats null | 5 |

### 6.5 Work Queue

| # | Test Case ID | §6.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M02 — Dinucleotide | ⚠ Weak | Strengthened: exact count, positions, all fields | ✅ Done |
| 2 | M03 — Trinucleotide | ⚠ Weak | Strengthened: exact count=2, both rotations verified | ✅ Done |
| 3 | M04 — Tetranucleotide | ⚠ Weak | Strengthened: exact count=2, AGAT+GATA verified | ✅ Done |
| 4 | M07 — RepeatType (×6) | ⚠ Weak | Strengthened: exact count=1, all fields verified | ✅ Done |
| 5 | M12 — Redundant unit | ❌ Missing | Implemented: AT×4 not ATAT×2 | ✅ Done |
| 6 | M15 — GGAT myoglobin | ⚠ Weak | Strengthened: exact count=1, GGAT×4, pos=2 | ✅ Done |
| 7 | S01 — Multiple repeats | ⚠ Weak | Strengthened: exact count=3, all results verified | ✅ Done |
| 8 | S05 — Non-standard chars | ❌ Missing | Implemented: DnaSequence rejects, string overload accepts | ✅ Done |
| 9 | S06 — Adjacent types | ❌ Missing | Implemented: A×5 + CAG×3 exact values | ✅ Done |
| 10 | 🔁 RepeatFinderTests dups | 🔁 Duplicate | Removed 14 duplicate tests | ✅ Done |

**Total items:** 10
**✅ Done:** 10 | **⛔ Blocked:** 0 | **Remaining:** 0

### 6.6 Post-Implementation Coverage

All MUST (M01-M16) and SHOULD (S01-S07) tests are ✅ Covered.
COULD tests: C02 and C03 covered (2026-09-30); C01 (performance) is benchmark-only.

---

### 6.7 Review 2026-09 (maximal-run semantics)

| ID | Test | Status |
|----|------|--------|
| M17 | `FindMicrosatellites_MaximalRuns_EachLocusReportedOnce` (6 cases: `ATATATA`, `ATATATAT` k=2, (CAG)×10+CA, `AAAAAACACACAC`, `ACACACGCGCGC`, `AAGATAGATAGATAGATAAA`) — values from brute-force maximal-repetition reference, MISA/pytrf-consistent | ✅ |
| M18 | `FindMicrosatellites_CancellableDnaOverload_InvalidParameters_Throw` | ✅ |

## 7. Open Questions

None — all behavior verified against external sources.

---

## 8. Validation Checklist

- [x] Evidence sources documented
- [x] Must tests defined with evidence
- [x] Existing tests audited
- [x] Coverage Classification completed (§6.2)
- [x] Invariants specified
- [x] Tests implemented
- [x] Zero assumptions — all design decisions backed by external sources
- [x] Duplicates removed (14 from RepeatFinderTests.cs)
- [x] Weak tests strengthened (6 tests hardened with exact values)
- [x] Tests passing (all REP-STR-001 fixtures green, 2026-09-30)
- [ ] Zero warnings (4 pre-existing in ApproximateMatcher_EditDistance_Tests.cs)

---

## 9. Approximate / Imperfect Tandem-Repeat Detection (TRF model — opt-in)

> **Updated 2026-09-30 (REP-APPROX-001, batch B04):** the detector follows the compiled TRF 4.10.0 model
> (wraparound DP, statistics between ADJACENT copies, k-tuple + sum-of-heads detection); the values below are
> the TRF-locked values of the current tests (the 2026-06 hand-derived "vs consensus" values — A2 94.4̄ %,
> A3 reported, A4 9.67 copies / 96.67 % — were wrong and are gone). Full spec: `tests/TestSpecs/REP-APPROX-001.md`;
> evidence: `docs/Evidence/REP-APPROX-001-Evidence.md`.


**Method under test:** `RepeatFinder.FindApproximateTandemRepeats(DnaSequence | string, int minPeriod, int maxPeriod, int minScore)`.

### 9.1 Evidence

| Source | Key Information |
|--------|-----------------|
| Benson G (1999), Nucleic Acids Res 27(2):573–580, https://doi.org/10.1093/nar/27.2.573 | Approximate tandem repeat = "two or more contiguous, *approximate* copies of a pattern". Reported statistics: period size, copy number (copies aligned with consensus), consensus size, percent matches between adjacent copies overall, percent indels between adjacent copies overall, alignment score. Scoring (match, mismatch, gap) = (+2, −7, −7); "Only those repeats scoring at least 50 … are reported". Consensus "by majority rule from the alignment". |
| TRF README / definitions (Benson-Genomics-Lab/TRF; tandem.bu.edu) | Recommended Match/Mismatch/Delta = 2/7/7; Minscore example 50; statistic definitions verbatim. |

### 9.2 Scoring constants (source-traceable)

| Constant | Value | Source |
|----------|-------|--------|
| Match weight | +2 | Benson (1999) |
| Mismatch penalty | −7 | Benson (1999) recommended |
| Indel penalty | −7 / gap column | Benson (1999) recommended (flat) |
| `DefaultApproximateMinScore` | 50 | Benson (1999) |

### 9.3 MUST cases (values from compiled TRF 4.10.0, `trf seq 2 7 7 80 10 <minscore> <maxperiod> -d -h`)

TRF row format: start end period copies consensus-size %matches %indels score (1-based inclusive indices; TRF
truncates percentages to integers — the API returns the exact ratios).

| ID | Sequence (periods, minScore) | Expected (TRF-locked) | Evidence |
|----|-------------------|-----------------------------|----------|
| A1 | `CACACACACA` (1–6, 10) | 1–10, period 2, copies 5.0, size 2, %matches **100**, %indels **0**, score **20**, consensus `CA`, entropy 1.00 | TRF row `1 10 2 5.0 2 100 0 20 50 50 0 0 1.00 CA` |
| A2 | `CAGCAGCAGTAGCAGCAG` (3–3, 10) | 1–18, period 3, copies 6.0, %matches **86.67 (= 1300/15; 13 matches, 2 mismatches between adjacent copies)**, %indels **0**, score **27**, `CAG`, %T = 100/18 | TRF row `1 18 3 6.0 3 86 0 27 33 27 33 5 1.80 CAG` + instrumented TRF counts |
| A2b | (same A2 sequence) perfect detector `FindMicrosatellites(...,1,6,3)` | only `CAG`×3 at pos 0 (fragmented) | contrast: perfect detector breaks the interrupted tract |
| A3 | `CACACATACACA` (2–2 and 1–500, 10) | **nothing reported** — the only run of ≥ 4 matches at distance 2 has 4 heads < sum-of-heads criterion 5 (d = 2, k = 4) | compiled TRF reports nothing |
| A4 | `CAGCAGCAGCAGCAGAGCAGCAGCAGCAG` (3–3, 10; one deletion) | 1–29, period 3, copies **10.0** (30 aligned consensus columns / 3), %matches **92.59 (= 2500/27)**, %indels **7.41 (= 200/27)**, score **51**, `CAG` (25 matches, 0 mismatches, 2 indels between adjacent copies) | TRF row `1 29 3 10.0 3 92 7 51 34 31 34 0 1.58 CAG` |
| A4b | A4 sequence (1–500, 10) | three overlapping rows: 1–29 p 3 score 51; 2–29 p 14 (2.0 copies, 100 %) score 56; 5–26 p 11 (2.0 copies, 100 %) score 44; at minScore 50 only (3, 51) and (14, 56) remain | TRF reports overlapping periods |
| A5 | `CACACACACA` at default minScore 50 | empty (score 20 < 50) | Benson "≥ 50 … reported" |
| A6 | A4 sequence (3–3) at default minScore | reported, score 51; minScore 52 → empty | Benson "≥ 50 … reported" |
| A7 | `""`, `null`, `A`, `ACG` (6–6) | empty | edge |
| A8 | `ACGTGCAT` (no repeat) | empty | edge |
| A9 | `minPeriod = 0` / `maxPeriod < minPeriod` / `maxPeriod > 2000` / `minScore < 1` (also on empty input) | `ArgumentOutOfRangeException` (eager) | parameter validation (TRF: MaxPeriod ≤ 2000, positive scores) |
| A10 | `null` DnaSequence | `ArgumentNullException` | parameter validation |
| A11 | determinism — same input twice / DnaSequence vs string → identical results | equal | determinism |
| A12 | all-`N` input; `N` inside a repeat | never reported; `N` counts as a mismatch | TRF matches only identical A/C/G/T |
| A13 | TRF README test sequences, homopolymer `A29 G A26`, lowercase/flanked input | published TRF tables reproduced (e.g. homopolymer `1 57 1 57.0 1 96 0 105 …`, entropy 0.13) | TRF README / compiled TRF |
| A14 | `TrfSumOfHeadsCriterion(d)` | TRF `sumdata80` table (d = 1 → 5, 29 → 9, 30 → 6, 159 → 69, 160 → 43, 2000 → 818) | TRF 4.10.0 `tr30dat.c` |

### 9.4 Coverage status

All A1–A14 implemented in `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_ApproximateTandemRepeats_Tests.cs`,
exact assertions (`Within(1e-9)`) on percentages / copy number against the TRF rows. ✅ Done = 14 / 14; Remaining = 0.

### 9.5 Residual (declared)

Detection is the TRF k-tuple + sum-of-heads trigger with TRF's own analysis component (WDP, consensus,
adjacent-copy statistics identical to compiled TRF on 1524/1524 candidates with pattern ≤ 20). The exhaustive
(start, period) window scan of the 2026-06 version was removed. TRF's apparent-size (simulated) criterion,
random-walk distance ranges, best-period list for d > 250 and narrow-band WDP for patterns > 20 are not
reproduced; measured agreement and the LIMITATIONS proposal are in `tests/TestSpecs/REP-APPROX-001.md` and
`docs/Validation/review-2026-09/B04.md` (REP-APPROX-001).

## 10. TRF Bernoulli statistical-significance measures (Benson 1999 — opt-in)

**Method under test:** `RepeatFinder.ComputeBernoulliStatistics(string repeatTract, int period, double expectedMatchProbability = 0.80)` → `TandemRepeatBernoulliStatistics`.

### 10.1 Evidence

Benson (1999) NAR 27(2):573–580; TRF desc/definitions pages (verbatim): "We model alignment of two tandem copies of a pattern of length n by a sequence of n independent Bernoulli trials"; PM = P(Heads) = "the average percent identity between the copies"; PI = "the average percentage of insertions and deletions between the copies"; statistics "between adjacent copies … not between the sequence and the consensus pattern"; defaults "PM = .80 and PI = .10".

### 10.2 MUST cases (exact hand-derived values; each trial = one alignment column between adjacent copies)

| ID | Input | Expected | Evidence |
|----|-------|----------|----------|
| B1 | `CACACACACA`, period 2 | 4 pairs, 8 trials, 8/0/0, PM **1.0**, PI **0**, E[matches] **8**, meets-0.80 **true** | PM = average % identity; perfect tract |
| B2 | `CAGCAGCAGTAGCAGCAG`, period 3 | 5 pairs, 15 trials, 13/2/0, PM **13/15**, PI **0**, E[matches] **13** | adjacent-copy PM (≠ 17/18 consensus) |
| B3 | `CACACATACACA`, period 2 | 10 trials, 8/2/0, PM **0.80**, E[matches] **8**, meets-0.80 **true** (inclusive) | PM on the default threshold |
| B4 | `ACACTGTG`, period 4 | the WDP aligns only one copy (`TGTG`): **0** pairs, **0** trials, PM **0**, meets-0.80 **false** | no two copies to compare |
| B5 | `CAGCAGCAGTAGCAGCAG`, period 3, PM-threshold 0.80 / 0.90 | meets 0.80 **true**, meets 0.90 **false** | custom PM threshold |
| B6 | `CAGCAGCAGCAGCAGAGCAGCAGCAGCAG` (deletion), period 3 | 9 pairs, 27 trials, **25/0/2**, PM **25/27**, PI **2/27** | TRF adjacent-copy counts (92 % / 7 %) |
| B7 | `CAGCAGCAGTAGCAGCAG`, period 3 | PM + mismatch-fraction + PI = **1.0** | Bernoulli outcomes partition the trials |
| B8 | `CAG`, period 3 | `ArgumentException` | model of two copies undefined for one copy |
| B9 | `null` / period 0 / period 2001 / PM 1.5 / PM NaN | `ArgumentNullException` / `ArgumentOutOfRangeException` ×4 | parameter validation |
| B10 | exposed defaults | `TrfDefaultMatchProbability` **0.80**, `TrfDefaultIndelProbability` **0.10** | Benson (1999) defaults |

### 10.3 Coverage status

All B1–B10 implemented in `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_ApproximateTandemRepeats_Tests.cs` (region "ComputeBernoulliStatistics"), exact assertions with `.Within(1e-9)` on the probabilities/expected matches. ✅ Done = 10 / 10; Remaining = 0. Invariants INV-09, INV-10.

## 11. MISA per-unit-size thresholds, compound SSRs, canonical motif classes (2026-09-30)

Sources opened: MISA `misa.pl` v1.0 (Thiel et al. 2003; raw GitHub mirror `cfljam/SSR_marker_design`) — `misa.ini`
`definition(unit_size,min_repeats): 1-10 2-6 3-5 4-5 5-5 6-5`, `interruptions(max_difference_for_2_SSRs): 100`,
compound assembly (types `c`, `c*`), `.statistics` table "Frequency of classified repeat types (considering sequence
complementary)"; Krait `src/motif.py` `StandardMotif` (lmdu/krait). Expected values are copied from real
`perl misa.pl` runs and from Krait's `standard()` in Python.

| ID | Test | Expected (reference) |
|----|------|----------|
| T01 | `MisaDefaultMinRepeats`, `MisaDefaultMaxInterruption` | {1:10, 2:6, 3:5, 4:5, 5:5, 6:5}, 100 (misa.ini) |
| T02 | `CCCCCCCCCGCCCCCCCCCC`, MISA thresholds | only `(C)10` at 11–20 (misa.pl `p1`) |
| T03 | 7-SSR test sequence, MISA thresholds | SSR list = misa.pl list (AC, AC, GT, TG, ACAT, A, T with coordinates) |
| T04 | uniform map ≡ `minRepeats` overload | 300 random sequences identical |
| T05 | map with unit lengths {3, 1} only | only those unit lengths, ordered by unit length |
| T06 | map validation (null, empty, unit < 1, copies < 2, summary unit > 6) | eager exceptions |
| T07 | `GetTandemRepeatSummary(dna, MisaDefaultMinRepeats)` | misa.pl class counts 2 / 4 / 1, total 7 |
| T08 | compound rows (6 cases: interrupted `c`, overlapping `c*`, exactly 100 interrupting bases, adjacent, 4-component mixed, adjacent+interrupted+overlapping) | misa.pl type, notation, size, start, end |
| T09 | 101 interrupting bases | two single SSRs (misa.pl `p2`, `p2`); `maxInterruption: 101` joins |
| T10 | `maxInterruption` 0 / 1 | misa.pl with `interruptions 0` / `1` |
| T11 | chaining compares with the previous component's end; compound end = last component's end | misa.pl loop semantics |
| T12 | equal starts keep input order (stable) | MISA ties = Perl hash order |
| T13 | documented detection differences vs MISA (overlap < p; greedy non-primitive consumption) | misa.pl `(AGAAA)8 46-85` vs maximal `(AAGAA)8 45-84`; `(TAAACT)6 131-166` vs `(CTTAAA)7 129-170` |
| T14 | `GetCanonicalMotifClass` (10 motifs) | misa.pl `.statistics` row names (AC/GT, A/T, AT/AT, ACAT/ATGT, …) |
| T15 | `GetStandardMotif` levels 0–4 (8 motifs) | Krait `StandardMotif(level).standard()` (fresh cache) |
| T16 | `GetCanonicalMotifFrequencies` / `GetStandardMotifFrequencies` on the test sequence | misa.pl classified table A/T 2, AC/GT 4, ACAT/ATGT 1; Krait A 2, AC 4, ATAC 1 |

Bulk cross-checks (scratch harness, not unit tests): 6 048 sequences in 6 `misa.ini` configurations (72 974 SSRs):
brute-force maximal primitive runs with per-size thresholds 0 mismatching sequences; `AssembleCompoundMicrosatellites`
fed with misa.pl's own SSR list (in misa.pl order) reproduces every misa.pl row (0 mismatching sequences, 12 330
compounds); with this library's SSR list, compound rows are identical wherever the SSR lists agree, and every SSR-list
difference is one of the two documented conventions (0 unexplained); MISA class names = misa.pl for all 5 356 primitive
motifs of 1–6 bp; Krait standard motif = Krait for all 5 460 motifs of 1–6 bp at levels 0–4.

