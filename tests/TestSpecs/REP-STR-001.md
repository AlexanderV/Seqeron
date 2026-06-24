# Test Specification: REP-STR-001

**Test Unit:** Microsatellite Detection (STR) — perfect (default) + approximate / imperfect / interrupted (opt-in, TRF model)
**Canonical Class:** `RepeatFinder`
**Primary Method:** `FindMicrosatellites` (perfect) + `FindApproximateTandemRepeats` (approximate, Benson 1999)
**Status:** Complete
**Last Updated:** 2026-06-24

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
| `FindMicrosatellites(DnaSequence, ..., CancellationToken)` | RepeatFinder | Cancellable | Smoke |
| `FindMicrosatellites(string, ..., CancellationToken)` | RepeatFinder | Cancellable | Smoke |

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
| C03 | Progress reporting | User feedback |

---

## 4. Test Data

### Biological Test Cases (Evidence-Based)

| Name | Sequence | Expected Result | Source |
|------|----------|-----------------|--------|
| Huntington CAG | `ATGCAGCAGCAGCAGCAGTGA` | CAG×5 at position 3 | Wikipedia: HD has CAG repeats |
| Dinucleotide CA | `AAACACACACACAAA` | CA×6 (or AC×6) | Wikipedia: common microsatellite |
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

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/RepeatFinder_Microsatellite_Tests.cs` — 34 tests
- **Supporting file:** `tests/Seqeron/Seqeron.Genomics.Tests/RepeatFinderTests.cs` — TandemRepeatSummary tests (4 tests, S07)
- **Cross-reference:** `tests/Seqeron/Seqeron.Genomics.Tests/PerformanceExtensionsTests.cs` — cancellation smoke (separate test unit)

### 6.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| **MUST Tests** | | |
| M01 — Mononucleotide detection | ✅ Covered | Exact: A×6, pos=4, len=6, type=Mononucleotide |
| M02 — Dinucleotide CA detection | ✅ Covered | Strengthened: count=2, AC×5 pos=2 + CA×5 pos=3 |
| M03 — Trinucleotide CAG detection | ✅ Covered | Strengthened: count=2, GCA×5 pos=2 + CAG×5 pos=3 |
| M04 — Tetranucleotide GATA detection | ✅ Covered | Strengthened: count=2, AGAT×4 pos=1 + GATA×4 pos=2 |
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
| S05 — Non-standard characters (N) | ✅ Covered | NEW: DnaSequence rejects N; string overload treats as regular char |
| S06 — Adjacent different repeat types | ✅ Covered | NEW: A×5 pos=0 + CAG×3 pos=5 |
| S07 — TandemRepeatSummary accuracy | ✅ Covered | 4 tests in RepeatFinderTests.cs |
| **COULD Tests** | | |
| C01 — Large sequence performance | ❌ Missing | Benchmark-only (SuffixTree.Benchmarks), not unit-testable |
| C02 — Cancellation mid-operation | ✅ Covered | Smoke test (deep cancellation in PerformanceExtensionsTests) |
| C03 — Progress reporting | ❌ Missing | Not implemented in API |
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
COULD tests: C01 (performance) is benchmark-only, C03 (progress) is not implemented. Neither blocks completion.

---

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
- [x] Tests passing (45/45)
- [ ] Zero warnings (4 pre-existing in ApproximateMatcher_EditDistance_Tests.cs)

---

## 9. Approximate / Imperfect Tandem-Repeat Detection (TRF model — opt-in)

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

### 9.3 MUST cases (exact hand-derived values)

| ID | Sequence (period) | Expected (verbatim-derived) | Evidence |
|----|-------------------|-----------------------------|----------|
| A1 | `CACACACACA`, period 2 | consensus `CA`, copies 5.0, %matches **100**, %indels **0**, score **20** | perfect-alignment control |
| A2 | `CAGCAGCAGTAGCAGCAG`, period 3 (one substitution at idx 9) | consensus `CAG`, copies 6.0, %matches **94.4̄ (= 17/18·100)**, %indels **0**, score **27 (= 17·2 − 7)** | Benson approximate def. |
| A2b | (same A2 sequence) perfect detector `FindMicrosatellites(...,minRepeats=3)` | reports only `CAG`×3 at pos 0 (fragmented) | contrast: perfect detector breaks the interrupted tract |
| A3 | `CACACATACACA`, period 2 (one substitution at idx 6) | consensus `CA`, copies 6.0, %matches **91.6̄ (= 11/12·100)**, %indels **0**, score **15 (= 11·2 − 7)** | Benson approximate def. |
| A4 | `CAGCAGCAGCAGCAGAGCAGCAGCAGCAG`, period 3 (one deletion) | consensus `CAG`, copies **29/3 = 9.6̄**, %matches **96.6̄ (= 29/30·100)**, %indels **3.3̄ (= 1/30·100)**, score **51 (= 29·2 − 7)** | Benson percent-indels statistic |
| A5 | `CACACACACA` at `minScore = 50` default | empty (score 20 < 50) | Benson "≥ 50 … reported" |
| A6 | A4 sequence at default `minScore` | reported (score 51 ≥ 50) | Benson "≥ 50 … reported" |
| A7 | `""` | empty | edge |
| A8 | `ACGTGCAT` (no repeat) | empty | edge |
| A9 | `minPeriod = 0` / `maxPeriod < minPeriod` | `ArgumentOutOfRangeException` | parameter validation |
| A10 | `null` DnaSequence | `ArgumentNullException` | parameter validation |
| A11 | determinism — same input twice → identical results | equal | determinism |

### 9.4 Coverage status

All A1–A11 implemented in `tests/Seqeron/Seqeron.Genomics.Tests/RepeatFinder_ApproximateTandemRepeats_Tests.cs`, exact assertions with `.Within(1e-9)` on the percentages/copy number. ✅ Done = 11 / 11; Remaining = 0.

### 9.5 Residual (honest)

The per-repeat **Bernoulli statistical-significance** measures are now reproduced (see §10). What remains unreproduced is TRF's probabilistic k-tuple distance-list **seeding** (the `R(d,k,pM)` 95% sum-of-heads percentile cut-off and the `W(d,pI)` random-walk band, whose values come from TRF's non-redistributable simulation tables) — a whole-genome-scale **performance** index, not a per-repeat-correctness gap; the deterministic exhaustive (start, period) scan already finds every candidate a seed would. See Evidence ASSUMPTION 1, 2 and 3.

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
| B4 | `ACACTGTG`, period 4 | 1 pair, PM **0**, meets-0.80 **false** | low-identity tract below PM = 0.80 |
| B5 | `CAGCAGCAGTAGCAGCAG`, period 3, PM-threshold 0.80 / 0.90 | meets 0.80 **true**, meets 0.90 **false** | custom PM threshold |
| B6 | `CAGCAGCAGCAGCAGAGCAGCAGCAGCAG` (deletion), period 3 | PI **> 0**; match+mismatch+indel = trials; PM+PI ≤ 1 | PI = average % indels |
| B7 | `CAGCAGCAGTAGCAGCAG`, period 3 | PM + mismatch-fraction + PI = **1.0** | Bernoulli outcomes partition the trials |
| B8 | `CAG`, period 3 | `ArgumentException` | model of two copies undefined for one copy |
| B9 | `null` / period 0 / PM 1.5 | `ArgumentNullException` / `ArgumentOutOfRangeException` ×2 | parameter validation |
| B10 | exposed defaults | `TrfDefaultMatchProbability` **0.80**, `TrfDefaultIndelProbability` **0.10** | Benson (1999) defaults |

### 10.3 Coverage status

All B1–B10 implemented in `tests/Seqeron/Seqeron.Genomics.Tests/RepeatFinder_ApproximateTandemRepeats_Tests.cs` (region "ComputeBernoulliStatistics"), exact assertions with `.Within(1e-9)` on the probabilities/expected matches. ✅ Done = 10 / 10; Remaining = 0. Invariants INV-09, INV-10.
