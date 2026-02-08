# Test Specification: REP-STR-001

**Test Unit:** Microsatellite Detection (STR)
**Canonical Class:** `RepeatFinder`
**Primary Method:** `FindMicrosatellites`
**Status:** Complete
**Last Updated:** 2026-01-22

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
| M11 | Invalid parameters throw | Defensive programming | Error handling |
| M12 | Redundant unit filtering | Implementation: "ATAT" → "AT"×2 | Avoid duplicate reporting |

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

### Current State (RepeatFinderTests.cs)

| Test Method | Coverage | Status | Action |
|-------------|----------|--------|--------|
| `FindMicrosatellites_MononucleotideRepeat_FindsRepeat` | M01 | ✓ Good | Keep |
| `FindMicrosatellites_DinucleotideRepeat_FindsRepeat` | M02 | ✓ Good | Keep |
| `FindMicrosatellites_TrinucleotideRepeat_CAGExpansion` | M03 | ✓ Good | Keep |
| `FindMicrosatellites_TetranucleotideRepeat_FindsRepeat` | M04 | ✓ Good | Keep |
| `FindMicrosatellites_NoRepeats_ReturnsEmpty` | Partial | Weak | Rename/fix |
| `FindMicrosatellites_MultipleDifferentRepeats_FindsAll` | S01 | ✓ Good | Keep |
| `FindMicrosatellites_MinRepeatsFilter_RespectsThreshold` | M06 | ✓ Good | Keep |
| `FindMicrosatellites_StringOverload_Works` | S02 | ✓ Good | Keep |
| `FindMicrosatellites_EmptySequence_ReturnsEmpty` | M05 | ✓ Good | Keep |
| `FindMicrosatellites_FullSequenceProperty_ReturnsCorrectSequence` | M08 | ✓ Good | Keep |
| `FindMicrosatellites_NullSequence_ThrowsException` | M10 | ✓ Good | Keep |
| `FindMicrosatellites_InvalidMinUnitLength_ThrowsException` | M11 | ✓ Good | Keep |
| `FindMicrosatellites_InvalidMaxUnitLength_ThrowsException` | M11 | ✓ Good | Keep |
| `FindMicrosatellites_InvalidMinRepeats_ThrowsException` | M11 | ✓ Good | Keep |

### Missing Tests (All Closed)
| Test | Status |
|------|--------|
| Hexanucleotide detection | ✅ Covered |
| RepeatType enum validation | ✅ Covered |
| Position range invariant | ✅ Covered |
| Invariant assertions (Assert.Multiple) | ✅ Covered |
| Case insensitivity | ✅ Covered |
| Cancellation smoke test | ✅ Covered (PerformanceExtensionsTests) |

### Consolidation Plan
1. ~~Create dedicated `RepeatFinder_Microsatellite_Tests.cs`~~ ✅ Done
2. Keep cancellation smoke test in PerformanceExtensionsTests
3. ~~Add missing Must tests~~ ✅ Done
4. ~~Strengthen assertions with invariant checks~~ ✅ Done

---

## 7. Assumptions

| # | Assumption | Justification |
|---|------------|---------------|
| A1 | Unit lengths 1-6 define microsatellites | Wikipedia states "1-6 or up to 10" - implementation uses 1-6 |
| A2 | Minimum repeats default 3 is sensible | Wikipedia: "typically repeated 5-50 times" - 3 is conservative |
| A3 | Redundant unit filtering is correct behavior | Avoids duplicate reporting of same repeat |

---

## 8. Open Questions

None - behavior is well-documented.

---

## 9. Validation Checklist

- [x] Evidence sources documented
- [x] Must tests defined with evidence
- [x] Existing tests audited
- [x] Consolidation plan created
- [x] Invariants specified
- [ ] Tests implemented
- [ ] Tests passing
- [ ] Zero warnings
