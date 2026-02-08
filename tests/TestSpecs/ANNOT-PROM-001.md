# Test Specification: ANNOT-PROM-001

**Test Unit**: Promoter Detection
**Algorithm Group**: Annotation
**Canonical Method**: `GenomeAnnotator.FindPromoterMotifs(string dnaSequence)`
**Date**: 2026-01-23
**Status**: Draft

---

## Summary

Tests for bacterial promoter motif detection, specifically identifying -35 box (TTGACA) and -10 box (TATAAT) consensus sequences and their variants.

---

## Evidence Sources

| Source | URL/Reference | Key Information |
|--------|---------------|-----------------|
| Wikipedia: Promoter (genetics) | https://en.wikipedia.org/wiki/Promoter_(genetics) | -35 box: TTGACA, -10 box: TATAAT, 17 bp optimal spacing |
| Wikipedia: Pribnow box | https://en.wikipedia.org/wiki/Pribnow_box | -10 element consensus TATAAT, nucleotide probabilities |
| Wikipedia: TATA box | https://en.wikipedia.org/wiki/TATA_box | Eukaryotic analog, Pribnow box is bacterial homolog |
| Pribnow (1975) | PNAS 72(3):784-788 | Original -10 box identification |
| Harley & Reynolds (1987) | NAR 15(5):2343-2361 | E. coli promoter analysis |

---

## Canonical Tests

### Method: `FindPromoterMotifs(string dnaSequence)`

**Returns**: `IEnumerable<(int position, string type, string sequence, double score)>`

---

## Must Tests (Required - Evidence-Based)

| ID | Test Name | Description | Evidence |
|----|-----------|-------------|----------|
| M01 | FindPromoterMotifs_FullMinus35Consensus_ReturnsCorrectHit | TTGACA returns position, type="-35 box", score=1.0 | Wikipedia: consensus is TTGACA |
| M02 | FindPromoterMotifs_FullMinus10Consensus_ReturnsCorrectHit | TATAAT returns position, type="-10 box", score=1.0 | Wikipedia: Pribnow box is TATAAT |
| M03 | FindPromoterMotifs_PartialMinus35_ReturnsLowerScore | TTGAC (5 bp) returns score < 1.0 | Implementation: score = length/6 |
| M04 | FindPromoterMotifs_PartialMinus10_ReturnsLowerScore | TATAA (5 bp) returns score < 1.0 | Implementation: score = length/6 |
| M05 | FindPromoterMotifs_NoMotifs_ReturnsEmpty | Sequence without motifs (e.g., CCCCCC) returns empty | Logic invariant |
| M06 | FindPromoterMotifs_EmptySequence_ReturnsEmpty | Empty string returns empty collection | Edge case |
| M07 | FindPromoterMotifs_MixedCase_HandlesCorrectly | Lowercase input should find motifs | Implementation: uses ToUpperInvariant |
| M08 | FindPromoterMotifs_MultipleMotifsOfSameType_ReturnsAll | Multiple -35 boxes in sequence returns all positions | Algorithm behavior |
| M09 | FindPromoterMotifs_BothMotifTypes_ReturnsBothTypes | Sequence with both -35 and -10 returns hits for both | Core functionality |
| M10 | FindPromoterMotifs_CorrectPositionReporting | Position is 0-based index of motif start | Algorithm invariant |

---

## Should Tests (Recommended)

| ID | Test Name | Description | Rationale |
|----|-----------|-------------|-----------|
| S01 | FindPromoterMotifs_AdjacentMotifs_ReportsAllPositions | Overlapping partial motifs reported | Completeness |
| S02 | FindPromoterMotifs_Score_EqualsLengthDividedBySix | Verify score calculation formula | Implementation verification |
| S03 | FindPromoterMotifs_AllMinus35Variants_Detected | TTGACA, TTGAC, TGACA, TTGA all detected | Full variant coverage |
| S04 | FindPromoterMotifs_AllMinus10Variants_Detected | TATAAT, TATAA, TAAAT, TATA all detected | Full variant coverage |

---

## Could Tests (Optional)

| ID | Test Name | Description | Rationale |
|----|-----------|-------------|-----------|
| C01 | FindPromoterMotifs_RealPromoterSequence_FindsExpectedMotifs | Test with known E. coli promoter | Biological validation |
| C02 | FindPromoterMotifs_LongSequence_PerformsEfficiently | Performance with genome-scale input | Performance baseline |

---

## Audit of Existing Tests

### Location: `GenomeAnnotatorTests.cs` (lines 156-186)

| Existing Test | Status | Assessment | Action |
|---------------|--------|------------|--------|
| `FindPromoterMotifs_FindsMinus35Box` | Weak | Only checks existence, not position/score | Strengthen |
| `FindPromoterMotifs_FindsMinus10Box` | Weak | Only checks existence, not position/score | Strengthen |
| `FindPromoterMotifs_NoMotifs_ReturnsEmpty` | Covered | Adequate for M05 | Keep, refactor |

### Consolidation Plan

1. **Create new canonical file**: `GenomeAnnotator_PromoterMotif_Tests.cs`
2. **Move and strengthen**: Extract promoter tests from `GenomeAnnotatorTests.cs`
3. **Remove from original**: Delete `#region FindPromoterMotifs Tests` from `GenomeAnnotatorTests.cs`
4. ~~**Add missing tests**: M01-M10, S01-S04~~ ✅ All added

---

## Open Questions / Decisions

| # | Question | Decision |
|---|----------|----------|
| 1 | Should tests validate -35/-10 spacing? | NO - Implementation does not enforce spacing; document as limitation |
| 2 | How to handle partial motif overlap? | All matches reported independently per implementation |

---

## Assumptions

| # | Assumption | Justification |
|---|------------|---------------|
| A01 | Score formula: length/6 is implementation-specific | Verified from source code; not from literature |

---

## Test File Structure

```
GenomeAnnotator_PromoterMotif_Tests.cs
├── Consensus Motif Detection
│   ├── M01: Full -35 box
│   ├── M02: Full -10 box
│   └── M09: Both types
├── Partial Motif Detection
│   ├── M03: Partial -35 box
│   ├── M04: Partial -10 box
│   └── S03/S04: All variants
├── Edge Cases
│   ├── M05: No motifs
│   ├── M06: Empty sequence
│   └── M07: Mixed case
├── Position and Score Validation
│   ├── M08: Multiple same-type motifs
│   ├── M10: Position correctness
│   └── S02: Score formula
└── Adjacent/Overlapping
    └── S01: Adjacent motifs
```
