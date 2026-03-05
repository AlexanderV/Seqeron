# Test Specification: ANNOT-PROM-001

**Test Unit**: Promoter Detection
**Algorithm Group**: Annotation
**Canonical Method**: `GenomeAnnotator.FindPromoterMotifs(string dnaSequence)`
**Date**: 2026-03-05
**Status**: Final

---

## Summary

Tests for bacterial promoter motif detection, specifically identifying -35 box (TTGACA) and -10 box (TATAAT) consensus sequences and their variants. Scoring uses position-specific nucleotide occurrence probabilities from E. coli promoter analysis.

---

## Evidence Sources

| Source | URL/Reference | Key Information |
|--------|---------------|-----------------|
| Wikipedia: Promoter (genetics) | https://en.wikipedia.org/wiki/Promoter_(genetics) | -35 box: TTGACA, -10 box: TATAAT, 17 bp optimal spacing, nucleotide probabilities |
| Wikipedia: Pribnow box | https://en.wikipedia.org/wiki/Pribnow_box | -10 element consensus TATAAT, alternative nucleotide probabilities |
| Wikipedia: TATA box | https://en.wikipedia.org/wiki/TATA_box | Eukaryotic analog, Pribnow box is bacterial homolog |
| Pribnow (1975) | PNAS 72(3):784-788 | Original -10 box identification |
| Harley & Reynolds (1987) | NAR 15(5):2343-2361 | E. coli promoter analysis, nucleotide occurrence probabilities |

---

## Scoring Mechanism

Scoring uses E. coli position-specific nucleotide occurrence probabilities
from Wikipedia "Promoter (genetics)" / Harley & Reynolds (1987).

**-35 box**: T(69%) T(79%) G(61%) A(56%) C(54%) A(54%) — total weight 3.73
**-10 box**: T(77%) A(76%) T(60%) A(61%) A(56%) T(82%) — total weight 4.12

**Formula**: `score = sum(matched position probabilities) / sum(all 6 consensus probabilities)`

Variants are consecutive consensus substrings (prefix-5bp, suffix-5bp, prefix-4bp):

| -35 box variant | Positions | Score |
|-----------------|-----------|-------|
| TTGACA | 1–6 (full) | 1.000 |
| TTGAC | 1–5 | 0.855 |
| TGACA | 2–6 | 0.815 |
| TTGA | 1–4 | 0.710 |

| -10 box variant | Positions | Score |
|-----------------|-----------|-------|
| TATAAT | 1–6 (full) | 1.000 |
| TATAA | 1–5 | 0.801 |
| ATAAT | 2–6 | 0.813 |
| TATA | 1–4 | 0.665 |

> **Note**: Wikipedia "Pribnow box" lists slightly different probabilities for the -10 element:
> T(82%) A(89%) T(52%) A(59%) A(49%) T(89%). Both come from E. coli studies;
> we use the Promoter (genetics) article values for consistency across both boxes.

---

## Canonical Tests

### Method: `FindPromoterMotifs(string dnaSequence)`

**Returns**: `IEnumerable<(int position, string type, string sequence, double score)>`

---

## Must Tests (Required — Evidence-Based)

| ID | Test Name | Description | Evidence |
|----|-----------|-------------|----------|
| M01 | FindPromoterMotifs_FullMinus35Consensus_ReturnsCorrectHit | TTGACA → position, type="-35 box", score=1.0 | Wikipedia: consensus is TTGACA |
| M02 | FindPromoterMotifs_FullMinus10Consensus_ReturnsCorrectHit | TATAAT → position, type="-10 box", score=1.0 | Wikipedia: Pribnow box is TATAAT |
| M03 | FindPromoterMotifs_PartialMinus35_ReturnsLowerScore | TTGAC (positions 1–5) → score ≈ 0.855 | Probability-weighted: (69+79+61+56+54)/373 |
| M04 | FindPromoterMotifs_PartialMinus10_ReturnsLowerScore | TATAA (positions 1–5) → score ≈ 0.801 | Probability-weighted: (77+76+60+61+56)/412 |
| M05 | FindPromoterMotifs_NoMotifs_ReturnsEmpty | Sequence without motifs (e.g., CCCCCC) → empty | Logic invariant |
| M06 | FindPromoterMotifs_EmptySequence_ReturnsEmpty | Empty string → empty collection | Edge case |
| M07 | FindPromoterMotifs_MixedCase_HandlesCorrectly | Lowercase input should find motifs | Uses ToUpperInvariant |
| M08 | FindPromoterMotifs_MultipleMotifsOfSameType_ReturnsAll | Multiple -35 boxes → all positions reported | Algorithm behavior |
| M09 | FindPromoterMotifs_BothMotifTypes_ReturnsBothTypes | Sequence with both -35 and -10 → hits for both | Core functionality |
| M10 | FindPromoterMotifs_CorrectPositionReporting | Position is 0-based index of motif start | Algorithm invariant |

---

## Should Tests (Recommended)

| ID | Test Name | Description | Rationale |
|----|-----------|-------------|-----------|
| S01 | FindPromoterMotifs_AdjacentMotifs_ReportsAllPositions | Overlapping partial motifs reported | Completeness |
| S02 | FindPromoterMotifs_Score_ReflectsPositionProbabilityWeights | Verify probability-weighted score formula | Literature-based verification |
| S03 | FindPromoterMotifs_AllMinus35Variants_Detected | TTGACA, TTGAC, TGACA, TTGA all detected with correct scores | Full variant coverage |
| S04 | FindPromoterMotifs_AllMinus10Variants_Detected | TATAAT, TATAA, ATAAT, TATA all detected with correct scores | Full variant coverage |

---

## Could Tests (Optional)

| ID | Test Name | Description | Rationale |
|----|-----------|-------------|-----------|
| C01 | FindPromoterMotifs_RealPromoterSequence_FindsExpectedMotifs | Test with known E. coli promoter | Biological validation |
| C02 | FindPromoterMotifs_LongSequence_PerformsEfficiently | Performance with genome-scale input | Performance baseline |

---

## Deviations from Literature

| Aspect | Literature | Implementation | Justification |
|--------|------------|----------------|---------------|
| Spacing validation | Optimal 17 bp between -35 and -10 | Not enforced | Independent motif search; spacing is a higher-level promoter prediction feature |
| Mismatch tolerance | 2–3 mismatches typical in real promoters | Only exact substring matches to consensus substrings | Exact matching of known conserved substrings avoids ambiguity; mismatch-tolerant detection would require a full PWM/HMM approach |

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
