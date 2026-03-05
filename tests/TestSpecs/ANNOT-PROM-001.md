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
| M01 | FullMinus35Consensus_ReturnsCorrectHit | TTGACA → exact 4 hits (full + 3 sub-variants), position, type="-35 box", score=1.0 | Wikipedia: consensus is TTGACA |
| M02 | FullMinus10Consensus_ReturnsCorrectHit | TATAAT → exact 4 hits (full + 3 sub-variants), position, type="-10 box", score=1.0 | Wikipedia: Pribnow box is TATAAT |
| M03 | PartialMinus35_ReturnsLowerScore | TTGAC (positions 1–5) → score = 0.855 | Probability-weighted: (69+79+61+56+54)/373 |
| M04 | PartialMinus10_ReturnsLowerScore | TATAA (positions 1–5) → score = 0.801 | Probability-weighted: (77+76+60+61+56)/412 |
| M05 | NoMotifs_ReturnsEmpty | Sequence without motifs (all C's) → empty | Logic invariant |
| M06 | EmptySequence_ReturnsEmpty | Empty string → empty collection | Edge case |
| M07 | MixedCase_HandlesCorrectly | Lowercase input → exact same position, type, score as uppercase | Uses ToUpperInvariant |
| M08 | MultipleMotifsOfSameType_ReturnsAll | Two TTGACA → exact count=2, exact positions | Algorithm behavior |
| M09 | BothMotifTypes_ReturnsBothTypes | TTGACA + TATAAT → exact 4 -35 hits + 4 -10 hits, both full scores=1.0 | Core functionality |

---

## Should Tests (Recommended)

| ID | Test Name | Description | Rationale |
|----|-----------|-------------|-----------|
| S01 | OverlappingMotifs_ReportsAllVariants | TATAAT → exact 4 hits with exact positions (0,0,1,0) | Variant decomposition |
| S02 | Score_ReflectsPositionProbabilityWeights | All 4 variants from TTGACA: exact scores, monotonic ordering 4bp<5bp<6bp | Literature-based verification |
| S03 | AllMinus35Variants_Detected | TTGACA, TTGAC, TGACA, TTGA — exact scores (1.000, 0.855, 0.815, 0.710) | Full variant coverage |
| S04 | AllMinus10Variants_Detected | TATAAT, TATAA, ATAAT, TATA — exact scores (1.000, 0.801, 0.813, 0.665) | Full variant coverage |

---

## Could Tests (Optional)

| ID | Test Name | Description | Rationale |
|----|-----------|-------------|-----------|
| C01 | RealisticPromoterSequence_FindsMotifs | Consensus -35/-10 with spacing → exact positions and scores | Biological plausibility |

---

## Deviations from Literature

| Aspect | Literature | Implementation | Justification |
|--------|------------|----------------|---------------|
| Spacing validation | Optimal 17 bp between -35 and -10 | Not enforced | Independent motif search; spacing is a higher-level promoter prediction feature |
| Mismatch tolerance | 2–3 mismatches typical in real promoters | Only exact substring matches to consensus substrings | Exact matching of known conserved substrings avoids ambiguity; mismatch-tolerant detection would require a full PWM/HMM approach |

---

## Test File Structure

```
GenomeAnnotator_PromoterMotif_Tests.cs (20 tests)
├── Consensus Motif Detection
│   ├── M01: Full -35 box (exact 4 hits, position, score)
│   ├── M02: Full -10 box (exact 4 hits, position, score)
│   └── M09: Both types (exact hit counts, full scores)
├── Partial Motif Detection
│   ├── M03: Partial -35 box (exact score 0.855)
│   ├── M04: Partial -10 box (exact score 0.801)
│   ├── S03: All -35 variants (4 TestCases with exact scores)
│   └── S04: All -10 variants (4 TestCases with exact scores)
├── Edge Cases
│   ├── M05: No motifs → empty
│   ├── M06: Empty sequence → empty
│   └── M07: Mixed case (exact position, type, score)
├── Score Validation
│   ├── M08: Multiple same-type motifs (exact count + positions)
│   └── S02: Score formula (all variants, monotonic ordering)
├── Overlapping
│   └── S01: Overlapping variants (exact count=4, exact positions)
└── Real-World Sequence
    └── C01: Consensus sequences with spacing (exact positions + scores)
```
