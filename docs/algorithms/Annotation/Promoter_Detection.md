# Promoter Detection (Bacterial Promoter Motif Identification)

## Overview

Bacterial promoter detection identifies conserved sequence motifs in DNA that signal RNA polymerase binding sites for transcription initiation. In bacteria, promoters contain two primary consensus elements recognized by the sigma factor subunit of RNA polymerase holoenzyme.

## Biological Background

### Bacterial Promoter Architecture

In bacteria, the promoter contains two short sequence elements upstream of the transcription start site:

| Element | Position | Consensus | Function |
|---------|----------|-----------|----------|
| **-35 box** | ~35 bp upstream | **TTGACA** | Initial recognition by σ70 |
| **-10 box** (Pribnow box) | ~10 bp upstream | **TATAAT** | DNA melting, strand separation |

**Key characteristics** (per Wikipedia):
- The optimal spacing between -35 and -10 elements is **17 bp**
- On average, only **3 to 4 of the 6 base pairs** in each consensus are found in natural promoters
- Artificial promoters with complete conservation transcribe at **lower frequencies** than those with mismatches
- AT-rich -10 box facilitates DNA strand separation due to weaker hydrogen bonding

### Nucleotide Occurrence Probability (E. coli)

**-10 box (Pribnow box):**
| Position | T | A | T | A | A | T |
|----------|---|---|---|---|---|---|
| Probability | 77% | 76% | 60% | 61% | 56% | 82% |

**-35 box:**
| Position | T | T | G | A | C | A |
|----------|---|---|---|---|---|---|
| Probability | 69% | 79% | 61% | 56% | 54% | 54% |

## Algorithm Description

### Input
- DNA sequence (string)

### Output
- Collection of tuples: `(position, type, sequence, score)`
  - `position`: 0-based index in the sequence
  - `type`: "-35 box" or "-10 box"
  - `sequence`: the matched motif sequence
  - `score`: confidence score (0.0 to 1.0)

### Detection Strategy

The implementation uses a pattern-matching approach with consensus variants:

1. **-35 box variants**: TTGACA (full), TTGAC, TGACA, TTGA (partial)
2. **-10 box variants**: TATAAT (full), TATAA, TAAAT, TATA (partial)

**Scoring**: `score = motif_length / 6.0` (normalized to canonical length)

### Complexity
- **Time**: O(n × m) where n = sequence length, m = number of motif variants
- **Space**: O(1) excluding output

## Implementation Notes

### Current Implementation (`GenomeAnnotator.FindPromoterMotifs`)

```csharp
// Searches for both full consensus and partial matches
// Returns hits for -35 box and -10 box independently
// Does NOT verify spacing between -35 and -10 elements
```

**Limitations**:
1. Searches for each motif independently (no paired -35/-10 validation)
2. Does not verify the 17 bp optimal spacing constraint
3. Partial motifs (4-5 bp) may produce false positives
4. Case-insensitive matching (sequence is uppercased)

### Deviations from Literature

| Aspect | Literature | Implementation |
|--------|------------|----------------|
| Spacing validation | Optimal 17 bp between -35 and -10 | **Not enforced** |
| Mismatch tolerance | 2-3 mismatches typical | Only exact substring matches to variants |
| Scoring | Based on similarity to consensus | Based on motif length only |

## Test Considerations

### Edge Cases
1. **Empty sequence**: Should return empty collection
2. **No motifs present**: Should return empty collection
3. **Overlapping motifs**: Both should be reported
4. **Case variations**: Should handle mixed case
5. **Adjacent motifs**: Multiple hits at consecutive positions
6. **Partial matches**: Short variants (4 bp) should have lower scores

### Expected Behavior
- Full consensus matches (6 bp) → score = 1.0
- Partial matches (4-5 bp) → score < 1.0 (proportional to length)
- Multiple independent hits can be returned for both -35 and -10 boxes

## References

1. Wikipedia: [Promoter (genetics)](https://en.wikipedia.org/wiki/Promoter_(genetics)) - Bacterial promoter structure
2. Wikipedia: [Pribnow box](https://en.wikipedia.org/wiki/Pribnow_box) - -10 element consensus
3. Wikipedia: [TATA box](https://en.wikipedia.org/wiki/TATA_box) - Eukaryotic analog
4. Pribnow, D. (1975). "Nucleotide sequence of an RNA polymerase binding site at an early T7 promoter." PNAS 72(3): 784-788.
5. Harley, C.B. & Reynolds, R.P. (1987). "Analysis of E. coli promoter sequences." Nucleic Acids Research 15(5): 2343-2361.

---

**Test Unit ID**: ANNOT-PROM-001
**Last Updated**: 2026-01-23
