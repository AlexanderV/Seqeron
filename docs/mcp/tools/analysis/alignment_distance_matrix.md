# alignment_distance_matrix

Pairwise distance matrix of a multiple alignment as EMBOSS 6.6.0 distmat (every correction and option, verified on 1 100 + 600 runs of the installed program, single-precision arithmetic as in distmat.c): matches = columns where both rows carry the same non-gap symbol, g = columns where either row has a gap; 'uncorrected' D = 1 - matches/(L - g + g*gapWeight); 'jukes-cantor' -b*ln(1 - D/b) with b = 3/4 (nucleotide) or 19/20 (protein); 'kimura' = Kimura 2-parameter (nucleotides) / Kimura protein -ln(1 - D - 0.2 D^2) (proteins); 'tamura', 'tajima-nei' and 'jin-nei' (Gamma distance; calculateA / parameterA) are nucleotide methods. ambiguous counts ambiguity codes (proteins: B/Z/J/X weights; nucleotides: distmat's ajBaseAlphacharCompare with its own defects, and a gap in the first sequence of a pair makes the program loop forever = tool error). position = 123 (all codon bases, default), 12, or 1 / 2 / 3 (nucleotides); tajima-nei and protein kimura with a position other than 123 read beyond the program's buffers (undefined) = tool error. Fractions, not percent; entries with a zero denominator or a saturated logarithm are null. Rows are read as EMBOSS reads an alignment (A4-3): letters and '*' are kept, '.' '~' '-' are gaps, '?' '#' '+' are the unknown residue, every other character (digits, blanks, '_', punctuation) is dropped.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `alignment_distance_matrix` |
| **Method ID** | `MultipleAlignmentStatistics.CalculateDistanceMatrix` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Pairwise distance matrix of a multiple alignment as EMBOSS 6.6.0 distmat (every correction and option, verified on 1 100 + 600 runs of the installed program, single-precision arithmetic as in distmat.c): matches = columns where both rows carry the same non-gap symbol, g = columns where either row has a gap; 'uncorrected' D = 1 - matches/(L - g + g*gapWeight); 'jukes-cantor' -b*ln(1 - D/b) with b = 3/4 (nucleotide) or 19/20 (protein); 'kimura' = Kimura 2-parameter (nucleotides) / Kimura protein -ln(1 - D - 0.2 D^2) (proteins); 'tamura', 'tajima-nei' and 'jin-nei' (Gamma distance; calculateA / parameterA) are nucleotide methods. ambiguous counts ambiguity codes (proteins: B/Z/J/X weights; nucleotides: distmat's ajBaseAlphacharCompare with its own defects, and a gap in the first sequence of a pair makes the program loop forever = tool error). position = 123 (all codon bases, default), 12, or 1 / 2 / 3 (nucleotides); tajima-nei and protein kimura with a position other than 123 read beyond the program's buffers (undefined) = tool error. Fractions, not percent; entries with a zero denominator or a saturated logarithm are null. Rows are read as EMBOSS reads an alignment (A4-3): letters and '*' are kept, '.' '~' '-' are gaps, '?' '#' '+' are the unknown residue, every other character (digits, blanks, '_', punctuation) is dropped.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MultipleAlignmentStatistics.cs#L220](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MultipleAlignmentStatistics.cs#L220)
- Algorithm: [Alignment_Statistics.md §5.6](../../../algorithms/Alignment/Alignment_Statistics.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | At least two aligned sequences of equal length. |
| `method` | string | No | 'uncorrected' (default), 'jukes-cantor', 'kimura', 'tamura', 'tajima-nei' or 'jin-nei'. |
| `residueType` | string | No | 'auto' (default; EMBOSS decides from the first row), 'nucleotide' or 'protein'. |
| `gapWeight` | number | No | distmat -gapweight (default 0; uncorrected and jukes-cantor). |
| `ambiguous` | boolean | No | distmat -ambiguous (default false). |
| `position` | integer | No | distmat -position: 123 (default), 12, or for nucleotides 1, 2, 3. |
| `calculateA` | boolean | No | distmat -calculatea: calculate the Jin-Nei parameter a instead of using parameterA (default false). |
| `parameterA` | number | No | distmat -parametera: the Jin-Nei Gamma parameter a (default 1.0). |
| `sequenceQualifiers` | string | No | EMBOSS sequence-input qualifiers as on the command line: -sbegin 3 -send -2 -sreverse -slower -supper -scircular -snucleotide -sprotein (begin / end 1-based, negative = from the end; -sreverse reverse-complements a nucleotide sequence, the range then refers to the original). Measured on EMBOSS 6.6.0 per program: distmat ignores -sbegin / -send; -sreverse reverse-complements every row of a nucleotide set; the set is nucleotide when its first row is (or -snucleotide); a protein set uses proteinMethod. |
| `proteinMethod` | string | No | distmat -protmethod for a protein set when sequenceQualifiers is given: same names as method (default: the same as method). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `method` | string | Method used. |
| `matrix` | array<array<number>> | Symmetric distance matrix (fractions, zero diagonal; null = undefined or saturated). |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least two sequences are required |
| 1002 | method must be 'uncorrected', 'jukes-cantor', 'kimura', 'tamura', 'tajima-nei' or 'jin-nei' |

## Examples

### Example 1: EMBOSS infoalign / distmat reference values

**Tool Call:**
```json
{
  "tool": "alignment_distance_matrix",
  "arguments": {
    "sequences": [
      "ACGT-ACGTAC",
      "ACGTTACG-AC",
      "AGGT-ACGTAC",
      "A-GTTTCGTAA"
    ]
  }
}
```

**Result (checked against EMBOSS 6.6.0):** matrix[0][2] = 0.1, matrix[1][2] = 0.1111 (EMBOSS distmat 6.6.0: 10.00, 11.11 percent).
