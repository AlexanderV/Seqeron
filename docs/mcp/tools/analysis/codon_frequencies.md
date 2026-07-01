# codon_frequencies

Codon usage frequencies in a reading frame.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `codon_frequencies` |
| **Method ID** | `SequenceStatistics.CalculateCodonFrequencies` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes codon usage frequencies over the non-overlapping in-frame triplets of a DNA
sequence, following the Kazusa CUTG definition `frequency = count(codon) / total
counted codons` (Nakamura et al. 2000). Reading starts at `readingFrame` (0, 1, or 2);
trailing bases that do not complete a triplet are ignored, and codons containing any
non-ACGT base are excluded from the count. Only observed codons appear as keys.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L688](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L688)
- Evidence: `docs/Evidence/SEQ-CODON-FREQ-001-Evidence.md`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `dnaSequence` | string | Yes | DNA sequence (min length 1) |
| `readingFrame` | integer | No | Reading frame 0, 1, or 2 (default 0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `frequencies` | object | Map of codon → frequency (count/total), summing to 1 over observed codons |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | Reading frame must be 0, 1, or 2 |

## Examples

### Example 1: Frame 0

**User Prompt:**
> Codon frequencies of "ATGATGAAA" in frame 0.

**Expected Tool Call:**
```json
{
  "tool": "codon_frequencies",
  "arguments": { "dnaSequence": "ATGATGAAA", "readingFrame": 0 }
}
```

**Response:**
```json
{ "frequencies": { "ATG": 0.6666666666666666, "AAA": 0.3333333333333333 } }
```
Non-overlapping triplets ATG, ATG, AAA ⇒ ATG=2/3, AAA=1/3.

### Example 2: Frame 1

**Input:** `{ "dnaSequence": "ATGATGAAA", "readingFrame": 1 }`
→ **Response:** `{ "frequencies": { "TGA": 1.0 } }`
From index 1 the triplets are TGA, TGA; trailing `AA` is ignored.

## Performance

- **Time Complexity:** O(n). **Space Complexity:** O(distinct codons).

## See Also

- [dinucleotide_frequencies](dinucleotide_frequencies.md)
