# dinucleotide_frequencies

Normalized adjacent-dinucleotide frequencies of a nucleotide sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `dinucleotide_frequencies` |
| **Method ID** | `SequenceStatistics.CalculateDinucleotideFrequencies` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the normalized frequency of each adjacent dinucleotide, `f_XY = count(XY) / (N − 1)`,
over the alphabet {A, T, G, C, U}. Dinucleotides containing a non-alphabet character are
excluded. Frequencies sum to 1 (Karlin genomic-signature convention). Counting is
case-insensitive; sequences shorter than 2 bp yield an empty map.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L602](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L602)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence (min length 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `frequencies` | object | Map of dinucleotide → frequency (sums to 1) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Homopolymer

**Input:** `{ "sequence": "AAAA" }`
→ AA,AA,AA over 3 positions → **Response:** `{ "frequencies": { "AA": 1.0 } }`

### Example 2: Alternating

**Input:** `{ "sequence": "ATAT" }`
→ AT,TA,AT → **Response:** `{ "frequencies": { "AT": 0.6667, "TA": 0.3333 } }`

## Performance

- **Time Complexity:** O(n). **Space Complexity:** O(distinct dinucleotides).

## See Also

- [dinucleotide_ratios](dinucleotide_ratios.md)
- [codon_frequencies](codon_frequencies.md)
