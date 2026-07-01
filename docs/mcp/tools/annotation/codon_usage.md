# codon_usage

Count codon occurrences (in-frame, 5'→3') for a coding sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `codon_usage` |
| **Method ID** | `GenomeAnnotator.GetCodonUsage` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Counts codon occurrences in a coding DNA sequence. The sequence is read in frame from position 0 in
non-overlapping 3-nt windows; only codons over the alphabet A/C/G/T are counted, input is
uppercased, and any partial trailing codon is ignored. Returns a map from each observed codon to its
occurrence count.

## Core Documentation Reference

- Source: [GenomeAnnotator.cs#L1101](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs#L1101)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `dnaSequence` | string | Yes | Coding DNA sequence (length should be a multiple of 3) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `usage` | object | Codon → occurrence count |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Repeated start codon

`ATGATGAAATAA` → codons ATG, ATG, AAA, TAA.

**Response:**
```json
{ "usage": { "ATG": 2, "AAA": 1, "TAA": 1 } }
```

### Example 2: Partial trailing codon ignored

`ATGAT` → only the full codon `ATG` is counted; the trailing `AT` is dropped.

**Response:**
```json
{ "usage": { "ATG": 1 } }
```

## Performance

- **Time Complexity:** O(n)
- **Space Complexity:** O(distinct codons)

## See Also

- [find_orfs](find_orfs.md) - ORF detection
- [coding_potential](coding_potential.md) - CPAT coding-potential score
