# kmer_entropy

Calculate k-mer entropy of a sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `kmer_entropy` |
| **Method ID** | `KmerAnalyzer.CalculateKmerEntropy` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates Shannon entropy based on k-mer frequencies in a sequence. Higher entropy values indicate more complexity and diversity in the k-mer distribution. This is useful for assessing sequence complexity at different scales by varying k.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L243](../../../../src/Seqeron/Seqeron.Genomics/KmerAnalyzer.cs#L243)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The sequence to analyze (min length: 1) |
| `k` | integer | No | K-mer length (default: 2, minimum: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `entropy` | number | K-mer entropy (bits) |
| `k` | integer | K-mer length used |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | K must be at least 1 |

## Examples

### Example 1: Dinucleotide entropy

**User Prompt:**
> What's the 2-mer entropy of "ATGCATGCAT"?

**Expected Tool Call:**
```json
{
  "tool": "kmer_entropy",
  "arguments": {
    "sequence": "ATGCATGCAT",
    "k": 2
  }
}
```

**Response:**
```json
{
  "entropy": 3.17,
  "k": 2
}
```

### Example 2: Low complexity sequence

**User Prompt:**
> Calculate k-mer entropy for "AAAAAAAAAA"

**Expected Tool Call:**
```json
{
  "tool": "kmer_entropy",
  "arguments": {
    "sequence": "AAAAAAAAAA",
    "k": 2
  }
}
```

**Response:**
```json
{
  "entropy": 0,
  "k": 2
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(4^k) for k-mer storage

## See Also

- [shannon_entropy](shannon_entropy.md) - Single-character entropy
- [linguistic_complexity](linguistic_complexity.md) - K-mer diversity measure
