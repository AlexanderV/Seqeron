# complexity_kmer_entropy

Calculate k-mer based entropy for DNA sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `complexity_kmer_entropy` |
| **Method ID** | `SequenceComplexity.CalculateKmerEntropy` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates Shannon entropy using k-mer frequencies in a DNA sequence. This provides a more nuanced complexity measure than single-base entropy by considering local sequence context. Higher entropy values indicate more diverse k-mer composition.

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L128](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/SequenceComplexity.cs#L128)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA sequence to analyze (must be valid DNA) |
| `k` | integer | No | K-mer size (default: 2 for dinucleotides, minimum: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `entropy` | number | K-mer entropy (bits) |
| `k` | integer | K-mer size used |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | K must be at least 1 |
| 2001 | Invalid DNA sequence |

## Examples

### Example 1: Dinucleotide entropy

**User Prompt:**
> Calculate 2-mer entropy of "ATGCATGCAT"

**Expected Tool Call:**
```json
{
  "tool": "complexity_kmer_entropy",
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
> What's the k-mer entropy of "AAAAAAAAAA"?

**Expected Tool Call:**
```json
{
  "tool": "complexity_kmer_entropy",
  "arguments": {
    "sequence": "AAAAAAAAAA"
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

- [kmer_entropy](kmer_entropy.md) - KmerAnalyzer version
- [complexity_shannon](complexity_shannon.md) - Single-base entropy
- [complexity_linguistic](complexity_linguistic.md) - Linguistic complexity
