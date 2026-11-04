# complexity_shannon

Calculate DNA Shannon entropy.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `complexity_shannon` |
| **Method ID** | `SequenceComplexity.CalculateShannonEntropy` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates Shannon entropy H = −Σ pᵢ·log₂ pᵢ (bits per base) of a DNA or RNA sequence over the alphabet A/C/G/T/U. Counting is case-insensitive and RNA U is counted as T (the same nucleotide class), so an RNA and its DNA transcript give the same value. Any other character (N, IUPAC ambiguity codes, gaps) is ignored, i.e. excluded from the counts and from the total. Maximum entropy is 2 bits (log2(4)), indicating equal distribution of all four nucleotides. Lower entropy values indicate biased composition or repetitive sequences.

Cross-check (scipy `scipy.stats.entropy(counts, base=2)`): `"ACGU"` → 2.0, `"AUAUAU"` → 1.0, `"AACG-"` → 1.5.

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L214](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs#L214)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA or RNA sequence to analyze (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `entropy` | number | Shannon entropy (0 to 2 bits for DNA) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Balanced sequence (max entropy)

**User Prompt:**
> Calculate Shannon entropy of "ATGCATGCATGC"

**Expected Tool Call:**
```json
{
  "tool": "complexity_shannon",
  "arguments": {
    "sequence": "ATGCATGCATGC"
  }
}
```

**Response:**
```json
{
  "entropy": 2.0
}
```

### Example 2: Low complexity sequence

**User Prompt:**
> What's the entropy of "AAAAAAAAAA"?

**Expected Tool Call:**
```json
{
  "tool": "complexity_shannon",
  "arguments": {
    "sequence": "AAAAAAAAAA"
  }
}
```

**Response:**
```json
{
  "entropy": 0
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1) for nucleotide counts

## See Also

- [shannon_entropy](shannon_entropy.md) - SequenceStatistics version
- [complexity_linguistic](complexity_linguistic.md) - Linguistic complexity
- [complexity_kmer_entropy](complexity_kmer_entropy.md) - K-mer entropy
