# kmer_distance

Alignment-free Euclidean distance between two sequences' k-mer composition.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `kmer_distance` |
| **Method ID** | `KmerAnalyzer.KmerDistance` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the **alignment-free k-mer distance**: the Euclidean distance between the
two sequences' normalized k-mer frequency vectors, taken over the union of k-mers
occurring in either sequence (a k-mer absent from one sequence contributes a 0
component). Identical sequences yield 0; sharing no k-mers yields the largest
distance. This is the frequency (relative-count) word-composition distance of
Zielezinski et al. (2017) / Vinga & Almeida (2003). Counting is case-insensitive.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L212](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L212)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `seq1` | string | Yes | First sequence (min length 1) |
| `seq2` | string | Yes | Second sequence (min length 1) |
| `k` | integer | Yes | k-mer length (> 0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `distance` | number | Non-negative Euclidean distance between the two frequency vectors |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | k must be positive |

## Examples

### Example 1: Identical sequences (distance = 0)

**User Prompt:**
> k-mer distance between "ACGTACGT" and itself for k=2.

**Expected Tool Call:**
```json
{
  "tool": "kmer_distance",
  "arguments": { "seq1": "ACGTACGT", "seq2": "ACGTACGT", "k": 2 }
}
```

**Response:**
```json
{ "distance": 0.0 }
```

### Example 2: Disjoint compositions (distance = √2)

**User Prompt:**
> k-mer distance between "AAAA" and "TTTT" for k=1.

**Expected Tool Call:**
```json
{
  "tool": "kmer_distance",
  "arguments": { "seq1": "AAAA", "seq2": "TTTT", "k": 1 }
}
```

**Response:**
```json
{ "distance": 1.4142135623730951 }
```
Frequency vectors {A:1} and {T:1} are orthogonal ⇒ √(1² + 1²) = √2.

## Performance

- **Time Complexity:** O(n) to build the two k-mer tables.
- **Space Complexity:** O(distinct k-mers).

## See Also

- [kmer_frequencies](kmer_frequencies.md) — the underlying frequency vectors
- [count_kmers](count_kmers.md) — raw counts
