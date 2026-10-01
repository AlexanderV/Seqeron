# analyze_kmers

Aggregate k-mer composition statistics for a sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `analyze_kmers` |
| **Method ID** | `KmerAnalyzer.AnalyzeKmers` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes comprehensive k-mer composition statistics for a sequence: the total
number of overlapping k-mers (`L − k + 1`), the number of **distinct** k-mers, the
maximum/minimum/average multiplicity, and the Shannon entropy of the k-mer
frequency distribution (`E_k = −Σ p·log₂ p` with `p = mult/(L−k+1)`, in bits).
Counting is case-insensitive. When `k` exceeds the sequence length no k-mer exists
and every statistic is 0. Field definitions follow Jellyfish `stats`
(`sub_commands/stats_main.cc`): Total, Distinct, Unique (= count 1, here
`singletonKmers`), Max_count; `uniqueKmers` is a legacy name for the **distinct**
count. Optional `lowerCount`/`upperCount` apply Jellyfish's `-L`/`-U` filters; all
statistics are then computed over the retained k-mers.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L680](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L680)
- Evidence: `docs/Evidence/KMER-STATS-001-Evidence.md`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence to analyze (min length 1) |
| `k` | integer | Yes | k-mer length (> 0) |
| `lowerCount` | integer | No | Ignore k-mers with count < lowerCount (Jellyfish `-L`; default 0) |
| `upperCount` | integer | No | Ignore k-mers with count > upperCount (Jellyfish `-U`; default unbounded) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `totalKmers` | integer | Number of k-mers including multiplicity (Jellyfish Total), `L − k + 1` when unfiltered |
| `uniqueKmers` | integer | Number of **distinct** k-mers (legacy name; = `distinctKmers`) |
| `distinctKmers` | integer | Number of distinct k-mers (Jellyfish Distinct) |
| `singletonKmers` | integer | Number of k-mers occurring exactly once (Jellyfish Unique) |
| `maxCount` | integer | Maximum k-mer multiplicity |
| `minCount` | integer | Minimum k-mer multiplicity |
| `averageCount` | number | Exact mean multiplicity `total/distinct` |
| `entropy` | number | Shannon entropy `−Σ p·log₂ p` (bits) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | k must be positive |

## Examples

### Example 1: Monomer table (Wikipedia GTAGAGCTGT, k=1)

**User Prompt:**
> Give me k-mer statistics of "GTAGAGCTGT" for k=1.

**Expected Tool Call:**
```json
{
  "tool": "analyze_kmers",
  "arguments": { "sequence": "GTAGAGCTGT", "k": 1 }
}
```

**Response:**
```json
{
  "totalKmers": 10,
  "uniqueKmers": 4,
  "maxCount": 4,
  "minCount": 1,
  "averageCount": 2.5,
  "entropy": 1.846439344671,
  "distinctKmers": 4,
  "singletonKmers": 1
}
```
Monomer table: G=4, T=3, A=2, C=1. Entropy of {0.4,0.3,0.2,0.1} = 1.84643934… bits.

### Example 2: All-distinct trimers (entropy = log₂ 8)

**User Prompt:**
> k-mer statistics of "GTAGAGCTGT" for k=3.

**Expected Tool Call:**
```json
{
  "tool": "analyze_kmers",
  "arguments": { "sequence": "GTAGAGCTGT", "k": 3 }
}
```

**Response:**
```json
{
  "totalKmers": 8,
  "uniqueKmers": 8,
  "maxCount": 1,
  "minCount": 1,
  "averageCount": 1.0,
  "entropy": 3.0,
  "distinctKmers": 8,
  "singletonKmers": 8
}
```
All 8 trimers distinct ⇒ H = log₂(8) = 3 bits exactly.

## Performance

- **Time Complexity:** O(n) to build the k-mer table.
- **Space Complexity:** O(distinct k-mers).

## See Also

- [count_kmers](count_kmers.md) — raw k-mer → count map
- [kmer_spectrum](kmer_spectrum.md) — frequency-of-frequencies
- [kmer_frequencies](kmer_frequencies.md) — normalized frequencies
