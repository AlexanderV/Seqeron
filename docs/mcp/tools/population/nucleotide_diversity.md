# nucleotide_diversity

Compute nucleotide diversity π from aligned sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Population |
| **Tool Name** | `nucleotide_diversity` |
| **Method ID** | `PopulationGeneticsAnalyzer.CalculateNucleotideDiversity` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Averages the per-site number of differences over all `C(n, 2)` sequence pairs: `π = (Σ_{i<j} diffs(i, j)) / (C(n, 2) · L)`, where L is the alignment length. Identical sequences give π = 0; two sequences differing at every position give π = 1. Fewer than 2 sequences returns 0.

## Core Documentation Reference

- Source: [PopulationGeneticsAnalyzer.cs#L205](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs#L205)
- Algorithm doc: [Diversity_Statistics.md](../../../algorithms/Population_Genetics/Diversity_Statistics.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | string[] | Yes | Aligned sequences of equal length |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `pi` | number | Nucleotide diversity π ∈ [0, 1] |

## Errors

None (fewer than 2 sequences returns 0).

## Examples

### Example 1: Two maximally different sequences

**User Prompt:**
> π for AAAA and TTTT?

**Expected Tool Call:**
```json
{
  "tool": "nucleotide_diversity",
  "arguments": { "sequences": ["AAAA", "TTTT"] }
}
```

**Response:**
```json
{ "pi": 1.0 }
```

4 differences / (1 pair × 4 sites) = 1.0.

### Example 2: Wikipedia Tajima's D example

**User Prompt:**
> π for the 5-sequence Wikipedia dataset.

**Expected Tool Call:**
```json
{
  "tool": "nucleotide_diversity",
  "arguments": {
    "sequences": ["00000000000000000000", "00100000000010000010", "00000000000010000010", "00000010000000000010", "00000010000010000010"]
  }
}
```

**Response:**
```json
{ "pi": 0.1 }
```

Total pairwise differences 20 over C(5,2)=10 pairs and 20 sites → k̂ = 2.0, π = 2.0/20 = 0.1.

## Performance

- **Time Complexity:** O(n²·L)
- **Space Complexity:** O(n·L)

## References

- Nei, M., Li, W.-H. (1979). *PNAS* 76:5269–5273.

## See Also

- [diversity_statistics](diversity_statistics.md), [wattersons_theta](wattersons_theta.md), [tajimas_d](tajimas_d.md)
