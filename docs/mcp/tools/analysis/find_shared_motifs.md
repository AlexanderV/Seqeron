# find_shared_motifs

k-mers shared across multiple DNA sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_shared_motifs` |
| **Method ID** | `MotifFinder.FindSharedMotifs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds **k-mers present in at least `minSequences`** of the input DNA sequences. For
each shared k-mer it reports the indices of the sequences that contain it and the
prevalence (fraction of input sequences containing it). A k-mer is counted once per
sequence regardless of internal repetition.

## Core Documentation Reference

- Source: [MotifFinder.cs#L585](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs#L585)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array of string | Yes | DNA sequences (≥ 1) |
| `k` | integer | No | k-mer length (default 6, ≥ 1) |
| `minSequences` | integer | No | Minimum sequences containing the motif (default 2, ≥ 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | `{ sequence, sequenceIndices[], prevalence }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least one sequence is required |
| 1001 | Invalid DNA sequence |

## Examples

### Example 1: Shared 4-mer (ATGC)

**User Prompt:**
> Which 4-mers are shared by "ATGCATGC" and "ATGCTTTT"?

**Expected Tool Call:**
```json
{
  "tool": "find_shared_motifs",
  "arguments": { "sequences": ["ATGCATGC", "ATGCTTTT"], "k": 4, "minSequences": 2 }
}
```

**Response:**
```json
{ "items": [ { "sequence": "ATGC", "sequenceIndices": [0, 1], "prevalence": 1.0 } ] }
```
ATGC appears in both sequences (prevalence 2/2 = 1.0).

### Example 2: No shared motif

**User Prompt:**
> Shared 4-mers of "AAAA" and "TTTT"?

**Expected Tool Call:**
```json
{
  "tool": "find_shared_motifs",
  "arguments": { "sequences": ["AAAA", "TTTT"], "k": 4, "minSequences": 2 }
}
```

**Response:**
```json
{ "items": [] }
```

## Performance

- **Time Complexity:** O(Σ nᵢ) over all sequences.
- **Space Complexity:** O(distinct k-mers).

## See Also

- [find_known_motifs](find_known_motifs.md) — search a known motif set in one sequence
- [discover_motifs](discover_motifs.md) — de novo motif discovery
