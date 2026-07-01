# create_pwm

Build a log-odds Position Weight Matrix from aligned DNA sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `create_pwm` |
| **Method ID** | `MotifFinder.CreatePwm` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Builds a log-odds Position Weight Matrix (PWM) from aligned, equal-length DNA sequences.
The matrix is 4×L with rows in order A, C, G, T. Each cell is
`log2( (count + pseudocount) / (N + 4·pseudocount) / 0.25 )`, i.e. the log-odds of the
smoothed observed frequency against a uniform 0.25 background. The consensus is the
per-column argmax base; `maxScore`/`minScore` are the best/worst achievable total scores.
All input sequences must be equal length and contain only A, C, G, T.

## Core Documentation Reference

- Source: [MotifFinder.cs#L213](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs#L213)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | string[] | Yes | Aligned DNA sequences of equal length |
| `pseudocount` | number | No | Smoothing pseudocount (default 0.25, ≥ 0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `matrix` | number[][] | Jagged 4×L log-odds matrix (rows A,C,G,T) |
| `length` | integer | Motif length L |
| `consensus` | string | Per-column argmax consensus |
| `maxScore` | number | Maximum achievable total score |
| `minScore` | number | Minimum achievable total score |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least one sequence is required |
| 1003 | Pseudocount must be non-negative |
| 1002 | All sequences must have the same length |
| 1002 | Invalid character … Only A, C, G, T are valid |

## Examples

### Example 1: Two identical motifs

**Input:** `{ "sequences": ["ACGT", "ACGT"], "pseudocount": 0.25 }`

At every column the present base has frequency `2.25/3 = 0.75` and score `log2(3) ≈ 1.585`;
absent bases score `-log2(3) ≈ -1.585`.

**Response (abridged):**
```json
{ "length": 4, "consensus": "ACGT", "maxScore": 6.3399, "minScore": -6.3399 }
```

### Example 2: Consensus of a variable column

**Input:** `{ "sequences": ["AAAA", "AAAT", "AAAA"] }`
→ Column 3 has A:2, T:1, so the consensus is `"AAAA"`.

## Performance

- **Time Complexity:** O(N·L). **Space Complexity:** O(L).

## See Also

- [scan_with_pwm](scan_with_pwm.md)
- [generate_consensus](generate_consensus.md)
