# scan_with_pwm

Scan a DNA sequence with a Position Weight Matrix.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `scan_with_pwm` |
| **Method ID** | `MotifFinder.ScanWithPwm` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans a DNA sequence with a **4×L Position Weight Matrix** whose rows are the bases
A, C, G, T. At each offset the score is the sum of `Matrix[baseIndex, position]` over
the window; matches scoring at or above `threshold` are returned with their position,
matched subsequence, PWM consensus and score (Wikipedia PWM scoring rule).

## Core Documentation Reference

- Source: [MotifFinder.cs#L282](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs#L282)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to scan (min length 1) |
| `pwm` | object | Yes | `{ matrix: 4×L jagged (rows A,C,G,T), length: L }` |
| `threshold` | number | No | Minimum score threshold (default 0.0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | Matches `{ position, matchedSequence, pattern, score }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |
| 1001 | Matrix must have 4 rows |

## Examples

### Example 1: ATGC-scoring PWM

**User Prompt:**
> Scan "ATGCATGC" with a PWM that scores ATGC highly.

**Expected Tool Call:**
```json
{
  "tool": "scan_with_pwm",
  "arguments": { "sequence": "ATGCATGC", "pwm": { "matrix": [ [1, -1, -1, -1], [-1, -1, -1, 1], [-1, -1, 1, -1], [-1, 1, -1, -1] ], "length": 4 }, "threshold": 0.0 }
}
```

**Response:**
```json
{ "items": [ { "position": 0, "matchedSequence": "ATGC", "score": 4.0 }, { "position": 4, "matchedSequence": "ATGC", "score": 4.0 } ] }
```
Each ATGC window sums four +1 cells = 4.

### Example 2: No match below threshold

**User Prompt:**
> Scan "AAAA" with the same PWM.

**Expected Tool Call:**
```json
{
  "tool": "scan_with_pwm",
  "arguments": { "sequence": "AAAA", "pwm": { "matrix": [ [1, -1, -1, -1], [-1, -1, -1, 1], [-1, -1, 1, -1], [-1, 1, -1, -1] ], "length": 4 }, "threshold": 0.0 }
}
```

**Response:**
```json
{ "items": [] }
```
AAAA scores 1 − 1 − 1 − 1 = −2, below the threshold.

## Performance

- **Time Complexity:** O((n − L)·L).
- **Space Complexity:** O(number of matches).

## See Also

- [create_pwm](create_pwm.md) — build a PWM from an alignment
- [generate_consensus](generate_consensus.md) — IUPAC consensus
- [find_exact_motif](find_exact_motif.md) — exact motif search
