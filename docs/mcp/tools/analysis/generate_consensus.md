# generate_consensus

IUPAC consensus sequence from aligned DNA sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `generate_consensus` |
| **Method ID** | `MotifFinder.GenerateConsensus` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Builds an **IUPAC consensus** from aligned, equal-length DNA sequences. At each
column, every base whose count strictly exceeds 25% of the sequence count is included,
and the set of included bases is mapped to its IUPAC ambiguity code (e.g. {A,T} → W,
{A,G} → R). A unanimous column yields the single base. Ties among "present" bases are
resolved by the IUPAC code for the whole set.

## Core Documentation Reference

- Source: [MotifFinder.cs#L339](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs#L339)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array of string | Yes | Aligned DNA sequences of equal length (≥ 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `consensus` | string | IUPAC consensus sequence |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least one sequence is required |

## Examples

### Example 1: Unanimous columns

**User Prompt:**
> Consensus of ["ATGC","ATGC","ATGC"].

**Expected Tool Call:**
```json
{
  "tool": "generate_consensus",
  "arguments": { "sequences": ["ATGC", "ATGC", "ATGC"] }
}
```

**Response:**
```json
{ "consensus": "ATGC" }
```

### Example 2: A/T ambiguity (→ W)

**User Prompt:**
> Consensus of ["AAAA","TTTT"].

**Expected Tool Call:**
```json
{
  "tool": "generate_consensus",
  "arguments": { "sequences": ["AAAA", "TTTT"] }
}
```

**Response:**
```json
{ "consensus": "WWWW" }
```
Each column has A and T each at 50% (> 25%), so both are included ⇒ IUPAC W.

## Performance

- **Time Complexity:** O(L · S) for S sequences of length L.
- **Space Complexity:** O(L).

## See Also

- [create_pwm](create_pwm.md) — position weight matrix from an alignment
- [scan_with_pwm](scan_with_pwm.md) — scan with a PWM
