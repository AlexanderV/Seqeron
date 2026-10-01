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
column, every base whose count strictly exceeds `inclusionThreshold` (default 25%) of the sequence count is included,
and the set of included bases is mapped to its IUPAC ambiguity code (e.g. {A,T} → W,
{A,G} → R). A unanimous column yields the single base. Ties among "present" bases are
resolved by the IUPAC code for the whole set. When no base exceeds the threshold (e.g. four
equally frequent bases), the column is the IUPAC code of the bases tied at the maximum
count (four-way tie ⇒ `N`); a column with no A/C/G/T (only gaps/`N`) ⇒ `N`.
`inclusionThreshold` = 0.25 calls `MotifFinder.GenerateConsensus(sequences)`; any other value in
[0, 1] calls the overload `MotifFinder.GenerateConsensus(sequences, inclusionThreshold)` (same
rule, bit-identical at 0.25). This 25% default is not a published rule — see
[generate_cavener_consensus](generate_cavener_consensus.md) (Cavener 1987) and
[generate_decipher_consensus](generate_decipher_consensus.md) (DECIPHER `ConsensusSequence`).

## Core Documentation Reference

- Source: [MotifFinder.cs#L441](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs#L441); threshold overload: [MotifFinder.AlignmentConsensus.cs#L532](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.AlignmentConsensus.cs#L532)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array of string | Yes | Aligned DNA sequences of equal length (≥ 1) |
| `inclusionThreshold` | number | No | Per-base inclusion threshold in [0, 1]: a base is included when its count > threshold × n (default 0.25) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `consensus` | string | IUPAC consensus sequence |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least one sequence is required |
| — | inclusionThreshold must be in [0, 1] (ArgumentOutOfRangeException) |
| — | Sequences cannot contain null elements (ArgumentException) |
| — | All sequences must have the same length (ArgumentException) |

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

### Example 3: Custom inclusion threshold

**Tool Call:**
```json
{
  "tool": "generate_consensus",
  "arguments": { "sequences": ["AAAA", "TTTT", "TTTT", "CCCC"], "inclusionThreshold": 0.2 }
}
```

**Response:**
```json
{ "consensus": "HHHH" }
```
A and C (25%) and T (50%) are all > 20% ⇒ {A,C,T} = H; with the default 0.25 only T passes ⇒ `TTTT`.

## Performance

- **Time Complexity:** O(L · S) for S sequences of length L.
- **Space Complexity:** O(L).

## See Also

- [generate_decipher_consensus](generate_decipher_consensus.md) — DECIPHER ConsensusSequence
- [create_pwm](create_pwm.md) — position weight matrix from an alignment
- [scan_with_pwm](scan_with_pwm.md) — scan with a PWM
