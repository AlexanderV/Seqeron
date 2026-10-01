# generate_dumb_consensus

Majority-threshold consensus (Biopython dumb_consensus semantics).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `generate_dumb_consensus` |
| **Method ID** | `MotifFinder.GenerateDumbConsensus` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Majority-threshold consensus with the semantics of Biopython `Bio.Align.AlignInfo.SummaryInfo.dumb_consensus(threshold, ambiguous, require_multiple)` (Biopython ≤ 1.85). Per column the residues other than the gaps `-` and `.` are counted (case-sensitive, any alphabet); if exactly one residue has the maximum count and that count divided by the number of non-gap residues is ≥ `threshold`, it is emitted; otherwise (ties, below threshold, all-gap column, or a single non-gap residue when `requireMultiple` is set) `ambiguous` is emitted. Cross-checked against Biopython 1.85 (708/708).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MotifFinder.AlignmentConsensus.cs#L454](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.AlignmentConsensus.cs#L454)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Aligned sequences of equal length |
| `threshold` | number | No | Required fraction of non-gap residues (default 0.7) |
| `ambiguous` | string | No | Single-character no-consensus symbol (default 'X') |
| `requireMultiple` | boolean | No | A column with one non-gap residue gives the ambiguous symbol (default false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `consensus` | string | Majority consensus |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least one sequence is required |
| 1002 | Ambiguous must be exactly one character |
| 1003 | Threshold cannot be NaN |
| — | Null element or unequal lengths (ArgumentException) |

## Examples

### Example 1: Threshold 0.7, ambiguous N

**User Prompt:**
> Biopython-style dumb consensus of ACGT/ATGT/ATGT with N for ambiguous columns.

**Tool Call:**
```json
{
  "tool": "generate_dumb_consensus",
  "arguments": {
    "sequences": [
      "ACGT",
      "ATGT",
      "ATGT"
    ],
    "threshold": 0.7,
    "ambiguous": "N"
  }
}
```

**Response:**
```json
{
  "consensus": "ANGT"
}
```

### Example 2: requireMultiple

**Tool Call:**
```json
{
  "tool": "generate_dumb_consensus",
  "arguments": {
    "sequences": [
      "A-",
      "-A",
      "--"
    ],
    "requireMultiple": true
  }
}
```

**Response:**
```json
{
  "consensus": "XX"
}
```

## Worked Example

Column 2 has C×1, T×2: T is 2/3 = 0.67 < 0.7 → N.

## See Also

- [generate_emboss_consensus](generate_emboss_consensus.md) — EMBOSS cons consensus
- `compute_consensus` — Majority consensus of reads (Alignment server)

## References

- Algorithm source: [Seqeron.Genomics.Analysis/MotifFinder.AlignmentConsensus.cs#L454](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.AlignmentConsensus.cs#L454)
- Binding: [AnalysisTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs)
