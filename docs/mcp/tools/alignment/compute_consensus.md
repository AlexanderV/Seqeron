# compute_consensus

Build a consensus from pre-aligned reads of equal length with Biopython `dumb_consensus` semantics: per column the gaps '-' and '.' are skipped (every other character, `N` included, is a residue); the unique most frequent residue is emitted if its fraction of the non-gap residues is >= `threshold` (default 0.5), otherwise the ambiguous symbol (default 'N'). Optional `requireMultiple` (Biopython `require_multiple`) and `ignoreCase` (default true; false = Biopython's exact comparison).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `compute_consensus` |
| **Method ID** | `SequenceAssembler.ComputeConsensus` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Build a consensus from pre-aligned reads of equal length with Biopython `dumb_consensus` semantics: per column the gaps '-' and '.' are skipped (every other character, `N` included, is a residue); the unique most frequent residue is emitted if its fraction of the non-gap residues is >= `threshold` (default 0.5), otherwise the ambiguous symbol (default 'N'). Optional `requireMultiple` (Biopython `require_multiple`) and `ignoreCase` (default true; false = Biopython's exact comparison).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `alignedReads` | array<string> | Yes | Pre-aligned reads of identical length. |
| `threshold` | number | No | Required fraction of the column's non-gap residues (default 0.5; Biopython default 0.7). NaN rejected. |
| `ambiguous` | string | No | One character for columns without consensus (default `N`; Biopython default `X`). |
| `requireMultiple` | boolean | No | A column with a single non-gap residue gives the ambiguous symbol (default false). |
| `ignoreCase` | boolean | No | Upper-case residues before counting (default true); false compares exactly as Biopython. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `consensus` | string | The consensus sequence (length equals the longest read). |

## Errors

| Code | Message |
|------|---------|
| 1001 | All aligned reads must have the same length |
| 1002 | Ambiguous must be exactly one character |
| 1003 | A read is null |
| 1004 | Threshold cannot be NaN |

## Examples

### Example 1: Majority vote

**Tool Call:**
```json
{
  "tool": "compute_consensus",
  "arguments": {
    "alignedReads": [
      "ACGT",
      "ACGT",
      "ACGA"
    ]
  }
}
```

**Response:**
```json
{
  "consensus": "ACGT"
}
```

### Example 2: Tie -> ambiguous

**Tool Call:**
```json
{
  "tool": "compute_consensus",
  "arguments": {
    "alignedReads": [
      "A",
      "C"
    ]
  }
}
```

**Response:**
```json
{
  "consensus": "N"
}
```

## Worked Example

For columns `A`,`C`,`G`,`{T,T,A}`, the fourth column has T at 2/3 (>= 0.5) so the consensus is `ACGT`. A 1:1 tie emits `N`. With `threshold: 0.7`, `ACGT, ATGT, ATGT` gives `ANGT` (T = 2/3 < 0.7), and `requireMultiple: true` turns `A-, AC` into `AN` — both identical to Biopython 1.85 `dumb_consensus`.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
