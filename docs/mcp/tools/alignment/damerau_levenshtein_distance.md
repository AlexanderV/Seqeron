# damerau_levenshtein_distance

Damerau–Levenshtein (unrestricted) or optimal string alignment distance between two sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `damerau_levenshtein_distance` |
| **Method ID** | `ApproximateMatcher.DamerauLevenshteinDistance` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Edit distance with unit-cost insertions, deletions, substitutions and transpositions of adjacent characters. `variant` = `unrestricted` (default) is the true Damerau–Levenshtein distance (Lowrance & Wagner 1975; `ApproximateMatcher.DamerauLevenshteinDistance`), a metric in which transposed characters may be separated by further edits (DL(CA, ABC) = 2). `variant` = `osa` is the optimal string alignment (restricted Damerau) distance (`ApproximateMatcher.OptimalStringAlignmentDistance`), where no substring is edited more than once (OSA(CA, ABC) = 3). Case-sensitive (ordinal). Cross-checked against rapidfuzz and jellyfish.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L751](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L751)
- Source (OptimalStringAlignmentDistance): [Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L705](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L705)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence. |
| `sequence2` | string | Yes | Second sequence. |
| `variant` | string | No | 'unrestricted' (true Damerau–Levenshtein, default) or 'osa' (optimal string alignment). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `distance` | integer | Distance (>= 0) |
| `variant` | string | Variant used: `unrestricted` or `osa` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1004 | Variant must be 'unrestricted' or 'osa' |

## Examples

### Example 1: True Damerau–Levenshtein

**User Prompt:**
> Damerau–Levenshtein distance between CA and ABC?

**Tool Call:**
```json
{
  "tool": "damerau_levenshtein_distance",
  "arguments": {
    "sequence1": "CA",
    "sequence2": "ABC"
  }
}
```

**Response:**
```json
{
  "distance": 2,
  "variant": "unrestricted"
}
```

### Example 2: Optimal string alignment

**Tool Call:**
```json
{
  "tool": "damerau_levenshtein_distance",
  "arguments": {
    "sequence1": "CA",
    "sequence2": "ABC",
    "variant": "osa"
  }
}
```

**Response:**
```json
{
  "distance": 3,
  "variant": "osa"
}
```

## Worked Example

CA → AC (transposition) → ABC (insertion) costs 2, but edits the transposed pair again, which OSA forbids; OSA needs 3 edits.

## See Also

- `edit_distance` — Levenshtein distance (Core server)
- [edit_alignment](edit_alignment.md) — Levenshtein alignment

## References

- Algorithm source: [Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L751](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L751)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
