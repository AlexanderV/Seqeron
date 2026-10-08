# damerau_levenshtein_distance

Damerau–Levenshtein (unrestricted) or optimal string alignment distance between two sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `damerau_levenshtein_distance` |
| **Method ID** | `ApproximateMatcher.DamerauLevenshteinDistance` |
| **Version** | 1.1.0 |
| **Stability** | Stable |

## Description

Edit distance with insertions, deletions, substitutions and transpositions of adjacent characters (unit costs by default). `variant` = `unrestricted` (default) is the true Damerau–Levenshtein distance (Lowrance & Wagner 1975; `ApproximateMatcher.DamerauLevenshteinDistance`), a metric in which transposed characters may be separated by further edits (DL(CA, ABC) = 2). `variant` = `osa` is the optimal string alignment (restricted Damerau) distance (`ApproximateMatcher.OptimalStringAlignmentDistance`), where no substring is edited more than once (OSA(CA, ABC) = 3). Case-sensitive (ordinal). Unit costs cross-checked against rapidfuzz and jellyfish.

Optional `insertionCost` / `deletionCost` / `substitutionCost` / `transpositionCost` (W_I, W_D, W_C, W_S; default 1) give the weighted distance of Lowrance & Wagner (1975) — sequence1 → sequence2, an insertion adds a sequence2 character, a deletion removes a sequence1 character. For `unrestricted` the Lowrance–Wagner recurrence is exact only when 2·W_S ≥ W_I + W_D, so cheaper transpositions are rejected (use `osa`). Weighted `osa` equals R stringdist `method='osa'`; weighted `unrestricted` equals the exhaustive minimum over all edit sequences (R stringdist `method='dl'` agrees only when W_D = W_I = W_S). The edit script itself: [damerau_alignment](damerau_alignment.md).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L1016](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L1016)
- Source (OptimalStringAlignmentDistance): [Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L961](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L961)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence. |
| `sequence2` | string | Yes | Second sequence. |
| `variant` | string | No | 'unrestricted' (true Damerau–Levenshtein, default) or 'osa' (optimal string alignment). |
| `insertionCost` | integer | No | Cost of inserting a sequence2 character (>= 0; default 1). |
| `deletionCost` | integer | No | Cost of deleting a sequence1 character (>= 0; default 1). |
| `substitutionCost` | integer | No | Cost of a substitution (>= 0; default 1). |
| `transpositionCost` | integer | No | Cost of swapping two adjacent characters (>= 0; default 1; 'unrestricted' needs 2·transpositionCost >= insertionCost + deletionCost). |

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
| 1004 | Costs must be >= 0; 'unrestricted' requires 2·transpositionCost >= insertionCost + deletionCost |

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

### Example 3: Weighted costs

**Tool Call:**
```json
{
  "tool": "damerau_levenshtein_distance",
  "arguments": {
    "sequence1": "CA",
    "sequence2": "ABC",
    "insertionCost": 2,
    "deletionCost": 3,
    "substitutionCost": 4,
    "transpositionCost": 3
  }
}
```

**Response:**
```json
{
  "distance": 5,
  "variant": "unrestricted"
}
```

(one transposition block: swap CA → AC (3) with B inserted between (2); `osa` gives 7.)

## Worked Example

CA → AC (transposition) → ABC (insertion) costs 2, but edits the transposed pair again, which OSA forbids; OSA needs 3 edits.

## See Also

- `edit_distance` — Levenshtein distance (Core server)
- [edit_alignment](edit_alignment.md) — Levenshtein alignment
- [damerau_alignment](damerau_alignment.md) — transposition-aware edit script

## References

- Algorithm source: [Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L1016](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L1016)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
