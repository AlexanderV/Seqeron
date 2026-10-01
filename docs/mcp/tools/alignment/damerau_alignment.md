# damerau_alignment

Optimal transposition-aware edit script (Lowrance–Wagner trace) turning sequence1 into sequence2.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `damerau_alignment` |
| **Method ID** | `ApproximateMatcher.GetDamerauLevenshteinAlignment` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Runs the Lowrance & Wagner (1975) dynamic program and traces back an optimal edit script. `variant` = `unrestricted` (default; `ApproximateMatcher.GetDamerauLevenshteinAlignment`) allows a **transposition block** a_k…a_i → b_l…b_j (a_k = b_j, a_i = b_l): the i−k−1 characters between the swapped pair are deleted, the pair is swapped, and b_(l+1..j−1) are inserted between them, at cost (i−k−1)·deletionCost + transpositionCost + (j−l−1)·insertionCost. `variant` = `osa` (`ApproximateMatcher.GetOptimalStringAlignment`) allows adjacent swaps only. The distance is the summed operation cost and equals [damerau_levenshtein_distance](damerau_levenshtein_distance.md) for the same variant and costs; replaying the operations on sequence1 yields sequence2 (verified on exhaustive and random cases, unrestricted optimality against an exhaustive Dijkstra search over all edit sequences). Traceback tie-break: diagonal (match/substitution), then transposition, then deletion, then insertion. Weighted costs as in damerau_levenshtein_distance (`unrestricted` requires 2·transpositionCost ≥ insertionCost + deletionCost). Case-sensitive. O(m·n) time and space.

**Script letters** (edlib letters as in [edit_alignment](edit_alignment.md)): `=` match, `X` substitution, `I` sequence1 character deleted, `D` sequence2 character inserted, `T` a transposition block followed by one `i` per character deleted and one `d` per character inserted inside the block (an adjacent swap is `T`).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/ApproximateMatcher.Damerau.cs#L151](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.Damerau.cs#L151)
- Source (GetOptimalStringAlignment): [Seqeron.Genomics.Alignment/ApproximateMatcher.Damerau.cs#L181](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.Damerau.cs#L181)
- Algorithm doc: [Edit_Distance.md](../../../algorithms/Pattern_Matching/Edit_Distance.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | Source sequence (s1). |
| `sequence2` | string | Yes | Target sequence (s2). |
| `variant` | string | No | 'unrestricted' (true Damerau–Levenshtein, default) or 'osa' (optimal string alignment). |
| `insertionCost` | integer | No | Cost of inserting a sequence2 character (>= 0; default 1). |
| `deletionCost` | integer | No | Cost of deleting a sequence1 character (>= 0; default 1). |
| `substitutionCost` | integer | No | Cost of a substitution (>= 0; default 1). |
| `transpositionCost` | integer | No | Cost of swapping two adjacent characters (>= 0; default 1; 'unrestricted' needs 2·transpositionCost >= insertionCost + deletionCost). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `distance` | integer | Summed cost of the operations (= the distance) |
| `variant` | string | `unrestricted` or `osa` |
| `script` | string | Compact script (letters above) |
| `transpositionCount` | integer | Number of transposition blocks |
| `operations` | array | `{kind, sourcePosition, sourceLength, targetPosition, targetLength, cost}`; kind = match / substitution / insertion / deletion / transposition; 0-based positions |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1004 | Variant must be 'unrestricted' or 'osa' |
| 1004 | Costs must be >= 0; 'unrestricted' requires 2·transpositionCost >= insertionCost + deletionCost |

## Examples

### Example 1: Transposition with an inserted character

**User Prompt:**
> Show the Damerau–Levenshtein edit script from CA to ABC.

**Tool Call:**
```json
{
  "tool": "damerau_alignment",
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
  "variant": "unrestricted",
  "script": "Td",
  "transpositionCount": 1,
  "operations": [
    { "kind": "transposition", "sourcePosition": 0, "sourceLength": 2, "targetPosition": 0, "targetLength": 3, "cost": 2 }
  ]
}
```

### Example 2: Adjacent swap (OSA)

**Tool Call:**
```json
{
  "tool": "damerau_alignment",
  "arguments": {
    "sequence1": "ACGT",
    "sequence2": "AGCT",
    "variant": "osa"
  }
}
```

**Response:**
```json
{
  "distance": 1,
  "variant": "osa",
  "script": "=T=",
  "transpositionCount": 1,
  "operations": [
    { "kind": "match", "sourcePosition": 0, "sourceLength": 1, "targetPosition": 0, "targetLength": 1, "cost": 0 },
    { "kind": "transposition", "sourcePosition": 1, "sourceLength": 2, "targetPosition": 1, "targetLength": 2, "cost": 1 },
    { "kind": "match", "sourcePosition": 3, "sourceLength": 1, "targetPosition": 3, "targetLength": 1, "cost": 0 }
  ]
}
```

## Worked Example

CA → ABC: the block a_1 a_2 = "CA" with a_1 = b_3 = C and a_2 = b_1 = A is swapped to "AC" and b_2 = B is inserted between them — one transposition plus one insertion, cost 2 (script `Td`). OSA cannot edit the swapped pair again and needs 3 operations.

## See Also

- [damerau_levenshtein_distance](damerau_levenshtein_distance.md) — the distance only
- [edit_alignment](edit_alignment.md) — Levenshtein alignment (no transpositions)

## References

- Lowrance R., Wagner R.A. (1975). An extension of the string-to-string correction problem. J. ACM 22(2):177–183.
- Algorithm source: [Seqeron.Genomics.Alignment/ApproximateMatcher.Damerau.cs#L151](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.Damerau.cs#L151)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
