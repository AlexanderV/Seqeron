# edit_alignment

Optimal global Levenshtein alignment (edit script, CIGAR, gapped strings) of a query against a target.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `edit_alignment` |
| **Method ID** | `ApproximateMatcher.GetEditAlignment` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Optimal global unit-cost (Levenshtein) alignment of `query` (rows) against `target` (columns) in edlib's convention: `=` match, `X` mismatch, `I` query character absent from the target, `D` target character absent from the query. Returns the distance (= edit distance), the per-column operations, the edlib EXTENDED CIGAR and STANDARD CIGAR (`M` for `=`/`X`), the gapped strings, the query-relative substitution positions and whether the path has indels. Case-sensitive. Default: full Wagner–Fischer matrix with a diagonal-first traceback (O(m·n) time and space; `ApproximateMatcher.GetEditAlignment`). With `linearSpace` = true, Hirschberg's divide-and-conquer algorithm (`ApproximateMatcher.GetEditAlignmentLinearSpace`, O(m+n) space) — same optimal distance, possibly a different co-optimal path.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L403](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L403)
- Source (GetEditAlignmentLinearSpace): [Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L446](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L446)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `query` | string | Yes | Query sequence (alignment rows). |
| `target` | string | Yes | Target sequence (alignment columns). |
| `linearSpace` | boolean | No | Use Hirschberg's linear-space algorithm instead of the full-matrix traceback (default false). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `distance` | integer | Edit distance (number of non-`=` operations) |
| `operations` | string | One of `=`, `X`, `I`, `D` per alignment column |
| `cigar` | string | edlib EXTENDED CIGAR, e.g. `3=1X1=1D1=` |
| `standardCigar` | string | edlib STANDARD CIGAR, e.g. `5M1D1M` |
| `alignedQuery` | string | Query with `-` in `D` columns |
| `alignedTarget` | string | Target with `-` in `I` columns |
| `substitutionPositions` | array<integer> | 0-based query indices of `X` columns |
| `hasIndels` | boolean | True when the path contains `I` or `D` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Query cannot be null or empty |
| 1003 | Target cannot be null or empty |

## Examples

### Example 1: survey vs surgery (edlib)

**User Prompt:**
> Align 'survey' to 'surgery' with unit edit costs.

**Tool Call:**
```json
{
  "tool": "edit_alignment",
  "arguments": {
    "query": "survey",
    "target": "surgery"
  }
}
```

**Response:**
```json
{
  "distance": 2,
  "operations": "===X=D=",
  "cigar": "3=1X1=1D1=",
  "standardCigar": "5M1D1M",
  "alignedQuery": "surve-y",
  "alignedTarget": "surgery",
  "substitutionPositions": [
    3
  ],
  "hasIndels": true
}
```

### Example 2: Linear space (Hirschberg)

**Tool Call:**
```json
{
  "tool": "edit_alignment",
  "arguments": {
    "query": "kitten",
    "target": "sitting",
    "linearSpace": true
  }
}
```

**Response:**
```json
{
  "distance": 3,
  "operations": "X===X=D",
  "cigar": "1X3=1X1=1D",
  "standardCigar": "6M1D",
  "alignedQuery": "kitten-",
  "alignedTarget": "sitting",
  "substitutionPositions": [
    0,
    4
  ],
  "hasIndels": true
}
```

## See Also

- [find_with_edits](find_with_edits.md) — Approximate matches with their alignments
- [damerau_levenshtein_distance](damerau_levenshtein_distance.md) — Edit distance with transpositions

## References

- Algorithm source: [Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L403](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L403)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
