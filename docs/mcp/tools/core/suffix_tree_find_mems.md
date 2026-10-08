# suffix_tree_find_mems

Find all maximal exact matches (MEMs) between a reference text and a query.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Core |
| **Tool Name** | `suffix_tree_find_mems` |
| **Method ID** | `SuffixTree.FindMaximalExactMatches` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds every maximal exact match (MEM) of length ≥ `minLength` between the reference `text` (on which the suffix tree is built) and `query`: each (positionInText, positionInQuery, length) that is left-maximal (the strings differ, or one ends, just before the match) and right-maximal (they differ, or one ends, just after it), with every reference occurrence. Identical to MUMmer 3 `mummer -maxmatch -l minLength` on the forward strand (cross-checked against the MUMmer 3.23 binary). Positions are 0-based; results are sorted by query position, then text position. Linear-time suffix-link streaming (findmaxmat.c scheme).

## Core Documentation Reference

- Source: [SuffixTree/SuffixTree.Algorithms.cs#L205](../../../../src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Algorithms.cs#L205)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `text` | string | Yes | The reference text (the suffix tree is built on it) |
| `query` | string | Yes | The query string matched against the reference |
| `minLength` | integer | No | Minimum match length (>= 1; MUMmer default 20) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `matches` | array<object> | MEMs `{positionInText, positionInQuery, length}`, sorted by query position then text position |

## Errors

| Code | Message |
|------|---------|
| 1001 | Text cannot be null or empty |
| 1002 | Query cannot be null or empty |
| 1003 | minLength must be >= 1 |

## Examples

### Example 1: MUMmer -maxmatch -l 3

**User Prompt:**
> What are the maximal exact matches of length ≥ 3 between GATTACAGATTACA and ATTACAGGATTACAT?

**Tool Call:**
```json
{
  "tool": "suffix_tree_find_mems",
  "arguments": {
    "text": "GATTACAGATTACA",
    "query": "ATTACAGGATTACAT",
    "minLength": 3
  }
}
```

**Response:**
```json
{
  "matches": [
    {
      "positionInText": 1,
      "positionInQuery": 0,
      "length": 7
    },
    {
      "positionInText": 8,
      "positionInQuery": 0,
      "length": 6
    },
    {
      "positionInText": 0,
      "positionInQuery": 7,
      "length": 7
    },
    {
      "positionInText": 7,
      "positionInQuery": 7,
      "length": 7
    }
  ]
}
```

## Worked Example

MUMmer 3.23 `mummer -maxmatch -l 3` prints `2 1 7 / 9 1 6 / 1 8 7 / 8 8 7` (1-based reference, query, length); converted to 0-based this is exactly the response above.

## See Also

- [suffix_tree_find_mums](suffix_tree_find_mums.md) — Maximal unique matches (MUMs)
- [suffix_tree_lcs](suffix_tree_lcs.md) — Longest common substring

## References

- Algorithm source: [SuffixTree/SuffixTree.Algorithms.cs#L205](../../../../src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Algorithms.cs#L205)
- Binding: [SuffixTreeCoreTools.cs](../../../../src/SuffixTree/Mcp/SuffixTree.Mcp.Core/Tools/SuffixTreeCoreTools.cs)
