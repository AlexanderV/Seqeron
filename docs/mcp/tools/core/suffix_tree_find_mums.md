# suffix_tree_find_mums

Find maximal unique matches (MUMs) between a reference text and a query.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Core |
| **Tool Name** | `suffix_tree_find_mums` |
| **Method ID** | `SuffixTree.FindMaximalUniqueMatches` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds maximal unique matches (MUMs) of length ≥ `minLength`: maximal exact matches whose string occurs exactly once in the reference `text` and — with `uniqueness` = `both` (default) — exactly once in `query` (MUMmer 3 `mummer -mum`). With `uniqueness` = `reference` only reference uniqueness is required (`mummer -mumreference`, MUMmer's default MUM-candidates). Forward strand, 0-based, sorted by query position then text position. Cross-checked against the MUMmer 3.23 binary (one documented difference: MUMmer's `-mum -l 1` sweep drops a length-1 MUM at reference position 0; the definition is kept).

## Core Documentation Reference

- Source: [SuffixTree/SuffixTree.Algorithms.cs#L213](../../../../src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Algorithms.cs#L213)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `text` | string | Yes | The reference text (the suffix tree is built on it) |
| `query` | string | Yes | The query string matched against the reference |
| `minLength` | integer | No | Minimum match length (>= 1; MUMmer default 20) |
| `uniqueness` | string | No | 'both' (unique in reference and query, -mum; default) or 'reference' (unique in reference only, -mumreference) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `matches` | array<object> | MUMs `{positionInText, positionInQuery, length}`, sorted by query position then text position |

## Errors

| Code | Message |
|------|---------|
| 1001 | Text cannot be null or empty |
| 1002 | Query cannot be null or empty |
| 1003 | minLength must be >= 1 |
| 1004 | Uniqueness must be 'both' or 'reference' |

## Examples

### Example 1: MUMmer -mum -l 4

**Tool Call:**
```json
{
  "tool": "suffix_tree_find_mums",
  "arguments": {
    "text": "TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCG",
    "query": "TTTCCTCATGCAAATAGCCTCATGCATCCGTAATGTAG",
    "minLength": 4,
    "uniqueness": "both"
  }
}
```

**Response:**
```json
{
  "matches": [
    {
      "positionInText": 0,
      "positionInQuery": 0,
      "length": 13
    },
    {
      "positionInText": 15,
      "positionInQuery": 10,
      "length": 4
    },
    {
      "positionInText": 25,
      "positionInQuery": 26,
      "length": 12
    }
  ]
}
```

### Example 2: MUMmer -mumreference -l 4

**Tool Call:**
```json
{
  "tool": "suffix_tree_find_mums",
  "arguments": {
    "text": "TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCG",
    "query": "TTTCCTCATGCAAATAGCCTCATGCATCCGTAATGTAG",
    "minLength": 4,
    "uniqueness": "reference"
  }
}
```

**Response:**
```json
{
  "matches": [
    {
      "positionInText": 0,
      "positionInQuery": 0,
      "length": 13
    },
    {
      "positionInText": 15,
      "positionInQuery": 10,
      "length": 4
    },
    {
      "positionInText": 3,
      "positionInQuery": 17,
      "length": 9
    },
    {
      "positionInText": 25,
      "positionInQuery": 26,
      "length": 12
    }
  ]
}
```

## Worked Example

The match of length 9 at (3, 17) is unique in the reference but its string `CCTCATGCA` occurs twice in the query, so it is a `-mumreference` match but not a `-mum` match.

## See Also

- [suffix_tree_find_mems](suffix_tree_find_mems.md) — All maximal exact matches (MEMs)

## References

- Algorithm source: [SuffixTree/SuffixTree.Algorithms.cs#L213](../../../../src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Algorithms.cs#L213)
- Binding: [SuffixTreeCoreTools.cs](../../../../src/SuffixTree/Mcp/SuffixTree.Mcp.Core/Tools/SuffixTreeCoreTools.cs)
