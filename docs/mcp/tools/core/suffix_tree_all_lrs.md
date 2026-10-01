# suffix_tree_all_lrs

Find every longest repeated substring (all ties) with all positions.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Core |
| **Tool Name** | `suffix_tree_all_lrs` |
| **Method ID** | `SuffixTree.FindAllLongestRepeatedSubstrings` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns every distinct longest repeated substring of `text` — all length ties, unlike `suffix_tree_lrs`, which returns a single representative — each with all 0-based start positions in ascending order (occurrences may overlap). Substrings are ordered by first occurrence; the list is empty when no character repeats. Identical across the in-memory and persistent suffix trees.

## Core Documentation Reference

- Source: [SuffixTree/SuffixTree.Algorithms.cs#L67](../../../../src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Algorithms.cs#L67)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `text` | string | Yes | The text to analyze |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `substrings` | array<object> | `{substring, positions}` for every longest repeated substring, ordered by first occurrence |
| `length` | integer | Common length of the substrings (0 when no character repeats) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Text cannot be null or empty |

## Examples

### Example 1: Three tied repeats

**User Prompt:**
> List all longest repeated substrings of abcxbcaxcab.

**Tool Call:**
```json
{
  "tool": "suffix_tree_all_lrs",
  "arguments": {
    "text": "abcxbcaxcab"
  }
}
```

**Response:**
```json
{
  "substrings": [
    {
      "substring": "ab",
      "positions": [
        0,
        9
      ]
    },
    {
      "substring": "bc",
      "positions": [
        1,
        4
      ]
    },
    {
      "substring": "ca",
      "positions": [
        5,
        8
      ]
    }
  ],
  "length": 2
}
```

### Example 2: Single repeat

**Tool Call:**
```json
{
  "tool": "suffix_tree_all_lrs",
  "arguments": {
    "text": "banana"
  }
}
```

**Response:**
```json
{
  "substrings": [
    {
      "substring": "ana",
      "positions": [
        1,
        3
      ]
    }
  ],
  "length": 3
}
```

## See Also

- [suffix_tree_lrs](suffix_tree_lrs.md) — One longest repeated substring
- [find_longest_repeat](find_longest_repeat.md) — Longest repeated region of a DNA sequence

## References

- Algorithm source: [SuffixTree/SuffixTree.Algorithms.cs#L67](../../../../src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Algorithms.cs#L67)
- Binding: [SuffixTreeCoreTools.cs](../../../../src/SuffixTree/Mcp/SuffixTree.Mcp.Core/Tools/SuffixTreeCoreTools.cs)
