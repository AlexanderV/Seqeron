# suffix_tree_count

Count the number of occurrences of a pattern in text using suffix tree.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Core |
| **Tool Name** | `suffix_tree_count` |
| **Method ID** | `SuffixTree.CountOccurrences` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Counts all occurrences of a pattern within text using a suffix tree. The algorithm uses precomputed leaf counts for O(m) time complexity where m is the pattern length. Overlapping occurrences are counted separately.

## Core Documentation Reference

- Source: [SuffixTree.Search.cs#L132](../../../../src/SuffixTree/SuffixTree.Search.cs#L132)
- Algorithm: Pattern matching with leaf count aggregation

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `text` | string | Yes | The text to search in (min length: 1) |
| `pattern` | string | Yes | The pattern to count |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `count` | integer | Number of pattern occurrences (>= 0) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Text cannot be null or empty |
| 1002 | Pattern cannot be null |

## Examples

### Example 1: Count overlapping occurrences

**User Prompt:**
> How many times does "ana" appear in "banana"?

**Expected Tool Call:**
```json
{
  "tool": "suffix_tree_count",
  "arguments": {
    "text": "banana",
    "pattern": "ana"
  }
}
```

**Response:**
```json
{
  "count": 2
}
```

### Example 2: Pattern not found

**User Prompt:**
> Count occurrences of "xyz" in "banana"

**Expected Tool Call:**
```json
{
  "tool": "suffix_tree_count",
  "arguments": {
    "text": "banana",
    "pattern": "xyz"
  }
}
```

**Response:**
```json
{
  "count": 0
}
```

## Performance

- **Time Complexity:** O(n) for tree construction + O(m) for counting
- **Space Complexity:** O(n) where n is text length
- **Note:** Uses precomputed leaf counts for instant retrieval

## See Also

- [suffix_tree_contains](suffix_tree_contains.md) - Check if pattern exists
- [suffix_tree_find_all](suffix_tree_find_all.md) - Find all occurrence positions
