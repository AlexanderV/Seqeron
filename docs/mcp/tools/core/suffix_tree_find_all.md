# suffix_tree_find_all

Find all positions where a pattern occurs in text using suffix tree.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Core |
| **Tool Name** | `suffix_tree_find_all` |
| **Method ID** | `SuffixTree.FindAllOccurrences` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds all positions where a pattern occurs in text using a suffix tree. Returns an array of starting positions (0-indexed) for each occurrence. Overlapping occurrences are reported separately.

## Core Documentation Reference

- Source: [SuffixTree.Search.cs#L115](../../../../src/SuffixTree/SuffixTree.Search.cs#L115)
- Algorithm: Pattern matching with suffix collection

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `text` | string | Yes | The text to search in (min length: 1) |
| `pattern` | string | Yes | The pattern to find |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `positions` | integer[] | Array of starting positions (0-indexed) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Text cannot be null or empty |
| 1002 | Pattern cannot be null |

## Examples

### Example 1: Find all occurrences

**User Prompt:**
> Where does "ana" appear in "banana"?

**Expected Tool Call:**
```json
{
  "tool": "suffix_tree_find_all",
  "arguments": {
    "text": "banana",
    "pattern": "ana"
  }
}
```

**Response:**
```json
{
  "positions": [1, 3]
}
```

### Example 2: Pattern not found

**User Prompt:**
> Find "xyz" in "banana"

**Expected Tool Call:**
```json
{
  "tool": "suffix_tree_find_all",
  "arguments": {
    "text": "banana",
    "pattern": "xyz"
  }
}
```

**Response:**
```json
{
  "positions": []
}
```

## Performance

- **Time Complexity:** O(n) for tree construction + O(m + k) for finding all occurrences
- **Space Complexity:** O(n) where n is text length
- **Note:** k is the number of occurrences found

## See Also

- [suffix_tree_contains](suffix_tree_contains.md) - Check if pattern exists
- [suffix_tree_count](suffix_tree_count.md) - Count occurrences
