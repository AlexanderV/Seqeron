# suffix_tree_contains

Check if a pattern exists in text using suffix tree.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Core |
| **Tool Name** | `suffix_tree_contains` |
| **Method ID** | `SuffixTree.Contains` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Searches for a pattern within text using a suffix tree data structure. The suffix tree enables O(m) time complexity for pattern searching, where m is the pattern length, regardless of the text size.

## Core Documentation Reference

- Source: [SuffixTree.Search.cs#L22](../../../../src/SuffixTree/SuffixTree.Search.cs#L22)
- Algorithm: Ukkonen's suffix tree with O(n) construction and O(m) search

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `text` | string | Yes | The text to search in (min length: 1) |
| `pattern` | string | Yes | The pattern to search for |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `found` | boolean | Whether the pattern was found in the text |

## Errors

| Code | Message |
|------|---------|
| 1001 | Text cannot be null or empty |
| 1002 | Pattern cannot be null |

## Examples

### Example 1: Find substring

**User Prompt:**
> Does the text "banana" contain the substring "ana"?

**Expected Tool Call:**
```json
{
  "tool": "suffix_tree_contains",
  "arguments": {
    "text": "banana",
    "pattern": "ana"
  }
}
```

**Response:**
```json
{
  "found": true
}
```

### Example 2: Pattern not found

**User Prompt:**
> Check if "xyz" appears in "banana"

**Expected Tool Call:**
```json
{
  "tool": "suffix_tree_contains",
  "arguments": {
    "text": "banana",
    "pattern": "xyz"
  }
}
```

**Response:**
```json
{
  "found": false
}
```

## Performance

- **Time Complexity:** O(n) for tree construction + O(m) for search
- **Space Complexity:** O(n) where n is text length
- **Allocations:** Minimal (tree is built once per invocation)

## See Also

- [suffix_tree_count](suffix_tree_count.md) - Count pattern occurrences
- [suffix_tree_find_all](suffix_tree_find_all.md) - Find all occurrence positions
