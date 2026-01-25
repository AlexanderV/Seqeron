# suffix_tree_lrs

Find the longest repeated substring in text.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Core |
| **Tool Name** | `suffix_tree_lrs` |
| **Method ID** | `SuffixTree.LongestRepeatedSubstring` |
| **Version** | 1.0.0 |

## Description

Finds the longest substring that appears at least twice in the text. Returns empty string if no repeats exist.

## Core Documentation Reference

- Source: [SuffixTree.Algorithms.cs#L14](../../../../SuffixTree/SuffixTree.Algorithms.cs#L14)

## Input/Output

**Input:** `{ text: string }`

**Output:** `{ substring: string, length: integer }`

## Example

**User Prompt:**
> What is the longest repeated substring in "banana"?

**Tool Call:**
```json
{ "tool": "suffix_tree_lrs", "arguments": { "text": "banana" } }
```

**Response:**
```json
{ "substring": "ana", "length": 3 }
```

## See Also

- [suffix_tree_lcs](suffix_tree_lcs.md) - Longest common substring between two texts
