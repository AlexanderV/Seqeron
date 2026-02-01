# suffix_tree_lcs

Find the longest common substring between two texts.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Core |
| **Tool Name** | `suffix_tree_lcs` |
| **Method ID** | `SuffixTree.LongestCommonSubstring` |
| **Version** | 1.0.0 |

## Description

Finds the longest substring that appears in both texts. Uses suffix tree for efficient O(n+m) time complexity.

## Core Documentation Reference

- Source: [SuffixTree.Algorithms.cs#L28](../../../../src/SuffixTree/SuffixTree.Algorithms.cs#L28)

## Input/Output

**Input:** `{ text1: string, text2: string }`

**Output:** `{ substring: string, length: integer }`

## Example

**User Prompt:**
> What is the longest common substring between "banana" and "panama"?

**Tool Call:**
```json
{ "tool": "suffix_tree_lcs", "arguments": { "text1": "banana", "text2": "panama" } }
```

**Response:**
```json
{ "substring": "ana", "length": 3 }
```

## See Also

- [suffix_tree_lrs](suffix_tree_lrs.md) - Longest repeated substring in single text
