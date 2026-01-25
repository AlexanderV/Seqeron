# suffix_tree_stats

Get statistics about a suffix tree built from text.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Core |
| **Tool Name** | `suffix_tree_stats` |
| **Method ID** | `SuffixTree.Properties` |
| **Version** | 1.0.0 |

## Description

Returns structural statistics about the suffix tree: total nodes, leaf count, maximum depth, and input text length.

## Core Documentation Reference

- Source: [SuffixTree.cs#L206](../../../../SuffixTree/SuffixTree.cs#L206)

## Input/Output

**Input:** `{ text: string }`

**Output:** `{ nodeCount: integer, leafCount: integer, maxDepth: integer, textLength: integer }`

## Example

**User Prompt:**
> What are the suffix tree statistics for "banana"?

**Tool Call:**
```json
{ "tool": "suffix_tree_stats", "arguments": { "text": "banana" } }
```

**Response:**
```json
{ "nodeCount": 10, "leafCount": 7, "maxDepth": 6, "textLength": 6 }
```
