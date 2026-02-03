# find_longest_repeat

Find the longest repeated region in a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Core |
| **Tool Name** | `find_longest_repeat` |
| **Method ID** | `GenomicAnalyzer.FindLongestRepeat` |
| **Version** | 1.0.0 |

## Description

Finds the longest repeated region in a DNA sequence. Useful for detecting tandem repeats, transposable elements, etc. Uses suffix tree for O(n) time complexity.

## Core Documentation Reference

- Source: [GenomicAnalyzer.cs#L20](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/GenomicAnalyzer.cs#L20)

## Input/Output

**Input:** `{ sequence: string }` (DNA: A, T, G, C)

**Output:** `{ repeat: string, positions: integer[], length: integer }`

## Example

**User Prompt:**
> Find repeated sequences in ATGATGATG

**Tool Call:**
```json
{ "tool": "find_longest_repeat", "arguments": { "sequence": "ATGATGATG" } }
```

**Response:**
```json
{ "repeat": "ATG", "positions": [0, 3, 6], "length": 3 }
```
