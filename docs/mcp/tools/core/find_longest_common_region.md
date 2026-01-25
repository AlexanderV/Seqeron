# find_longest_common_region

Find the longest common region between two DNA sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Core |
| **Tool Name** | `find_longest_common_region` |
| **Method ID** | `GenomicAnalyzer.FindLongestCommonRegion` |
| **Version** | 1.0.0 |

## Description

Finds the longest common substring between two DNA sequences and returns positions in both. Useful for identifying conserved regions and gene homology. Time complexity: O(n+m).

## Core Documentation Reference

- Source: [GenomicAnalyzer.cs#L178](../../../../SuffixTree.Genomics/GenomicAnalyzer.cs#L178)

## Input/Output

**Input:** `{ sequence1: string, sequence2: string }`

**Output:** `{ region: string, position1: integer, position2: integer, length: integer }`

## Example

**User Prompt:**
> Find common region between ATGCATGC and CATGCAT

**Tool Call:**
```json
{ "tool": "find_longest_common_region", "arguments": { "sequence1": "ATGCATGC", "sequence2": "CATGCAT" } }
```

**Response:**
```json
{ "region": "ATGC", "position1": 0, "position2": 1, "length": 4 }
```
