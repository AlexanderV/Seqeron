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

Finds the longest common substring between two DNA sequences and returns the **first (lowest) 0-based occurrence** in each. On a length tie the substring whose first match in `sequence2` ends earliest is returned. Useful for identifying conserved regions and gene homology. Time complexity: O(n+m).

## Core Documentation Reference

- Source: [GenomicAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs)

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
{ "region": "ATGCAT", "position1": 0, "position2": 1, "length": 6 }
```
