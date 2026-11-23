# find_longest_common_region

Find the longest common region between two DNA sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Core |
| **Tool Name** | `find_longest_common_region` |
| **Method ID** | `GenomicAnalyzer.FindLongestCommonRegion` |
| **Version** | 1.1.0 |

## Description

Finds the longest common substring between two DNA sequences and returns the **first (lowest) 0-based occurrence** in each. On a length tie the substring whose first match in `sequence2` ends earliest is returned. Useful for identifying conserved regions and gene homology. Time complexity: O(n+m).

The sequences may hold N / IUPAC codes and lower case (the string overloads: only A/C/G/T match, any other symbol ends a region). `strands` (`Forward` default, `ReverseComplement`, `Both`) also searches the reverse complement of `sequence2` (`isReverseComplement`, `position2` = the forward-strand start of the window that spells it). `maxMismatches > 0` returns the k-mismatch longest common substring (`GenomicAnalyzer.FindLongestCommonRegionWithMismatches`; `region` is the window of `sequence1`, `region2` that of `sequence2`, `mismatches` their Hamming distance). The tool observes the cancellation token before and after the search; the search itself is not interruptible.

## Core Documentation Reference

- Source: [GenomicAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs)

## Input/Output

**Input:** `{ sequence1: string, sequence2: string, strands?: "Forward" | "ReverseComplement" | "Both", maxMismatches?: integer }`

**Output:** `{ region: string, position1: integer, position2: integer, length: integer, isReverseComplement?: boolean, mismatches?: integer, region2?: string }`

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
