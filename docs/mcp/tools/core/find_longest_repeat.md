# find_longest_repeat

Find the longest repeated region in a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Core |
| **Tool Name** | `find_longest_repeat` |
| **Method ID** | `GenomicAnalyzer.FindLongestRepeat` |
| **Version** | 1.1.0 |

## Description

Finds the longest repeated region in a DNA sequence: the longest substring occurring at least twice (the deepest internal node of the suffix tree; occurrences may overlap, so `ATGATGATG` gives `ATGATG` at 0 and 3). O(n) after the suffix-tree build. When several different substrings tie for the maximum length, the one whose first occurrence is leftmost is returned (every tied substring: `suffix_tree_all_lrs`). Returns an empty `repeat` when nothing repeats. The sequence may hold N / IUPAC codes and lower case (string overload of `GenomicAnalyzer.FindLongestRepeat`: only A/C/G/T match, any other symbol ends a repeat); `bothStrands: true` reports the longest repeat class over both strands (`GenomicAnalyzer.FindAllLongestRepeatsBothStrands`, first class) with `reversePositions`.

## Core Documentation Reference

- Source: [GenomicAnalyzer.cs#L29](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs#L19)

## Input/Output

**Input:** `{ sequence: string, bothStrands?: boolean }` (nucleotides: A, C, G, T; N and IUPAC codes accepted and never repeat)

**Output:** `{ repeat: string, positions: integer[], length: integer, reversePositions?: integer[] }`

## Example

**User Prompt:**
> Find repeated sequences in ATGATGATG

**Tool Call:**
```json
{ "tool": "find_longest_repeat", "arguments": { "sequence": "ATGATGATG" } }
```

**Response:**
```json
{ "repeat": "ATGATG", "positions": [0, 3], "length": 6 }
```
