# find_longest_common_region_mismatches

Longest common substring with at most k mismatches (k-mismatch LCS) of two DNA sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_longest_common_region_mismatches` |
| **Method ID** | `GenomicAnalyzer.FindLongestCommonRegionWithMismatches` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Longest common substring with at most k mismatches (k-mismatch LCS, Flouri et al. 2015; Thankachan et al. 2016): the longest pair of equal-length windows, one of each sequence, that differ in at most maxMismatches positions (Hamming, no gaps); the window may start or end on a mismatch. Ties: smallest start in sequence2, then smallest in sequence1; with maxMismatches = 0 this is find_longest_common_region. strands Forward (default) | ReverseComplement | Both compares sequence1 with the reverse complement of sequence2 (positionInSecond is the forward-strand start; sequence2 is the window of the reverse-complemented sequence2); Both prefers the forward strand on a tie. A symbol other than A/C/G/T (N included) is a mismatch against everything. O(n*m) time. Length 0 and empty strings when no window exists.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/GenomicAnalyzer.CommonRegions.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.CommonRegions.cs)
- Algorithm doc: [Common_Region_Detection.md](../../../algorithms/Sequence_Comparison/Common_Region_Detection.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First nucleotide sequence: A/C/G/T plus IUPAC codes such as N (case-insensitive). |
| `sequence2` | string | Yes | Second nucleotide sequence. |
| `maxMismatches` | integer | Yes | Largest number of mismatching positions in the window pair (>= 0). |
| `strands` | string | No | Which strands of sequence2 (case-insensitive): Forward (default) \| ReverseComplement \| Both. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sequence1` | string | Window of sequence1 |
| `sequence2` | string | Equal-length window of sequence2 (or of its reverse complement) |
| `positionInFirst` | integer | 0-based start in sequence1 |
| `positionInSecond` | integer | 0-based forward-strand start in sequence2 |
| `length` | integer | Window length (0 when none) |
| `mismatches` | integer | Mismatching positions in the window pair |
| `isReverseComplement` | boolean | The window pair is against the reverse complement |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence (A/C/G/T and IUPAC codes only) |
| 1002 | maxMismatches cannot be negative |

## Examples

### Example 1: One mismatch

**User Prompt:**
> Longest common region of CACAAAAC and AACCCCAACAGCAC with at most one mismatch.

**Expected Tool Call:**
```json
{"tool": "find_longest_common_region_mismatches", "arguments": {"sequence1": "CACAAAAC", "sequence2": "AACCCCAACAGCAC", "maxMismatches": 1}}
```

**Response:**
```json
{"sequence1": "CACAA", "sequence2": "CCCAA", "positionInFirst": 0, "positionInSecond": 3, "length": 5, "mismatches": 1, "isReverseComplement": false}
```

Brute-force oracle row `5:0:3:1`.

## See Also

- [find_common_regions](find_common_regions.md) - exact common regions (modes, strands)
- [find_longest_common_region](../core/find_longest_common_region.md) - exact LCS
