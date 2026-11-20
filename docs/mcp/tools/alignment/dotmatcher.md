# dotmatcher

EMBOSS dotmatcher: dot plot segments of two sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `dotmatcher` |
| **Method ID** | `DotMatcher.Compute` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 dotmatcher (600 DNA / protein pairs equal): the dot plot of two sequences by windowed substitution-matrix scores. On every diagonal a window of windowSize residue pairs slides along; a maximal run of windows whose summed matrix score reaches threshold is one segment from the first residue of its first window to the last residue of its last window. Returns the segments in the order the program draws them, 1-based inclusive positions on both sequences (the program's own scan treats the SECOND sequence as its first: the first diagonals start on the first residue of sequence2). matrix = BLOSUM62 | PAM250 | EDNAFULL | any EMBOSS matrix name; default EDNAFULL when sequence1 consists of nucleotide codes, EBLOSUM62 otherwise; a character that is not a letter of the matrix scores 0. Sequences not longer than the window give no segments. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program). The plot is refused with an ArgumentException above 4 194 304 segments (SequenceAligner.DotMatcherMaxSegments; a threshold of 0 or less over a long pair): raise the threshold or the window.

## Core Documentation Reference

- Source: [DotMatcher.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/DotMatcher.cs)
- Algorithm doc: [Local_Alignment_Smith_Waterman.md](../../../algorithms/Alignment/Local_Alignment_Smith_Waterman.md)
- Oracle: [generate_dotmatcher_oracle.py](../../../../scripts/alignment/generate_dotmatcher_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence (the plot's y axis in the program). |
| `sequence2` | string | Yes | Second sequence (the plot's x axis). |
| `windowSize` | integer | No | Window size (dotmatcher -windowsize, default 10). Default: `10`. |
| `threshold` | integer | No | Score a window must reach (dotmatcher -threshold, default 23). Default: `23`. |
| `matrix` | string | No | Scoring matrix name (BLOSUM62, EDNAFULL, ...); default by sequence type. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `segments` | array | Segments `{ start1, start2, end1, end2 }` (1-based inclusive on both sequences) in the order the program draws them. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Unknown matrix '...' |

## Examples

### Example 1: Two similar DNA sequences with an insertion

**User Prompt:**
> Dot plot of ACGTACGTAACCGGTTACGTAAACCCGGGTTT against ACGTACGGTAACCGGTTACGTAAACCCGGGTT, window 6, threshold 20.

**Expected Tool Call:**
```json
{"tool": "dotmatcher", "arguments": {"sequence1": "ACGTACGTAACCGGTTACGTAAACCCGGGTTT", "sequence2": "ACGTACGGTAACCGGTTACGTAAACCCGGGTT", "windowSize": 6, "threshold": 20}}
```

**Response:**
```json
{"segments": [{"start1": 1, "start2": 1, "end1": 8, "end2": 8}, {"start1": 5, "start2": 1, "end1": 10, "end2": 6}, {"start1": 11, "start2": 5, "end1": 18, "end2": 12}, {"start1": 18, "start2": 7, "end1": 24, "end2": 13}, {"start1": 21, "start2": 9, "end1": 27, "end2": 15}, {"start1": 22, "start2": 9, "end1": 28, "end2": 15}, {"start1": 24, "start2": 11, "end1": 30, "end2": 17}, {"start1": 25, "start2": 11, "end1": 31, "end2": 17}, {"start1": 17, "start2": 1, "end1": 22, "end2": 6}, {"start1": 1, "start2": 1, "end1": 8, "end2": 8}, {"start1": 5, "start2": 6, "end1": 31, "end2": 32}, {"start1": 2, "start2": 15, "end1": 12, "end2": 25}, {"start1": 8, "start2": 22, "end1": 14, "end2": 28}, {"start1": 8, "start2": 23, "end1": 14, "end2": 29}, {"start1": 10, "start2": 25, "end1": 16, "end2": 31}, {"start1": 10, "start2": 26, "end1": 16, "end2": 32}, {"start1": 1, "start2": 18, "end1": 6, "end2": 23}]}
```

dotmatcher -windowsize 6 -threshold 20 -graph data draws the same 17 lines (x on the second sequence, y on the first).

## See Also

- [local_align_top_n](local_align_top_n.md) - non-intersecting local alignments
- [find_common_regions](../analysis/find_common_regions.md) - exact common regions
