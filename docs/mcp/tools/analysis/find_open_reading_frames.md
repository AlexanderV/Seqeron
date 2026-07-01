# find_open_reading_frames

Open reading frames (ORFs) in all six reading frames.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_open_reading_frames` |
| **Method ID** | `GenomicAnalyzer.FindOpenReadingFrames` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds **open reading frames** in all six frames (3 forward + 3 reverse-complement).
Each ORF runs from a start codon `ATG` to its first in-frame stop codon
(`TAA`/`TAG`/`TGA`), inclusive, and is reported only if its nucleotide length is at
least `minLength`. Every ATG is considered independently, so nested ORFs sharing a
stop are all reported (canonical Rosalind ORF semantics). The reported sequence length
is divisible by 3.

## Core Documentation Reference

- Source: [GenomicAnalyzer.cs#L449](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs#L449)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `minLength` | integer | No | Minimum ORF length in nucleotides (default 100, ≥ 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | ORFs: `{ sequence, position, frame, isReverseComplement, length, codonCount }` |

`frame` is 1–3; `position` is 0-based in the scanned strand; `codonCount` = length / 3.

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |
| 1002 | Minimum ORF length must be at least 1 |

## Examples

### Example 1: Minimal forward ORF (ATG-TAA)

**User Prompt:**
> Find ORFs of length ≥ 6 in "ATGTAA".

**Expected Tool Call:**
```json
{
  "tool": "find_open_reading_frames",
  "arguments": { "sequence": "ATGTAA", "minLength": 6 }
}
```

**Response:**
```json
{ "items": [ { "sequence": "ATGTAA", "position": 0, "frame": 1, "isReverseComplement": false, "length": 6, "codonCount": 2 } ] }
```

### Example 2: One extra codon (ATG-AAA-TAG)

**User Prompt:**
> ORFs of length ≥ 9 in "ATGAAATAG".

**Expected Tool Call:**
```json
{
  "tool": "find_open_reading_frames",
  "arguments": { "sequence": "ATGAAATAG", "minLength": 9 }
}
```

**Response:**
```json
{ "items": [ { "sequence": "ATGAAATAG", "position": 0, "frame": 1, "isReverseComplement": false, "length": 9, "codonCount": 3 } ] }
```

## Performance

- **Time Complexity:** O(n) per frame (6 frames).
- **Space Complexity:** O(number of ORFs).

## See Also

- [find_motif](find_motif.md) — locate a known motif
- [codon_frequencies](codon_frequencies.md) — codon usage in a reading frame
