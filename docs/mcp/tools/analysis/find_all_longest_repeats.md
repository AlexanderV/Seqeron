# find_all_longest_repeats

Every longest repeated substring of a DNA sequence (all length ties), optionally over both strands.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_all_longest_repeats` |
| **Method ID** | `GenomicAnalyzer.FindAllLongestRepeats` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EVERY longest repeated substring of a DNA sequence (all length ties; find_longest_repeat returns the leftmost one) with all its 0-based positions (overlapping occurrences counted, ascending), the deepest internal nodes of the suffix tree (Gusfield 1997 sec. 7.1); ordered by first occurrence. bothStrands=true: the longest repeat classes over both strands (a string and its reverse complement are one class, the smaller string is reported; reversePositions = forward-strand starts of the windows that spell its reverse complement), ordered by string. IUPAC input is accepted: only A/C/G/T can repeat, an N run is never a repeat. Items is empty (length 0) when nothing repeats.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/GenomicAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs)
- Algorithm doc: [Repeat_Detection.md](../../../algorithms/Repeat_Analysis/Repeat_Detection.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence: A/C/G/T plus IUPAC codes such as N (case-insensitive); non-A/C/G/T symbols never repeat. |
| `bothStrands` | boolean | No | Search both strands (default false). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | Every longest repeat: `{ sequence, positions[], length, count, reversePositions[]? }`, ordered by first occurrence (by string with `bothStrands`) |
| `length` | integer | Their common length (0 when nothing repeats) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence (A/C/G/T and IUPAC codes only) |

## Examples

### Example 1: Tandem ATGATG

**User Prompt:**
> All longest repeats of "ATGATG".

**Expected Tool Call:**
```json
{"tool": "find_all_longest_repeats", "arguments": {"sequence": "ATGATG"}}
```

**Response:**
```json
{"items": [{"sequence": "ATG", "positions": [0, 3], "length": 3, "count": 2}], "length": 3}
```

### Example 2: Both strands

**User Prompt:**
> Longest repeat classes of ACGGTTCAAAAGAACCGT on both strands.

**Expected Tool Call:**
```json
{"tool": "find_all_longest_repeats", "arguments": {"sequence": "ACGGTTCAAAAGAACCGT", "bothStrands": true}}
```

**Response:**
```json
{"items": [{"sequence": "ACGGTTC", "positions": [0], "length": 7, "count": 2, "reversePositions": [11]}], "length": 7}
```

`ACGGTTC` at 0 reappears as its reverse complement `GAACCGT` at 11.

## See Also

- [find_repeats](find_repeats.md) - every repeated substring / maximal repeats
- [find_longest_repeat](../core/find_longest_repeat.md) - one longest repeat (leftmost)
