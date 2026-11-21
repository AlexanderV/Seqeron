# mummer_gaps

MUMmer 3.23 gaps: the longest consistent set of matches, then the other matches, in mgaps format.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `mummer_gaps` |
| **Method ID** | `MummerGaps.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MUMmer 3.23 gaps: reads a list of matches between two strings (start1 start2 length, as mummer prints them without the header lines; text that is not a number is refused - the program loops forever on it - and an empty list too - the program crashes) and prints first the longest consistent set of matches (ordered along both strings; overlaps are cut off: the overlap column shows the cut) after '> fileName [reverse] Consistent matches', then all the other matches after '> fileName [reverse] Other matches', in the format of mgaps: start1 start2 length overlap gap1 gap2 (the first line of a block has none - -; the other matches show a gap only after a match they chain to). Ties in the second start keep their input order. reverseComplement (-r) writes 'reverse' into the header lines, which mummer-annotate reads. Byte for byte equal to the installed gaps on 259 random match lists (+ the empty-list crash).

## Core Documentation Reference

- Source: [MummerGaps.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MummerGaps.cs)
- Algorithm doc: [Mummer_Package_Programs.md](../../../algorithms/Alignment/Mummer_Package_Programs.md)
- Oracle: [generate_mummer_gaps_oracle.py](../../../../scripts/alignment/generate_mummer_gaps_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `matches` | string | Yes | The matches: white-space separated numbers, three per match. |
| `fileName` | string | No | The first argument of the program, printed in the header lines (the file name of the first sequence). (default: reference.fa) |
| `reverseComplement` | boolean | No | The second argument -r: the header lines read 'reverse' (default false). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | What gaps writes to standard output. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Matches cannot be null or empty |
| 1001 | gaps: the input contains '...': the program loops forever on text that is not a number. |
| 1001 | gaps: no matches (the program dereferences a null array and crashes). |

## Examples

### Example 1: Consistent and other matches

**User Prompt:**
> Order these matches.

**Expected Tool Call:**
```json
{"tool": "mummer_gaps", "arguments": {"matches": "       3         9        20\n      40        50        15\n      70        80        30\n      30        20         5\n", "fileName": "g1.fa", "reverseComplement": true}}
```

**Response:**
```json
{"text": "> g1.fa reverse Consistent matches\n       3        9     20    none      -      -\n      40       50     15    none     17     21\n      70       80     30    none     15     15\n> g1.fa reverse Other matches\n      30       20      5    none      -      -\n"}
```

gaps g1.fa -r: three matches chain along both strings (gaps 17 / 21 and 15 / 15), the fourth is a match on its own.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
