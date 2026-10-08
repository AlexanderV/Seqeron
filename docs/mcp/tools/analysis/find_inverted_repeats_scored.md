# find_inverted_repeats_scored

Find gapped, mismatch-tolerant inverted repeats scored like EMBOSS einverted.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_inverted_repeats_scored` |
| **Method ID** | `RepeatFinder.FindInvertedRepeatsScored` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds inverted repeats with EMBOSS `einverted` (Durbin & Thierry-Mieg 1993; EMBOSS 6.6.0): a local
alignment of the sequence against its reverse complement in which a Watson–Crick pair scores `+matchScore`, any
other pair `mismatchScore` and each gap `−gapPenalty`. Every repeat whose score reaches `threshold` is reported in
einverted's report order with both alignment rows and the match line. Coordinates are **0-based and inclusive**
(einverted prints 1-based); alignment rows use the upper-cased input (einverted prints lower case). Symbols other
than A/C/G/T never pair. `maxRepeatLength` bounds the extent from the start of the repeat to the end of its
inverted copy (einverted `-maxrepeat`; memory O(min(maxRepeatLength, n)²)). Cross-checked against the einverted
binary (Evidence REP-INV-001). Use `find_inverted_repeats` for exact/maximal stems (EMBOSS palindrome).

## Core Documentation Reference

- Source: [RepeatFinder.cs#L4066](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs#L4066)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (case-insensitive; symbols other than A/C/G/T never pair) |
| `gapPenalty` | integer | No | Gap penalty (einverted -gap) (default 12) |
| `threshold` | integer | No | Minimum reported score (einverted -threshold) (default 50) |
| `matchScore` | integer | No | Score of a Watson-Crick pair (einverted -match) (default 3) |
| `mismatchScore` | integer | No | Score of any other pair (einverted -mismatch) (default -4) |
| `maxRepeatLength` | integer | No | Maximum extent of the repeat incl. its inverted copy (einverted -maxrepeat) (default 2000) |

## Output Schema

`items`: `leftArmStart, leftArmEnd, rightArmStart, rightArmEnd` (0-based inclusive), `score, matches, mismatches, gaps`, `leftArmAlignment, matchLine, rightArmAlignment`, `leftArmLength, rightArmLength, loopLength, percentMatches`

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | gapPenalty, threshold, matchScore must be >= 0; mismatchScore <= 0; maxRepeatLength >= 2 |

## Examples

### Example 1: einverted defaults on a planted imperfect repeat ("Score 60: 28/31 (90%) matches, 1 gaps")

**Input:** `{"sequence": "CCCAACCCATGCGTACGTTAGCCTAGGATCCATTTTTTTTTGGATACTAGGCAACGTACGCATGGGAAGGG"}`

**Output:**

```json
{"items": [{"leftArmStart": 0, "leftArmEnd": 31, "rightArmStart": 40, "rightArmEnd": 70, "score": 60, "matches": 28, "mismatches": 3, "gaps": 1, "leftArmAlignment": "CCCAACCCATGCGTACGTTAGCCTAGGATCCA", "matchLine": "|||  |||||||||||||| |||||| |||||", "rightArmAlignment": "GGGAAGGGTACGCATGCAA-CGGATCATAGGT", "leftArmLength": 32, "rightArmLength": 31, "loopLength": 8, "percentMatches": 90.3225806451613}]}
```

### Example 2: einverted -gap 8 -threshold 30 (Score 64 = 84 - 12 - 8)

**Input:** `{"sequence": "CCCAACCCATGCGTACGTTAGCCTAGGATCCATTTTTTTTTGGATACTAGGCAACGTACGCATGGGAAGGG", "gapPenalty": 8, "threshold": 30}`

**Output:**

```json
{"items": [{"leftArmStart": 0, "leftArmEnd": 31, "rightArmStart": 40, "rightArmEnd": 70, "score": 64, "matches": 28, "mismatches": 3, "gaps": 1, "leftArmAlignment": "CCCAACCCATGCGTACGTTAGCCTAGGATCCA", "matchLine": "|||  |||||||||||||| |||||| |||||", "rightArmAlignment": "GGGAAGGGTACGCATGCAA-CGGATCATAGGT", "leftArmLength": 32, "rightArmLength": 31, "loopLength": 8, "percentMatches": 90.3225806451613}]}
```

### Example 3: No repeat above threshold

**Input:** `{"sequence": "ACGTACGTAC"}`

**Output:**

```json
{"items": []}
```

## Performance

- O(n · min(maxRepeatLength, n)) time; O(min(maxRepeatLength, n)²) memory.

## See Also

- [find_inverted_repeats](find_inverted_repeats.md)
- [find_reverse_complement_repeats](find_reverse_complement_repeats.md)
- [find_degenerate_repeats](find_degenerate_repeats.md)
