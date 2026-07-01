# find_inverted_repeats

Inverted repeats / hairpin candidates in a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_inverted_repeats` |
| **Method ID** | `RepeatFinder.FindInvertedRepeats` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds **inverted repeats**: pairs of arms where the right arm equals the reverse
complement of the left arm, separated by a loop of length
`minLoopLength..maxLoopLength`. Arms are at least `minArmLength` long. Such structures
can fold into hairpins; `canFormHairpin` is true when the loop is ≥ 3 nt.
`totalLength = 2 × armLength + loopLength`.

## Core Documentation Reference

- Source: [RepeatFinder.cs#L702](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs#L702)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `minArmLength` | integer | No | Minimum arm length (default 4, ≥ 2) |
| `maxLoopLength` | integer | No | Maximum loop length (default 50) |
| `minLoopLength` | integer | No | Minimum loop length (default 3, ≥ 0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | `{ leftArmStart, rightArmStart, armLength, loopLength, leftArm, rightArm, loop, canFormHairpin, totalLength }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | minArmLength must be ≥ 2 |

## Examples

### Example 1: GGGG-AAA-CCCC hairpin

**User Prompt:**
> Find inverted repeats in "GGGGAAACCCC".

**Expected Tool Call:**
```json
{
  "tool": "find_inverted_repeats",
  "arguments": { "sequence": "GGGGAAACCCC", "minArmLength": 4, "maxLoopLength": 50, "minLoopLength": 3 }
}
```

**Response:**
```json
{ "items": [ { "leftArmStart": 0, "rightArmStart": 7, "armLength": 4, "loopLength": 3, "leftArm": "GGGG", "rightArm": "CCCC", "loop": "AAA", "canFormHairpin": true, "totalLength": 11 } ] }
```
CCCC is the reverse complement of GGGG, with a 3-nt AAA loop.

### Example 2: No inverted repeat

**User Prompt:**
> Inverted repeats in "AAAAAAAA"?

**Expected Tool Call:**
```json
{
  "tool": "find_inverted_repeats",
  "arguments": { "sequence": "AAAAAAAA", "minArmLength": 4, "maxLoopLength": 50, "minLoopLength": 3 }
}
```

**Response:**
```json
{ "items": [] }
```
A has complement T, so a poly-A tract has no reverse-complement arms.

## Performance

- **Time Complexity:** O(n² · maxLoopLength) worst case.
- **Space Complexity:** O(number of inverted repeats).

## See Also

- [find_palindromes](find_palindromes.md) — zero-loop inverted repeats
- [find_rna_inverted_repeats](find_rna_inverted_repeats.md) — RNA hairpin stems
- [find_stem_loops](find_stem_loops.md) — RNA stem-loops with energy
