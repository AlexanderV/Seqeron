# find_inverted_repeats

Inverted repeats / hairpin candidates in a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_inverted_repeats` |
| **Method ID** | `RepeatFinder.FindInvertedRepeats` |
| **Version** | 1.1.0 |
| **Stability** | Stable |

## Description

Finds **inverted repeats** (exact stems by default): pairs of arms where the right arm equals the
reverse complement of the left arm, separated by a loop of length
`minLoopLength..maxLoopLength`. Arms are at least `minArmLength` long. Only **maximal** stems are
reported: a stem lying inside another stem in both arms (its sub-stems and slipped re-pairings) is
dropped, the reporting rule of EMBOSS `palindrome` (`-nummismatches 0 -overlap Y`); a stem is
extended inward only while the loop stays ≥ `minLoopLength`. Only A/C/G/T pair (N never pairs).
Optional: `maxMismatches` (EMBOSS `palindrome -nummismatches`: stems start at a pairing outer pair, extend inward
until the (k+1)-th mismatch, inner-end mismatches trimmed; the count is returned in `mismatches`), `maxArmLength`
(`palindrome -maxpallen`: longer stems are neither reported nor split) and `allowWobble` (G·U / G·T pairs).
No gaps (bulges) — scored gapped inverted repeats are `RepeatFinder.FindInvertedRepeatsScored` (EMBOSS `einverted`, C# API).
Such structures can fold into hairpins; `canFormHairpin` is true when the loop is ≥ 3 nt.
`totalLength = 2 × armLength + loopLength`. Items are ordered by `leftArmStart`, then `rightArmStart`.

## Core Documentation Reference

- Source: [RepeatFinder.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs) (`FindInvertedRepeats`)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `minArmLength` | integer | No | Minimum arm length (default 4, ≥ 2) |
| `maxLoopLength` | integer | No | Maximum loop length (default 50, ≥ `minLoopLength`) |
| `minLoopLength` | integer | No | Minimum loop length (default 3, ≥ 0) |
| `maxMismatches` | integer | No | Maximum mismatched pairs inside a stem (default 0, ≥ 0) |
| `maxArmLength` | integer | No | Maximum arm length (default 2147483647 = unbounded, ≥ `minArmLength`) |
| `allowWobble` | boolean | No | Allow G·U (G·T) wobble pairs (default false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | `{ leftArmStart, rightArmStart, armLength, loopLength, leftArm, rightArm, loop, canFormHairpin, totalLength, mismatches }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | minArmLength must be ≥ 2; minLoopLength must be ≥ 0; maxLoopLength must be ≥ minLoopLength; maxMismatches must be ≥ 0; maxArmLength must be ≥ minArmLength |

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
{ "items": [ { "leftArmStart": 0, "rightArmStart": 7, "armLength": 4, "loopLength": 3, "leftArm": "GGGG", "rightArm": "CCCC", "loop": "AAA", "canFormHairpin": true, "totalLength": 11, "mismatches": 0 } ] }
```
CCCC is the reverse complement of GGGG, with a 3-nt AAA loop.

### Example 2: Only the maximal stem is reported

**User Prompt:**
> Find inverted repeats in "GGGGGAAACCCCC".

**Expected Tool Call:**
```json
{
  "tool": "find_inverted_repeats",
  "arguments": { "sequence": "GGGGGAAACCCCC", "minArmLength": 4, "maxLoopLength": 50, "minLoopLength": 3 }
}
```

**Response:**
```json
{ "items": [ { "leftArmStart": 0, "rightArmStart": 8, "armLength": 5, "loopLength": 3, "leftArm": "GGGGG", "rightArm": "CCCCC", "loop": "AAA", "canFormHairpin": true, "totalLength": 13, "mismatches": 0 } ] }
```
The arm-4 sub-stems GGGG/CCCC (1, 8) and (0, 9) lie inside this stem and are not reported
(EMBOSS palindrome gives the same single stem).

### Example 3: No inverted repeat

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

### Example 4: One mismatch allowed

**User Prompt:**
> Find inverted repeats in "GAATTCAGGAAAACCTCAATTC" allowing one mismatch, arms ≥ 4, loop 0–10.

**Expected Tool Call:**
```json
{
  "tool": "find_inverted_repeats",
  "arguments": { "sequence": "GAATTCAGGAAAACCTCAATTC", "minArmLength": 4, "maxLoopLength": 10, "minLoopLength": 0, "maxMismatches": 1 }
}
```

**Response:**
```json
{ "items": [ { "leftArmStart": 0, "rightArmStart": 13, "armLength": 9, "loopLength": 4, "leftArm": "GAATTCAGG", "rightArm": "CCTCAATTC", "loop": "AAAA", "canFormHairpin": true, "totalLength": 22, "mismatches": 1 } ] }
```
Same stem as EMBOSS `palindrome -nummismatches 1` (gaattcagg / cttaactcc, one interior G·C→G·T mismatch); with
`maxMismatches` 0 nothing is found.

## Performance

- **Time Complexity:** O(n · W) candidate scan (W = maxLoopLength − minLoopLength + 1, capped by n) plus
  stem extension and an O(W · arm) containment check per candidate stem.
- **Space Complexity:** O(number of inverted repeats).
- With `maxMismatches > 0` or a bounded arm: O(n · (W/2 + k + arm)) diagonal walk plus a binary-search nesting test.

## See Also

- [find_palindromes](find_palindromes.md) — zero-loop inverted repeats
- [find_rna_inverted_repeats](find_rna_inverted_repeats.md) — RNA hairpin stems
- [find_stem_loops](find_stem_loops.md) — RNA stem-loops with energy
