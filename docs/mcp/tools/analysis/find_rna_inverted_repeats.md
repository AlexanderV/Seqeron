# find_rna_inverted_repeats

Antiparallel complementary regions (RNA hairpin stems).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_rna_inverted_repeats` |
| **Method ID** | `RnaSecondaryStructure.FindInvertedRepeats` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds the **maximal antiparallel stems** — a left arm that pairs with a right arm across a loop of
`minSpacing..maxSpacing` nt. These are potential RNA hairpin stems. Pairs are A·U (A·T), G·C and, by default,
G·U wobble (ViennaRNA `CanPair`; `T` is read as `U`, ambiguity codes never pair). A stem is dropped only when it
lies inside another stem in both arms; alternative stems that share bases are all reported. Coordinates are
0-based inclusive; results are ordered by `start1`, then `start2`. The search is the same engine as the DNA
`find_inverted_repeats` (EMBOSS `palindrome` stems); only the pairing rule differs. Gapped/scored (einverted)
repeats: `find_inverted_repeats_scored`.

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs) (`FindInvertedRepeats`), delegating to `RepeatFinder` ([Inverted_Repeats.md](../../../algorithms/RnaStructure/Inverted_Repeats.md))

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | RNA sequence (min length 1) |
| `minLength` | integer | No | Minimum arm length (default 4) |
| `minSpacing` | integer | No | Minimum spacing between arms (default 3) |
| `maxSpacing` | integer | No | Maximum spacing between arms (default 100) |
| `maxMismatches` | integer | No | Maximum interior mismatches per stem (default 0; IUPACpal k / palindrome `-nummismatches`) |
| `allowWobble` | boolean | No | Count G·U as a pair (default true); false = A·U / G·C only |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | `{ start1, end1, start2, end2, length }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: GCGC hairpin stem

**User Prompt:**
> Find RNA inverted repeats in "GCGCAAAAAAGCGC".

**Expected Tool Call:**
```json
{ "tool": "find_rna_inverted_repeats", "arguments": { "sequence": "GCGCAAAAAAGCGC", "minLength": 4, "minSpacing": 3, "maxSpacing": 100 } }
```

**Response:**
```json
{ "items": [ { "start1": 0, "end1": 3, "start2": 10, "end2": 13, "length": 4 } ] }
```
The 5' GCGC pairs antiparallel with the 3' GCGC across a 6-nt loop.

### Example 2: No complementary arms

**User Prompt:**
> RNA inverted repeats in a poly-A tract?

**Expected Tool Call:**
```json
{ "tool": "find_rna_inverted_repeats", "arguments": { "sequence": "AAAAAAAAAA" } }
```

**Response:**
```json
{ "items": [] }
```

## Performance

- **Time Complexity:** O(n² · maxSpacing) worst case.
- **Space Complexity:** O(number of repeats).

## See Also

- [find_stem_loops](find_stem_loops.md)
- [find_inverted_repeats](find_inverted_repeats.md) — DNA variant
