# find_stem_loops

Enumerate hairpin stem-loop candidates in an RNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_stem_loops` |
| **Method ID** | `RnaSecondaryStructure.FindStemLoops` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Enumerates **hairpin stem-loop candidates** in an RNA sequence. Each result reports
the stem (5'/3' bounds, length, base pairs, energy), the loop (type, bounds, size,
sequence), the dot-bracket notation, and the Turner 2004 total free energy. Stems must
be at least `minStemLength` base pairs and loops within `[minLoopSize, maxLoopSize]`.
G-U wobble pairs are allowed in the stem when `allowWobble` is true.

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs#L458](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs#L458)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `rnaSequence` | string | Yes | RNA sequence (min length 1) |
| `minStemLength` | integer | No | Minimum stem length (default 3) |
| `minLoopSize` | integer | No | Minimum loop size (default 3) |
| `maxLoopSize` | integer | No | Maximum loop size (default 10) |
| `allowWobble` | boolean | No | Allow G-U wobble pairs in the stem (default true) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | Stem-loops: `{ start, end, stem, loop, totalFreeEnergy, dotBracketNotation }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Simple 3bp hairpin

**User Prompt:**
> Find stem-loops in "GGGAAAACCC" (3bp stem, 4nt loop).

**Expected Tool Call:**
```json
{
  "tool": "find_stem_loops",
  "arguments": { "rnaSequence": "GGGAAAACCC", "minStemLength": 3, "minLoopSize": 4, "maxLoopSize": 4 }
}
```

**Response:**
```json
{ "items": [ { "start": 0, "end": 9, "dotBracketNotation": "(((....)))" } ] }
```
GGG pairs with CCC around a 4-nt AAAA loop.

### Example 2: No hairpin

**User Prompt:**
> Find stem-loops in "AAAAAAAA".

**Expected Tool Call:**
```json
{
  "tool": "find_stem_loops",
  "arguments": { "rnaSequence": "AAAAAAAA" }
}
```

**Response:**
```json
{ "items": [] }
```

## Performance

- **Time Complexity:** O(n²) over start positions and stem lengths.
- **Space Complexity:** O(number of stem-loops).

## See Also

- [predict_rna_structure](predict_rna_structure.md) — full structure prediction
- [minimum_free_energy](minimum_free_energy.md) — Zuker MFE
- [find_rna_inverted_repeats](find_rna_inverted_repeats.md) — complementary arms
