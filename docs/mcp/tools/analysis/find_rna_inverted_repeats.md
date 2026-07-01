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

Finds **antiparallel complementary regions** — a left arm whose reverse complement
matches a right arm across a loop of `minSpacing..maxSpacing` nt. These are potential
RNA hairpin stems. Each result reports the two arm spans and the arm length. This is
distinct from the DNA `find_inverted_repeats` (RNA complement rules, arm/loop semantics).

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs#L2193](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs#L2193)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | RNA sequence (min length 1) |
| `minLength` | integer | No | Minimum arm length (default 4) |
| `minSpacing` | integer | No | Minimum spacing between arms (default 3) |
| `maxSpacing` | integer | No | Maximum spacing between arms (default 100) |

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
