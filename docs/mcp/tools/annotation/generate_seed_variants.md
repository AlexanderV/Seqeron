# generate_seed_variants

Enumerate single-nucleotide variants of a seed sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `generate_seed_variants` |
| **Method ID** | `MiRnaAnalyzer.GenerateSeedVariants` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the original seed followed by every single-position substitution over the RNA alphabet A/C/G/U.
For a seed of length L this yields `1 + 3L` sequences (the original plus 3 substitutions at each of L
positions). Useful for enumerating near-seed target patterns. The `includeWobble` flag is reserved for a
future wobble-aware expansion and currently does not change the output.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L2715](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L2715)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `seedSequence` | string | Yes | Seed nucleotide sequence |
| `includeWobble` | boolean | No | Reserved (currently unused); default true |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `variants` | string[] | Original seed followed by all single-substitution variants |

## Errors

| Code | Message |
|------|---------|
| 1001 | Seed sequence cannot be null or empty |

## Examples

### Example 1: Two-nucleotide seed

**Response:**
```json
{ "variants": ["AC", "CC", "GC", "UC", "AA", "AG", "AU"] }
```

The original `AC`, then position-0 substitutions (`CC`, `GC`, `UC`) and position-1 substitutions
(`AA`, `AG`, `AU`).

### Example 2: Seven-nucleotide seed

A length-7 seed such as `GAGGUAG` yields `1 + 7·3 = 22` variants, the first of which is the original.

## Performance

- **Time Complexity:** O(L) variants for seed length L
- **Space Complexity:** O(L)

## See Also

- [mirna_seed_sequence](mirna_seed_sequence.md) — extract a seed
- [find_similar_mirnas](find_similar_mirnas.md) — near-seed database search
