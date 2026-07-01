# mirna_seed_sequence

Extract the canonical seed region (positions 2-8) from a miRNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `mirna_seed_sequence` |
| **Method ID** | `MiRnaAnalyzer.GetSeedSequence` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the miRNA seed region — nucleotides 2–8 (0-based `Substring(1, 7)`), upper-cased. The seed is the
primary determinant of miRNA target specificity (Bartel 2009). The sequence must be at least 8 nt long.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L96](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L96)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `miRnaSequence` | string | Yes | miRNA nucleotide sequence (≥ 8 nt) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `seed` | string | Seed region (positions 2–8), upper-cased |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | miRNA sequence must be at least 8 nt long to extract seed |

## Examples

### Example 1: let-7a seed

**User Prompt:**
> What's the seed of let-7a (UGAGGUAGUAGGUUGUAUAGUU)?

**Response:**
```json
{ "seed": "GAGGUAG" }
```

### Example 2: Lower-case input

```json
{ "tool": "mirna_seed_sequence", "arguments": { "miRnaSequence": "ugagguaguag" } }
```

**Response:**
```json
{ "seed": "GAGGUAG" }
```

## Performance

- **Time Complexity:** O(1)
- **Space Complexity:** O(1)

## See Also

- [create_mirna](create_mirna.md) — build a full miRNA record with seed metadata
- [compare_seed_regions](compare_seed_regions.md) — seed-family comparison
