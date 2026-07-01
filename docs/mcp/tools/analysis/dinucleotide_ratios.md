# dinucleotide_ratios

Observed/expected dinucleotide relative abundance (odds ratios).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `dinucleotide_ratios` |
| **Method ID** | `SequenceStatistics.CalculateDinucleotideRatios` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the dinucleotide relative abundance (odds ratio) `ρ_XY = f_XY / (f_X · f_Y)`, where
`f_XY` is the dinucleotide frequency and `f_X`, `f_Y` are single-base frequencies. `ρ = 1`
means no bias; `ρ > 1` over-representation, `ρ < 1` under-representation (Karlin & Burge 1995).
This is the standard measure behind CpG-island detection. When a constituent base is absent the
expected frequency is 0 and the ratio is reported as 0. Sequences shorter than 2 bp yield an
empty map.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L641](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L641)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence (min length 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `ratios` | object | Map of dinucleotide → observed/expected ratio |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Alternating AT

**Input:** `{ "sequence": "ATAT" }`
→ f_A = f_T = 0.5; f_AT = 2/3, f_TA = 1/3 →
`ρ_AT = (2/3)/0.25 = 2.6667`, `ρ_TA = (1/3)/0.25 = 1.3333`.

**Response:**
```json
{ "ratios": { "AT": 2.6667, "TA": 1.3333 } }
```

### Example 2: Homopolymer

**Input:** `{ "sequence": "AAAA" }`
→ f_A = 1, f_AA = 1 → `ρ_AA = 1 / (1·1) = 1.0`.

## Performance

- **Time Complexity:** O(n). **Space Complexity:** O(distinct dinucleotides).

## See Also

- [dinucleotide_frequencies](dinucleotide_frequencies.md)
