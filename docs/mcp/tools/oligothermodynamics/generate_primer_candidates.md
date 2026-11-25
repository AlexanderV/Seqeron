# generate_primer_candidates

Enumerate and evaluate every primer candidate in a region.

## Overview

| Property | Value |
|----------|-------|
| **Server** | OligoThermodynamics |
| **Tool Name** | `generate_primer_candidates` |
| **Method ID** | `ThermodynamicPrimerDesigner.GeneratePrimerCandidates` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

For every start position `s` in `[region_start, region_end)` and every admissible length `len ∈ [parameters.MinLength, parameters.MaxLength]` with `s + len ≤ region_end`, emits an evaluated `PrimerCandidate`. Candidates are returned in **generation order** (increasing start, then increasing length), **not** sorted by score, and both valid and invalid candidates are included. When `forward = false`, each candidate sequence is the reverse complement of the template substring at `[s, s+len)`.

## Core Documentation Reference

- Source: [Seqeron.Genomics.OligoThermodynamics/ThermodynamicPrimerDesigner.Core.cs#L2177](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.OligoThermodynamics/ThermodynamicPrimerDesigner.Core.cs#L2177)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `template` | string | Yes | DNA template sequence. |
| `region_start` | integer | Yes | 0-based start of the search region (inclusive). |
| `region_end` | integer | Yes | 0-based end of the search region (exclusive). |
| `forward` | boolean | No | True for forward-strand candidates, false for reverse-complement candidates. (default: true) |
| `parameters` | object | No | Optional primer design parameters. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `candidates` | array | Evaluated `PrimerCandidate` records in generation order. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Template cannot be null or empty |
| 1002 | Region start must be non-negative |
| 1003 | Region end must not exceed the template length |
| 1004 | Region start must be strictly less than region end |

## Examples

### Example 1: Forward candidates over [0, 60)

With default parameters (lengths 18–25) the standard 258 bp template yields **316** candidates; the first is start 0, length 18.

**Input:** `{ "region_start": 0, "region_end": 60, "forward": true }`

**Response (abridged):** `{ "candidates": [ { "position": 0, "length": 18, ... }, ... ] }` (316 items)

### Example 2: Reverse candidates are reverse complements

For template `AACCGGTT…` (40 bp), `forward = false`, the first candidate at position 0 length 18 has sequence `TTAACCGGTTAACCGGTT` = reverse complement of `AACCGGTTAACCGGTTAA`.

## See Also

- [design_primers](design_primers.md), [evaluate_primer](evaluate_primer.md)
