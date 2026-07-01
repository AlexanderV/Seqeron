# find_degenerate_motif

Find IUPAC-degenerate motif matches in a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_degenerate_motif` |
| **Method ID** | `MotifFinder.FindDegenerateMotif` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds all occurrences of a motif expressed with IUPAC ambiguity codes (N, R, Y, S, W, K, M, B,
D, H, V). Each position is matched against the code's base set (e.g. R = A/G, N = any). Every
match reports its 0-based position, the matched substring, the uppercased pattern, and a score
of 1.0. The pattern is validated: an unknown IUPAC code raises an error. Input must be valid DNA.

## Core Documentation Reference

- Source: [MotifFinder.cs#L90](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs#L90)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `motif` | string | Yes | IUPAC motif pattern |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | Matches: `position, matchedSequence, pattern, score` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |
| 1002 | Invalid IUPAC code … in motif pattern |

## Examples

### Example 1: Purine-led motif

**Input:** `{ "sequence": "AGGTAG", "motif": "RGG" }`
→ R = A/G; only position 0 "AGG" matches →
`{ "items": [ { "position": 0, "matchedSequence": "AGG", "pattern": "RGG", "score": 1.0 } ] }`

### Example 2: Wildcard N

**Input:** `{ "sequence": "ACG", "motif": "N" }`
→ matches every position → items at 0, 1, 2.

## Performance

- **Time Complexity:** O(n · m) for pattern length m. **Space Complexity:** O(matches).

## See Also

- [find_exact_motif](find_exact_motif.md)
- [find_regulatory_elements](find_regulatory_elements.md)
