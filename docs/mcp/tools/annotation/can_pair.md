# can_pair

Test whether two RNA bases can pair.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `can_pair` |
| **Method ID** | `MiRnaAnalyzer.CanPair` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns whether two nucleotides can form a base pair: Watson-Crick `A-U` and `G-C`, or the `G-U` wobble
pair. Input is case-insensitive and DNA `T` is treated as RNA `U`.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L326](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L326)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `base1` | string | Yes | First base (single character) |
| `base2` | string | Yes | Second base (single character) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `canPair` | boolean | true for A-U, G-C, or G-U |

## Errors

| Code | Message |
|------|---------|
| 1001 | base1 must be a single character |
| 1001 | base2 must be a single character |

## Examples

### Example 1: Watson-Crick pair

```json
{ "tool": "can_pair", "arguments": { "base1": "A", "base2": "U" } }
```

**Response:**
```json
{ "canPair": true }
```

### Example 2: Non-pairing bases

```json
{ "canPair": false }
```

(for `A` and `G`)

## Performance

- **Time Complexity:** O(1)
- **Space Complexity:** O(1)

## See Also

- [is_wobble_pair](is_wobble_pair.md) — G-U wobble test only
- [rna_reverse_complement](rna_reverse_complement.md) — reverse complement
