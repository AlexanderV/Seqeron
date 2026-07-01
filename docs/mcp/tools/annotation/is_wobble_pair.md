# is_wobble_pair

Test whether two bases form a G-U wobble pair.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `is_wobble_pair` |
| **Method ID** | `MiRnaAnalyzer.IsWobblePair` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns whether two nucleotides form a `G-U` wobble pair (either order). Input is case-insensitive and DNA
`T` is treated as RNA `U`. Watson-Crick pairs return `false` — see [can_pair](can_pair.md) for any valid
pairing.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L340](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L340)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `base1` | string | Yes | First base (single character) |
| `base2` | string | Yes | Second base (single character) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `isWobble` | boolean | true only for G-U or U-G |

## Errors

| Code | Message |
|------|---------|
| 1001 | base1 must be a single character |
| 1001 | base2 must be a single character |

## Examples

### Example 1: G-U wobble

```json
{ "tool": "is_wobble_pair", "arguments": { "base1": "G", "base2": "U" } }
```

**Response:**
```json
{ "isWobble": true }
```

### Example 2: Watson-Crick is not wobble

```json
{ "isWobble": false }
```

(for `G` and `C`)

## Performance

- **Time Complexity:** O(1)
- **Space Complexity:** O(1)

## See Also

- [can_pair](can_pair.md) — any valid base pair including G-U
