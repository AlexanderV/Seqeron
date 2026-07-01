# parse_dot_bracket

Parse dot-bracket notation into base-pair coordinates.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `parse_dot_bracket` |
| **Method ID** | `RnaSecondaryStructure.ParseDotBracket` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Parses **dot-bracket** notation into a list of base-pair coordinate tuples using one
stack per bracket family (`()`, `[]`, `{}`, `<>`, and upper/lower letter pairs).
Unpaired symbols (`.`, `,`, `-`, `:`, `_`, `~`) are skipped. Each pair is
`(openingIndex, closingIndex)`.

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs#L2068](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs#L2068)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `dotBracket` | string | Yes | Dot-bracket secondary-structure string (min length 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `pairs` | array | Base-pair coordinates `{ position1, position2 }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Dot-bracket string cannot be null or empty |

## Examples

### Example 1: Nested pair

**User Prompt:**
> Parse the base pairs of "(())".

**Expected Tool Call:**
```json
{
  "tool": "parse_dot_bracket",
  "arguments": { "dotBracket": "(())" }
}
```

**Response:**
```json
{ "pairs": [ { "position1": 1, "position2": 2 }, { "position1": 0, "position2": 3 } ] }
```
The inner pair (1,2) is closed first, then the outer pair (0,3).

### Example 2: Fully unpaired

**User Prompt:**
> Parse "....".

**Expected Tool Call:**
```json
{
  "tool": "parse_dot_bracket",
  "arguments": { "dotBracket": "...." }
}
```

**Response:**
```json
{ "pairs": [] }
```

## Performance

- **Time Complexity:** O(n).
- **Space Complexity:** O(n).

## See Also

- [validate_dot_bracket](validate_dot_bracket.md) — balance validation
- [detect_pseudoknots](detect_pseudoknots.md) — find crossing pairs
