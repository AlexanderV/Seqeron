# validate_dot_bracket

Validate a dot-bracket secondary-structure string.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `validate_dot_bracket` |
| **Method ID** | `RnaSecondaryStructure.ValidateDotBracket` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Validates that every bracket symbol in a **dot-bracket** (or extended WUSS) string is
balanced and family-matched: each closing symbol matches an earlier unmatched opening
symbol of the same family (`()`, `[]`, `{}`, `<>`, letter pairs), and no opener is left
unclosed. A mispairing such as `(]` is rejected even though the total count is balanced.

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs#L2114](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs#L2114)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `dotBracket` | string | Yes | Dot-bracket secondary-structure string (min length 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `result` | boolean | True when the structure is balanced |

## Errors

| Code | Message |
|------|---------|
| 1001 | Dot-bracket string cannot be null or empty |

## Examples

### Example 1: Balanced

**User Prompt:**
> Is "(())" a valid dot-bracket structure?

**Expected Tool Call:**
```json
{
  "tool": "validate_dot_bracket",
  "arguments": { "dotBracket": "(())" }
}
```

**Response:**
```json
{ "result": true }
```

### Example 2: Unbalanced

**User Prompt:**
> Is "(()" valid?

**Expected Tool Call:**
```json
{
  "tool": "validate_dot_bracket",
  "arguments": { "dotBracket": "(()" }
}
```

**Response:**
```json
{ "result": false }
```

## Performance

- **Time Complexity:** O(n).
- **Space Complexity:** O(n) for the stacks.

## See Also

- [parse_dot_bracket](parse_dot_bracket.md) — extract base pairs
- [predict_rna_structure](predict_rna_structure.md) — produce a dot-bracket structure
