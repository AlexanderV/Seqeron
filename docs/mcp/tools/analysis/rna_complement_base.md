# rna_complement_base

RNA complement of a single base.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `rna_complement_base` |
| **Method ID** | `RnaSecondaryStructure.GetComplement` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the **RNA complement** of a single base: A↔U, G↔C (IUPAC ambiguity codes are
also mapped). T is treated as U's DNA equivalent and complements to A.

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs#L449](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs#L449)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `base` | string | Yes | RNA base (length-1 string) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `complement` | string | Complement base |

## Errors

| Code | Message |
|------|---------|
| 1001 | Expected a single character (length-1 string) |

## Examples

### Example 1: A → U

**User Prompt:**
> What is the RNA complement of A?

**Expected Tool Call:**
```json
{
  "tool": "rna_complement_base",
  "arguments": { "base": "A" }
}
```

**Response:**
```json
{ "complement": "U" }
```

### Example 2: G → C

**User Prompt:**
> RNA complement of G.

**Expected Tool Call:**
```json
{
  "tool": "rna_complement_base",
  "arguments": { "base": "G" }
}
```

**Response:**
```json
{ "complement": "C" }
```

## Performance

- **Time Complexity:** O(1).
- **Space Complexity:** O(1).

## See Also

- [can_pair](can_pair.md) — whether two RNA bases pair
- [base_pair_type](base_pair_type.md) — Watson-Crick / wobble classification
