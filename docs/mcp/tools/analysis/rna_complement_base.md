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

Returns the **RNA complement** of a single base, matching Biopython 1.88 `complement_rna`:
A→U, U→A, G↔C, and DNA T→A (T is read as U's DNA equivalent). IUPAC ambiguity codes are
complemented in the RNA alphabet (IUPAC NC-IUB 1984): R↔Y, K↔M, B↔V, D↔H; S, W and N map to
themselves. Input is case-insensitive; recognized codes are returned upper-case; any other character
(e.g. a gap) passes through unchanged. Core method: `SequenceExtensions.GetRnaComplementBase`
(via `RnaSecondaryStructure.GetComplement`).

Cross-check: `complement_rna("ACGTURYSWKMBDHVN")` = `UGCAAYRSWMKVHDBN`.

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
