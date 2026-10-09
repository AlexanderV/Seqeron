# complement_base

Get the complement of a single nucleotide base (DNA alphabet by default, RNA alphabet with `rna=true`).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `complement_base` |
| **Method ID** | `SequenceExtensions.GetComplementBase` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the complementary base for a given nucleotide (IUPAC-complete). Input is case-insensitive; recognized codes are returned upper-case; any other character passes through unchanged.

| Mode | Core method | A | T | U | G | C | Reference |
|------|-------------|---|---|---|---|---|-----------|
| `rna=false` (default) | `SequenceExtensions.GetComplementBase` | T | A | A | C | G | Biopython 1.88 `complement` |
| `rna=true` | `SequenceExtensions.GetRnaComplementBase` | U | A | A | C | G | Biopython 1.88 `complement_rna` |

Ambiguity codes are complemented identically in both modes (IUPAC NC-IUB 1984): R↔Y, K↔M, B↔V, D↔H; S, W and N map to themselves.
Cross-check: `complement("ACGTURYSWKMBDHVN")` = `TGCAAYRSWMKVHDBN`, `complement_rna("ACGTURYSWKMBDHVN")` = `UGCAAYRSWMKVHDBN`.
The default mode always emits the DNA alphabet, so `A` → `T` even for an RNA input; pass `rna=true` to get `A` → `U`.

## Core Documentation Reference

- Source: [SequenceExtensions.cs#L83](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs#L83)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `nucleotide` | string | Yes | The nucleotide base (A, C, G, T, U or an IUPAC ambiguity code), exactly 1 character |
| `rna` | boolean | No | Emit the RNA alphabet (A→U, Biopython `complement_rna`) instead of DNA (A→T). Default `false` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `complement` | string | The complementary base |
| `original` | string | The input base |

## Errors

| Code | Message |
|------|---------|
| 1001 | Must provide exactly one nucleotide character |

## Examples

### Example 1: DNA complement

**User Prompt:**
> What's the complement of "A"?

**Expected Tool Call:**
```json
{
  "tool": "complement_base",
  "arguments": {
    "nucleotide": "A"
  }
}
```

**Response:**
```json
{
  "complement": "T",
  "original": "A"
}
```

### Example 2: RNA complement

**User Prompt:**
> Get complement of "U"

**Expected Tool Call:**
```json
{
  "tool": "complement_base",
  "arguments": {
    "nucleotide": "U"
  }
}
```

**Response:**
```json
{
  "complement": "A",
  "original": "U"
}
```

### Example 3: RNA alphabet

**User Prompt:**
> What base pairs with "A" in RNA?

**Expected Tool Call:**
```json
{
  "tool": "complement_base",
  "arguments": {
    "nucleotide": "A",
    "rna": true
  }
}
```

**Response:**
```json
{
  "complement": "U",
  "original": "A"
}
```

## Performance

- **Time Complexity:** O(1)
- **Space Complexity:** O(1)

## See Also

- [dna_reverse_complement](dna_reverse_complement.md) - Full sequence reverse complement
- [rna_complement_base](../analysis/rna_complement_base.md) - RNA complement of a single base (Analysis server)
- [nucleotide_composition](nucleotide_composition.md) - Sequence composition analysis
