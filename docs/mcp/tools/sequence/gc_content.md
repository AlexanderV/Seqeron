# gc_content

Calculate GC content of a DNA/RNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `gc_content` |
| **Method ID** | `SequenceExtensions.CalculateGcContentFast` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates the GC content (percentage of G and C nucleotides) of a DNA/RNA sequence. GC content is an important metric for primer design, genome analysis, and sequence characterization.

Formula: GC% = (G + C) / (A + C + G + T + U) × 100, case-insensitive. IUPAC ambiguity codes (including S and W), N, gaps and every other character are excluded from both the numerator and the denominator; a sequence with no A/C/G/T/U returns 0. The percentage and both counts come from the Core canonical primitives (`SequenceExtensions.CalculateGcContentFast` and `SequenceExtensions.CountGcAndValidNucleotides`), so the tool output is identical to the C# API.

Cross-check (Biopython 1.88 `gc_fraction(seq, "remove")` × 100): `"ATGCGATCGATCG"` → 53.846…, `"ACGT-NNNN"` → 50, `"acguACGU"` → 50. Difference: Biopython counts S as G/C and W as A/T; this tool excludes them like every other ambiguity code (`"ACGTSS"`: Biopython 66.67, tool 50 — use `SequenceExtensions.CalculateGcFraction(seq, GcAmbiguityMode.Remove)` in C# for exact Biopython parity).

## Core Documentation Reference

- Source: [SequenceExtensions.cs#L78](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs#L78)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA or RNA sequence (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `gcContent` | number | GC content percentage (0-100) |
| `gcCount` | integer | Number of G and C nucleotides |
| `totalCount` | integer | Number of valid A/C/G/T/U nucleotides (the denominator); IUPAC codes, N and gaps are not counted, so this is not the sequence length when such characters are present |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: 50% GC content

**User Prompt:**
> What's the GC content of "ATGC"?

**Expected Tool Call:**
```json
{
  "tool": "gc_content",
  "arguments": {
    "sequence": "ATGC"
  }
}
```

**Response:**
```json
{
  "gcContent": 50,
  "gcCount": 2,
  "totalCount": 4
}
```

### Example 2: GC-rich sequence

**User Prompt:**
> Calculate GC content of "GCGCGC"

**Expected Tool Call:**
```json
{
  "tool": "gc_content",
  "arguments": {
    "sequence": "GCGCGC"
  }
}
```

**Response:**
```json
{
  "gcContent": 100,
  "gcCount": 6,
  "totalCount": 6
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1)

## See Also

- [nucleotide_composition](nucleotide_composition.md) - Full composition analysis
- [summarize_sequence](summarize_sequence.md) - Comprehensive statistics
