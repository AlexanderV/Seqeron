# simulate_bisulfite_conversion

Simulate in-silico sodium-bisulfite conversion of a DNA strand.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `simulate_bisulfite_conversion` |
| **Method ID** | `EpigeneticsAnalyzer.SimulateBisulfiteConversion` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Models bisulfite treatment (Frommer et al. 1992): unmethylated cytosines are deaminated to uracil, which
reads and amplifies as thymine, while 5-methylcytosines are non-reactive and remain cytosines. Each
unprotected `C`/`c` is replaced with `T`/`t` (case preserved); cytosines whose 0-based index is listed in
`methylatedPositions` stay cytosines; all non-cytosine bases pass through unchanged. Only the supplied
strand is converted (the protocol is strand-specific). The result has the same length as the input.

## Core Documentation Reference

- Source: [EpigeneticsAnalyzer.cs#L382](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs#L382)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence (min length: 1) |
| `methylatedPositions` | integer[] \| null | No | 0-based indices of protected cytosines; null = fully unmethylated |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `converted` | string | Bisulfite-converted strand (same length as input) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Fully unmethylated

**User Prompt:**
> Bisulfite-convert ACGCGT with no methylation.

**Response:**
```json
{ "converted": "ATGTGT" }
```

Both cytosines (indices 1 and 3) convert to thymine.

### Example 2: Protecting a methylated cytosine

```json
{ "tool": "simulate_bisulfite_conversion", "arguments": { "sequence": "ACGCGT", "methylatedPositions": [1] } }
```

**Response:**
```json
{ "converted": "ACGTGT" }
```

The protected C at index 1 survives; the C at index 3 still converts.

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(n)

## See Also

- [methylation_from_bisulfite](methylation_from_bisulfite.md) — call methylation levels from bisulfite reads
- [find_methylation_sites](find_methylation_sites.md) — candidate methylation contexts
