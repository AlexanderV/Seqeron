# summarize_sequence

Generate comprehensive summary statistics for a DNA/RNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `summarize_sequence` |
| **Method ID** | `SequenceStatistics.SummarizeNucleotideSequence` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Generates comprehensive summary statistics for a DNA/RNA sequence in a single call. Includes length, nucleotide composition, GC content (fraction 0–1 over A/C/G/T/U), Shannon entropy (bits), linguistic complexity, basic melting temperature (Wallace rule below 14 A/C/G/T/U bases, else 64.9 + 41·(GC − 16.4)/N; U read as T) and the single-stranded linear average molecular weight (Biopython `molecular_weight`; RNA table when the sequence has U and no T, DNA otherwise; unknown symbols skipped). With `weightAmbiguousBasesInTm` = true the Tm weights IUPAC ambiguity codes as Biopython non-strict `Tm_Wallace` / `Tm_GC` (e.g. `acgtNNSWacgt-ry` → 42.0 instead of 24.0).

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L2745](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L2745)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA or RNA sequence to analyze (min length: 1) |
| `weightAmbiguousBasesInTm` | boolean | No | Weight IUPAC codes in the Tm as Biopython `Tm_Wallace`/`Tm_GC(strict=False)` (default `false`: A/C/G/T/U only) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `length` | integer | Sequence length |
| `gcContent` | number | GC content fraction (0–1) |
| `entropy` | number | Shannon entropy (bits) |
| `complexity` | number | Linguistic complexity |
| `meltingTemperature` | number | Melting temperature (°C) |
| `composition` | object | Nucleotide counts (A, T, G, C, U, N) |
| `molecularWeight` | number | Single-stranded linear average molecular weight (Da) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: DNA sequence summary

**User Prompt:**
> Summarize the sequence "ATGCGATCGATCG"

**Expected Tool Call:**
```json
{
  "tool": "summarize_sequence",
  "arguments": {
    "sequence": "ATGCGATCGATCG"
  }
}
```

**Response:**
```json
{
  "length": 13,
  "gcContent": 0.46153846153846156,
  "entropy": 1.9877733714879842,
  "complexity": 0.7037037037037037,
  "meltingTemperature": 40,
  "composition": { "A": 3, "T": 3, "G": 4, "C": 3, "U": 0, "N": 0 },
  "molecularWeight": 4054.5834
}
```

## Performance

- **Time Complexity:** O(n × k) where n is sequence length, k is max k-mer size
- **Space Complexity:** O(4^k) for k-mer storage

## See Also

- [nucleotide_composition](nucleotide_composition.md) - Basic composition only
- [shannon_entropy](shannon_entropy.md) - Entropy calculation
- [linguistic_complexity](linguistic_complexity.md) - Complexity calculation
- [melting_temperature](melting_temperature.md) - Tm calculation
