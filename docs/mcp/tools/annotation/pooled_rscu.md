# pooled_rscu

Relative Synonymous Codon Usage (RSCU) pooled over several in-frame coding sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `pooled_rscu` |
| **Method ID** | `GenomeAnnotator.GetCodonUsage` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Pools codon counts over all `sequences` (each read in frame from position 0, non-overlapping; DNA or RNA,
case-insensitive, U read as T; a triplet with another symbol is skipped without shifting the frame; partial
trailing codon ignored, as in `codon_usage`) and computes RSCU = n·x / Σx per synonymous family of the chosen
NCBI translation table (Sharp, Tuohy & Mosurski 1986; CodonW `-totals`), delegating to
`CodonUsageAnalyzer.CalculateRscu`. Stop codons are excluded and an unobserved family gives 0.

## Core Documentation Reference

- Source: [GenomeAnnotator.cs#L2292](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs#L2292) — `GetCodonUsage(IEnumerable<string>, GeneticCode)` overload

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | string[] | Yes | In-frame coding DNA or RNA sequences (non-empty list) |
| `geneticCode` | integer | No | NCBI translation table number (default 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `rscu` | object | Sense codon (upper-case DNA) -> RSCU |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequences cannot be null or empty |
| 1001 | Unknown genetic code table: N. Supported NCBI tables: ... |

## Examples

### Example 1: Leucine family pooled over two sequences

`["CTTCTTCTG", "CTTTTA"]` -> CTT x3, CTG x1, TTA x1; Leu has 6 codons: CTT = 6·3/5 = 3.6, CTG = TTA = 1.2, other Leu codons 0.

**Response (Leu codons):**
```json
{ "rscu": { "CTT": 3.6, "CTG": 1.2, "TTA": 1.2, "CTC": 0.0, "CTA": 0.0, "TTG": 0.0 } }
```

## Performance

- **Time Complexity:** O(total length)
- **Space Complexity:** O(64)

## See Also

- [codon_usage](codon_usage.md) - raw codon counts of one sequence
