# codon_usage_statistics

Aggregate codon-usage statistics for a coding sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `codon_usage_statistics` |
| **Method ID** | `CodonUsageAnalyzer.GetStatistics` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Produces a one-shot codon-usage report: per-codon counts, RSCU (Relative Synonymous Codon Usage), ENC (Effective Number of Codons, Wright), total codons, GC% at codon positions 1/2/3, GC3s (GC% at the third position of synonymous codons only, excluding Met/Trp/stop), and overall GC = (GC1+GC2+GC3)/3.

## Core Documentation Reference

- Source: [CodonUsageAnalyzer.cs#L378](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs#L378)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Coding DNA sequence (frame 0), non-empty. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `codonCounts` | object | Codon → count. |
| `rscu` | object | Codon → RSCU. |
| `enc` | number | Effective Number of Codons (20..61). |
| `totalCodons` | integer | Number of complete codons. |
| `gc1` / `gc2` / `gc3` | number | GC% at codon position 1 / 2 / 3. |
| `gc3s` | number | GC% at the third position of synonymous codons. |
| `overallGc` | number | (GC1+GC2+GC3)/3. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Three-codon gene `ATGTTTAAA`

Codons ATG, TTT, AAA. GC1=0, GC2=0, GC3=33.33% (only the G in ATG), overall GC=11.11%. GC3s=0 (synonymous codons TTT, AAA have T/A at position 3). RSCU[TTT]=2.0 (Phe family, only TTT used).

## Performance

- **Time Complexity:** O(n) in sequence length.
- **Space Complexity:** O(k) in distinct codons.

## See Also

- [rscu](rscu.md), [effective_number_of_codons](effective_number_of_codons.md), [count_codons](count_codons.md)
