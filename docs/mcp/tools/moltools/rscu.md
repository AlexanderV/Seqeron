# rscu

Relative Synonymous Codon Usage per codon.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `rscu` |
| **Method ID** | `CodonUsageAnalyzer.CalculateRscu` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes RSCU (Sharp, Tuohy & Mosurski 1986) for each codon: `RSCU = observed / (total_synonymous / n_synonymous)`, i.e. observed count divided by the count expected if all synonymous codons of that amino acid were used equally. `RSCU = 1` means no bias, `> 1` over-represented, `< 1` under-represented. Codons of amino acids not present in the sequence get RSCU 0.

## Core Documentation Reference

- Source: [CodonUsageAnalyzer.cs#L80](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs#L80)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Coding DNA sequence (frame 0), non-empty. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `rscu` | object | Codon → RSCU value. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Biased usage

`TTTTTTTTC` → TTT×2, TTC×1; Phe expected = 3/2 = 1.5 → `RSCU[TTT] = 1.333`, `RSCU[TTC] = 0.667`.

### Example 2: Unbiased usage

`TTTTTC` → TTT×1, TTC×1 → `RSCU = 1.0` for each.

## See Also

- [count_codons](count_codons.md), [codon_usage_statistics](codon_usage_statistics.md), [codon_adaptation_index](codon_adaptation_index.md)
