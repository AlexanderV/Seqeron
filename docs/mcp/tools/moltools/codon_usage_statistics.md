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
| `genetic_code` | integer | No | NCBI translation table number (default `1` = Standard; e.g. `2` vertebrate mitochondrial, `11` bacterial/plastid). Defines the synonymous families and stops (CodonW `-code`). One of 1–6, 9–16, 21–33. |

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
| 1002 | Unknown genetic code table: N. Supported NCBI tables: 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33. |

## Examples

### Example 1: Three-codon gene `ATGTTTAAA`

Codons ATG, TTT, AAA. GC1=0, GC2=0, GC3=33.33% (only the G in ATG), overall GC=11.11%. GC3s=0 (synonymous codons TTT, AAA have T/A at position 3). RSCU[TTT]=2.0 (Phe family, only TTT used).

### Example 2: Non-standard genetic code (`genetic_code = 2`)

`ATAATG`: under the default table 1 ATA (Ile, T at position 3) is the only synonymous codon, so GC3s = 0; under NCBI table 2 Met = {ATA, ATG} is a two-codon family, so both codons are synonymous and GC3s = 50% (CodonW 1.4.4 `-gc3s`: 0.000 / 0.500). Counts and GC1/GC2/GC3 are unchanged.

## Performance

- **Time Complexity:** O(n) in sequence length.
- **Space Complexity:** O(k) in distinct codons.

## See Also

- [rscu](rscu.md), [effective_number_of_codons](effective_number_of_codons.md), [count_codons](count_codons.md)
