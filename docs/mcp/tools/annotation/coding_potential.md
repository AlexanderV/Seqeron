# coding_potential

Score the coding potential of a DNA sequence using the CPAT hexamer usage-bias log-likelihood.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `coding_potential` |
| **Method ID** | `GenomeAnnotator.CalculateCodingPotential` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Implements the CPAT hexamer usage-bias score (`FrameKmer.kmer_ratio`, frame 0; Wang et al. 2013).
The sequence is scanned in-frame in 6-nt windows stepping by 3. For each hexamer present in both the
coding and non-coding frequency tables the score adds `ln(coding/noncoding)`; a hexamer found only in
the coding table adds `+1`, one found only in the non-coding table subtracts `1`; hexamers missing
from either table (or with both counts zero) are skipped. The result is the sum divided by the number
of scored hexamers. **Positive** values indicate a coding sequence, **negative** values a non-coding
one. Sequences shorter than one hexamer, or with no scorable hexamer, return `0`.

## Core Documentation Reference

- Source: [GenomeAnnotator.cs#L794](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs#L794)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to score |
| `codingHexamerFrequencies` | object | Yes | In-frame 6-mer → count/proportion from a coding training set |
| `noncodingHexamerFrequencies` | object | Yes | In-frame 6-mer → count/proportion from a non-coding training set |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `score` | number | Mean per-hexamer log-likelihood ratio (positive = coding, negative = non-coding) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Two in-frame hexamers

`ATGAAACCC` with coding `{ATGAAA:8, AAACCC:2}` and non-coding `{ATGAAA:2, AAACCC:4}` →
mean of `ln(8/2)` and `ln(2/4)` = 0.3465735902799726.

**Response:**
```json
{ "score": 0.34657359027997264 }
```

### Example 2: Single in-both hexamer

`ATGAAA` with coding `{ATGAAA:4}` and non-coding `{ATGAAA:1}` → `ln(4)` = 1.3862943611198906.

## Performance

- **Time Complexity:** O(n) hexamer scan
- **Space Complexity:** O(1) beyond the input tables

## See Also

- [find_orfs](find_orfs.md) - ORF detection
- [predict_genes](predict_genes.md) - ORF-based gene prediction
