# cpg_observed_expected

Compute the CpG observed/expected ratio for a sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `cpg_observed_expected` |
| **Method ID** | `EpigeneticsAnalyzer.CalculateCpGObservedExpected` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the CpG observed/expected (O/E) ratio, a core Gardiner-Garden & Frommer criterion for CpG-island
detection:

```
O/E = observed_CpG / expected_CpG
expected_CpG = (C_count * G_count) / sequence_length
```

Regions where CpG dinucleotides occur roughly as often as expected from base composition (O/E ≈ 1) are
CpG islands; most of the genome is CpG-depleted (O/E well below 1). The ratio is `0` when the sequence is
shorter than 2 bp or contains no C/G (expected = 0).

## Core Documentation Reference

- Source: [EpigeneticsAnalyzer.cs#L270](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs#L270)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `ratio` | number | Observed/expected CpG ratio; `0` when expected is 0 |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: CpG-rich sequence

**User Prompt:**
> What's the CpG O/E ratio of CGCGCG?

**Response:**
```json
{ "ratio": 2.0 }
```

Length 6, C=3, G=3, observed CpG = 3, expected = (3×3)/6 = 1.5, ratio = 3 / 1.5 = 2.0.

### Example 2: No CpG dinucleotide

```json
{ "tool": "cpg_observed_expected", "arguments": { "sequence": "ATGCATGC" } }
```

**Response:**
```json
{ "ratio": 0.0 }
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1)

## See Also

- [find_cpg_islands](find_cpg_islands.md) — CpG islands using the O/E criterion
- [find_cpg_sites](find_cpg_sites.md) — CpG dinucleotide positions
