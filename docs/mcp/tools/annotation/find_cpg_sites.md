# find_cpg_sites

Return all CpG dinucleotide start positions in a sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_cpg_sites` |
| **Method ID** | `EpigeneticsAnalyzer.FindCpGSites` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans a nucleotide sequence (case-insensitively) and returns the 0-based start index of every CpG
dinucleotide — a cytosine (`C`) immediately followed by a guanine (`G`). CpG sites are the principal
substrate for DNA methylation in vertebrates and the building blocks of CpG islands.

## Core Documentation Reference

- Source: [EpigeneticsAnalyzer.cs#L148](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs#L148)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence to scan (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `positions` | integer[] | 0-based start index of each CG dinucleotide |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Two CpG sites

**User Prompt:**
> Where are the CpG sites in ACGTACGT?

**Expected Tool Call:**
```json
{ "tool": "find_cpg_sites", "arguments": { "sequence": "ACGTACGT" } }
```

**Response:**
```json
{ "positions": [1, 5] }
```

The C at index 1 is followed by G at 2; the C at index 5 is followed by G at 6.

### Example 2: No CpG sites

```json
{ "tool": "find_cpg_sites", "arguments": { "sequence": "AAATTT" } }
```

**Response:**
```json
{ "positions": [] }
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(k) for k CpG sites

## See Also

- [find_cpg_islands](find_cpg_islands.md) — CpG islands via Gardiner-Garden & Frommer criteria
- [cpg_observed_expected](cpg_observed_expected.md) — CpG observed/expected ratio
- [find_methylation_sites](find_methylation_sites.md) — CpG/CHG/CHH methylation contexts
