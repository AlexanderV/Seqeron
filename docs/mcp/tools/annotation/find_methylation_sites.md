# find_methylation_sites

Find candidate cytosine methylation sites and classify their sequence context.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_methylation_sites` |
| **Method ID** | `EpigeneticsAnalyzer.FindMethylationSites` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans a nucleotide sequence and reports every cytosine that falls into a recognised methylation context —
**CpG** (C followed by G), **CHG** (C, H, G) or **CHH** (C, H, H), where H = A/C/T. Each site carries its
0-based position, its context type, and an up-to-3-base context window (2 bases for a terminal CpG). Because
sequence-only input has no bisulfite read evidence, `methylationLevel` and `coverage` are `0`; measured
levels come from [methylation_from_bisulfite](methylation_from_bisulfite.md) or caller-supplied sites.

## Core Documentation Reference

- Source: [EpigeneticsAnalyzer.cs#L229](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs#L229)
- Context rules: [EpigeneticsAnalyzer.cs#L179](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs#L179)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence to scan (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sites[].position` | integer | 0-based cytosine position |
| `sites[].type` | string | `CpG`, `CHG`, or `CHH` |
| `sites[].context` | string | Up-to-3-base context window starting at the cytosine |
| `sites[].methylationLevel` | number | `0` for sequence-only input |
| `sites[].coverage` | integer | `0` for sequence-only input |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Mixed contexts

**User Prompt:**
> Classify the methylation contexts in CAGCTTCG.

**Response:**
```json
{
  "sites": [
    { "position": 0, "type": "CHG", "context": "CAG", "methylationLevel": 0, "coverage": 0 },
    { "position": 3, "type": "CHH", "context": "CTT", "methylationLevel": 0, "coverage": 0 },
    { "position": 6, "type": "CpG", "context": "CG", "methylationLevel": 0, "coverage": 0 }
  ]
}
```

### Example 2: No classifiable cytosine

```json
{ "tool": "find_methylation_sites", "arguments": { "sequence": "ATATAT" } }
```

**Response:**
```json
{ "sites": [] }
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(k) for k methylation sites

## See Also

- [find_cpg_sites](find_cpg_sites.md) — CpG dinucleotide positions only
- [methylation_from_bisulfite](methylation_from_bisulfite.md) — measured methylation from reads
- [methylation_profile](methylation_profile.md) — aggregate sites into a profile
