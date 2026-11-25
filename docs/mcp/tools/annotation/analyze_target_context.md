# analyze_target_context

Compute AU content and positional context score around a miRNA target site.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `analyze_target_context` |
| **Method ID** | `MiRnaAnalyzer.AnalyzeTargetContext` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Analyses the local context of a miRNA target site within a 3'UTR (the `mRnaSequence`). It returns the
TargetScan 5.0 position-weighted local A/U fraction of the flanking window (Grimson et al. 2007; 30 nt each
side, site excluded, 7mer-m8 weights), flags whether the site starts within the first 15 nt (1-based start < 15,
`too_close`) or ends in the last 15% of the sequence, and returns the UTR-context part of the fitted context score:
`round3(-0.5 x auContent + 0.269) + round3(0.000091 x min(d5, d3, 1500) - 0.037)` in log2 units (negative =
repressed; 0 when near the start). The site-type and 3' pairing terms need the miRNA (library
`ScoreTargetSiteGrimson`).

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L2591](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L2591)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `mRnaSequence` | string | Yes | mRNA nucleotide sequence (min length: 1) |
| `targetStart` | integer | Yes | 0-based inclusive start of the target site |
| `targetEnd` | integer | Yes | 0-based inclusive end of the target site |
| `contextWindow` | integer | No | Flanking context window (nt) on each side (default: 30) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `auContent` | number | Fraction of A/U bases in the context window (0-1) |
| `nearStart` | boolean | Site starts within the first 15 nt of the 3'UTR (1-based start < 15) |
| `nearEnd` | boolean | Site ends within the last 15% of the transcript |
| `contextScore` | number | Local-AU + position contribution of the TargetScan 5.0 context score (log2 units, negative = repressed); 0 when `nearStart` |

## Errors

| Code | Message |
|------|---------|
| 1001 | mRNA sequence cannot be null or empty |
| 1002 | Target indices are out of range for the sequence |

## Examples

### Example 1: AU-rich mid-transcript site

**User Prompt:**
> Analyse the target context at positions 40–50 of a 100-nt poly-A mRNA.

**Expected Tool Call:**
```json
{
  "tool": "analyze_target_context",
  "arguments": { "mRnaSequence": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "targetStart": 40, "targetEnd": 50, "contextWindow": 30 }
}
```

**Response:**
```json
{ "auContent": 1.0, "nearStart": false, "nearEnd": false, "contextScore": -0.264 }
```

### Example 2: GC-rich mid-transcript site

**Response:**
```json
{ "auContent": 0.0, "nearStart": false, "nearEnd": false, "contextScore": 0.236 }
```

## Performance

- **Time Complexity:** O(w) where w is the context window width
- **Space Complexity:** O(w)

## See Also

- [site_accessibility](site_accessibility.md) - Local secondary-structure accessibility
- [find_mirna_target_sites](find_mirna_target_sites.md) - Scan an mRNA for target sites
