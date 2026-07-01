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

Analyses the local context of a miRNA target site within an mRNA. It measures the AU content of
a flanking window (AU-rich context favours miRNA targeting), flags whether the site lies near the
transcript start (first 15%) or end (last 15%), and returns a composite context favourability
score in `[0,1]`: `0.5 × auContent`, plus a `0.3` bonus when the site is neither near-start nor
near-end, capped at `1`.

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
| `nearStart` | boolean | Site starts within the first 15% of the transcript |
| `nearEnd` | boolean | Site ends within the last 15% of the transcript |
| `contextScore` | number | Favourability score in [0,1] |

## Errors

| Code | Message |
|------|---------|
| 1001 | mRNA sequence cannot be null or empty |
| 1002 | Target indices are out of range for the sequence |

## Examples

### Example 1: AU-rich mid-transcript site

**User Prompt:**
> Analyse the target context at positions 8–11 of a 20-nt poly-A mRNA.

**Expected Tool Call:**
```json
{
  "tool": "analyze_target_context",
  "arguments": { "mRnaSequence": "AAAAAAAAAAAAAAAAAAAA", "targetStart": 8, "targetEnd": 11, "contextWindow": 30 }
}
```

**Response:**
```json
{ "auContent": 1.0, "nearStart": false, "nearEnd": false, "contextScore": 0.8 }
```

### Example 2: GC-rich mid-transcript site

**Response:**
```json
{ "auContent": 0.0, "nearStart": false, "nearEnd": false, "contextScore": 0.3 }
```

## Performance

- **Time Complexity:** O(w) where w is the context window width
- **Space Complexity:** O(w)

## See Also

- [site_accessibility](site_accessibility.md) - Local secondary-structure accessibility
- [find_mirna_target_sites](find_mirna_target_sites.md) - Scan an mRNA for target sites
