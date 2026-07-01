# find_branch_points

Find intron branch-point candidates.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_branch_points` |
| **Method ID** | `SpliceSitePredictor.FindBranchPoints` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans a sequence (T→U normalised) for branch-point candidates by scoring each 7-mer window against a
YNYURAC-like position weight matrix. For each window scoring at least `minScore`, the branch adenosine is
reported at `position = windowStart + 5`, with the 7-mer motif, PWM score, and confidence. The search can
be restricted to `[searchStart, searchEnd)` (`searchEnd = -1` scans to the end). The branch point is the
adenosine that attacks the 5' splice site during lariat formation.

## Core Documentation Reference

- Source: [SpliceSitePredictor.cs#L383](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs#L383)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | RNA/DNA sequence to scan |
| `searchStart` | integer | No | Inclusive start of search window (default 0) |
| `searchEnd` | integer | No | End of search window; -1 = to end (default -1) |
| `minScore` | number | No | Minimum PWM score [0,1] (default 0.5) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sites[].position` | integer | 0-based branch-adenosine index (window start + 5) |
| `sites[].type` | string | `Branch` |
| `sites[].motif` | string | 7-mer at the window start |
| `sites[].score` | number | PWM match score |
| `sites[].confidence` | number | Confidence |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Branch consensus present

A sequence containing a `UACUAAC`-like branch motif yields at least one `Branch` site (with `minScore = 0`).

### Example 2: Too short

A sequence shorter than 7 nt cannot contain a branch window:

**Response:**
```json
{ "sites": [] }
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(k) for k sites

## See Also

- [find_acceptor_sites](find_acceptor_sites.md) — 3' acceptor sites
- [predict_introns](predict_introns.md) — introns with flanking sites and branch points
