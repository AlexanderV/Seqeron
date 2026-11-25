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

Scans a sequence (upper-cased, T→U) for branch-point candidates. Every motif placement whose start index
lies in [`searchStart`, `searchEnd`] (inclusive; `searchEnd = -1` = to the end) is scored; placements with
a non-ACGTU character are skipped.

- **`model = "human"` (default):** 5-nt `yUnAy` (positions -3..+1, branch A at 0; Gao et al. 2008),
  conservation-weighted match fraction with the Gao et al. Table 4 frequencies (y@-3 0.790, U@-2 0.746,
  A@0 0.923, y@+1 0.751; -1 uninformative); perfect `yUnAy` = 1.0.
- **`model = "yeast"`:** 7-nt `UACUAAC` (S. cerevisiae; branch A = sixth nt), fraction of matching positions.

`position` is the 0-based index of the branch-adenosine position (motif start + 3 for human, + 5 for
yeast); confidence = clamp((score − 0.5)/0.5, 0, 1).

## Core Documentation Reference

- Source: [SpliceSitePredictor.cs#L768](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs#L768)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | RNA/DNA sequence to scan |
| `searchStart` | integer | No | Smallest motif start index (inclusive, default 0) |
| `searchEnd` | integer | No | Largest motif start index (inclusive); -1 = to end (default -1) |
| `minScore` | number | No | Minimum score [0,1] (default 0.5) |
| `model` | string | No | `human` (yUnAy, default) or `yeast` (UACUAAC) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sites[].position` | integer | 0-based branch-adenosine index |
| `sites[].type` | string | `Branch` |
| `sites[].motif` | string | Scored motif (5 nt human, 7 nt yeast; RNA alphabet) |
| `sites[].score` | number | Match score [0,1] |
| `sites[].confidence` | number | Confidence |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Unknown model '...'; expected 'human' or 'yeast'. |

## Examples

### Example 1: Human yUnAy

`GGGGGUACUAACGGGGGUUUUUCAGG`, `minScore = 0.75`:

**Response:**
```json
{ "sites": [ { "position": 10, "type": "Branch", "motif": "CUAAC", "score": 1.0 },
             { "position": 23, "type": "Branch", "motif": "UUCAG", "score": 0.7660436137071651 } ] }
```

### Example 2: Yeast UACUAAC

Same sequence, `model = "yeast"`, `minScore = 0.75` → one site, position 10, motif `UACUAAC`, score 1.0.

### Example 3: Too short

A sequence shorter than the motif (5 nt human) has no placement: `{ "sites": [] }`.

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(k) for k sites

## See Also

- [find_acceptor_sites](find_acceptor_sites.md) — 3' acceptor sites
- [predict_introns](predict_introns.md) — introns with flanking sites and branch points
