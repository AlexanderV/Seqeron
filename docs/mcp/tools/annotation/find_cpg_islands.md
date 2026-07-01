# find_cpg_islands

Identify CpG islands using the Gardiner-Garden & Frommer criteria.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_cpg_islands` |
| **Method ID** | `EpigeneticsAnalyzer.FindCpGIslands` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Slides a window of `minLength` bases across the sequence and flags windows that satisfy the classic
Gardiner-Garden & Frommer (1987) CpG-island criteria: GC fraction ≥ `minGc` **and** CpG observed/expected
ratio ≥ `minCpGRatio`. Contiguous passing windows are merged; a merged run is reported only if it is at
least `minLength` long and the whole merged region still meets both thresholds. Each island reports its
0-based inclusive start, exclusive end, GC fraction, and CpG O/E ratio.

Defaults follow the standard definition: `minLength = 200`, `minGc = 0.5`, `minCpGRatio = 0.6`.

## Core Documentation Reference

- Source: [EpigeneticsAnalyzer.cs#L297](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs#L297)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence (min length: 1) |
| `minLength` | integer | No | Minimum island length in nt (default 200) |
| `minGc` | number | No | Minimum GC fraction 0..1 (default 0.5) |
| `minCpGRatio` | number | No | Minimum CpG O/E ratio (default 0.6) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `islands[].start` | integer | 0-based inclusive start |
| `islands[].end` | integer | Exclusive end |
| `islands[].gcContent` | number | GC fraction over the island (0..1) |
| `islands[].cpGRatio` | number | CpG observed/expected ratio over the island |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Perfect CpG island

**User Prompt:**
> Find CpG islands in a 200 bp run of CG repeats.

**Response:**
```json
{ "islands": [ { "start": 0, "end": 200, "gcContent": 1.0, "cpGRatio": 2.0 } ] }
```

Every window is 100% GC with CpG O/E = 2.0, forming one island spanning 0–200.

### Example 2: Below minimum length

A 100 bp sequence is shorter than the default 200 bp minimum:

**Response:**
```json
{ "islands": [] }
```

## Performance

- **Time Complexity:** O(n · minLength) worst case (per-window GC/O/E recomputation)
- **Space Complexity:** O(minLength)

## See Also

- [cpg_observed_expected](cpg_observed_expected.md) — CpG O/E ratio for a single sequence
- [find_cpg_sites](find_cpg_sites.md) — CpG dinucleotide positions
