# assemble_breakpoint_sequence

Heuristically assemble a breakpoint-junction sequence from split-read clipped fragments.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `assemble_breakpoint_sequence` |
| **Method ID** | `StructuralVariantAnalyzer.AssembleBreakpointSequence` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Given the split reads that support a structural-variant breakpoint, returns a representative
junction sequence. The current heuristic selects the **clipped sequence of the read with the
largest `clipLength`** (the longest soft-clip carries the most junction evidence). If no reads are
supplied, the result is `null`. `minOverlap` is reserved for future overlap-layout assembly and
does not change the selection today.

## Core Documentation Reference

- Source: [StructuralVariantAnalyzer.cs#L1252](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs#L1252)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `splitReads` | array | Yes | Split reads `{ readId, chromosome, primaryPosition, supplementaryPosition, clipLength, clippedSequence }` |
| `minOverlap` | integer | No | Minimum overlap (nt) between assembled fragments (default `10`) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sequence` | string \| null | Assembled breakpoint sequence, or `null` when no reads were supplied |

## Errors

| Code | Message |
|------|---------|
| 1001 | splitReads cannot be null |

## Examples

### Example 1: Longest split read wins

Three split reads with clip lengths 20, 50, 30 — the 50 nt read's clipped sequence is returned.

**Response:**
```json
{ "sequence": "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC" }
```

### Example 2: No reads

**Response:**
```json
{ "sequence": null }
```

## Performance

- **Time Complexity:** O(n) over the split reads
- **Space Complexity:** O(n)

## See Also

- [find_split_reads](find_split_reads.md) - Detect split reads from soft-clips
- [cluster_split_reads](cluster_split_reads.md) - Cluster split reads into breakpoints
- [find_microhomology](find_microhomology.md) - Microhomology at a junction
