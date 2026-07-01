# find_split_reads

Find split reads from soft-clipped CIGAR alignments.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_split_reads` |
| **Method ID** | `StructuralVariantAnalyzer.FindSplitReads` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Parses each alignment's CIGAR string for soft-clips (`S` operations) and emits a split read for every clip
at least `minClipLength` long. Each split read records the primary alignment position, the inferred
supplementary position (the aligned length past the primary for a right clip), the clip length, and the
clipped sequence. Split reads pinpoint breakpoints at single-base resolution.

## Core Documentation Reference

- Source: [StructuralVariantAnalyzer.cs#L399](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs#L399)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `alignments` | array | Yes | Reads `{ readId, chromosome, position, cigar, sequence }` |
| `minClipLength` | integer | No | Minimum soft-clip length to call a split read (default 20) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `reads[].readId` | string | Read identifier |
| `reads[].chromosome` | string | Chromosome |
| `reads[].primaryPosition` | integer | Primary alignment position |
| `reads[].supplementaryPosition` | integer | Inferred supplementary position |
| `reads[].clipLength` | integer | Soft-clip length |
| `reads[].clippedSequence` | string | Clipped sequence |

## Errors

| Code | Message |
|------|---------|
| 1001 | alignments cannot be null |

## Examples

### Example 1: Left soft-clip

A `30S70M` read yields one split read with a 30 nt clip at the primary position (clipped sequence = first
30 bases).

### Example 2: Fully matched read

A `100M` read has no soft-clips:

**Response:**
```json
{ "reads": [] }
```

## Performance

- **Time Complexity:** O(Σ|cigar|) over all reads
- **Space Complexity:** O(k) for k split reads

## See Also

- [cluster_split_reads](cluster_split_reads.md) — cluster split reads into breakpoints
- [find_discordant_pairs](find_discordant_pairs.md) — read-pair SV evidence
