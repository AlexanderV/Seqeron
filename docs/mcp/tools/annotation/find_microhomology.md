# find_microhomology

Find the longest microhomology at a breakpoint junction.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_microhomology` |
| **Method ID** | `StructuralVariantAnalyzer.FindMicrohomology` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the longest sequence (up to `maxLength` nt) that is both a suffix of the left flank and a prefix of
the right flank — the microhomology at an SV breakpoint junction. Comparison is case-insensitive and the
returned sequence is upper-cased. Microhomology at a junction is a signature of microhomology-mediated
break-induced replication / end-joining.

## Core Documentation Reference

- Source: [StructuralVariantAnalyzer.cs#L1269](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs#L1269)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `leftFlank` | string | Yes | Left flanking sequence (5' side) |
| `rightFlank` | string | Yes | Right flanking sequence (3' side) |
| `maxLength` | integer | No | Maximum microhomology length to search (default 20) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `microhomologyLength` | integer | Length of the shared microhomology (0 if none) |
| `sequence` | string | The microhomology sequence (upper-cased) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Left flank cannot be null or empty |
| 1001 | Right flank cannot be null or empty |

## Examples

### Example 1: Shared CGT microhomology

`AAACGT` (suffix `CGT`) meets `CGTTTT` (prefix `CGT`):

**Response:**
```json
{ "microhomologyLength": 3, "sequence": "CGT" }
```

### Example 2: No microhomology

```json
{ "microhomologyLength": 0, "sequence": "" }
```

## Performance

- **Time Complexity:** O(maxLength · L) worst case
- **Space Complexity:** O(maxLength)

## See Also

- [cluster_split_reads](cluster_split_reads.md) — breakpoint detection
- [assemble_breakpoint_sequence](assemble_breakpoint_sequence.md) — junction assembly
