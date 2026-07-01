# find_restriction_sites

Find restriction sites for named enzymes on both strands.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `find_restriction_sites` |
| **Method ID** | `RestrictionAnalyzer.FindSites` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans a DNA sequence on both strands for the recognition sites of the named built-in enzymes. IUPAC degenerate codes in recognition sequences are matched against ACGT input. Each hit reports its forward-strand position, strand, cut position, and matched sequence. A palindromic enzyme reports both a forward and a reverse site at the same position.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L207](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L207)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to scan (non-empty). |
| `enzyme_names` | string[] | Yes | Enzyme names to look for (≥1). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sites` | array | Sites, each `{position, isForwardStrand, cutPosition, recognizedSequence}`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | At least one enzyme name is required |
| 4001 | Unknown enzyme name |

## Examples

### Example 1: EcoRI

`AAAGAATTCAAA` + `EcoRI` → one forward site at position 3 (cut 4) and one reverse site at position 3 (cut 8).

### Example 2: No site → `sites: []`.

## See Also

- [find_all_restriction_sites](find_all_restriction_sites.md), [restriction_digest](restriction_digest.md), [restriction_map](restriction_map.md)
