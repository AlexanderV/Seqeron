# find_all_restriction_sites

Find sites for every built-in restriction enzyme.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `find_all_restriction_sites` |
| **Method ID** | `RestrictionAnalyzer.FindAllSites` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Performs a comprehensive restriction scan: every enzyme in the built-in database is searched against both strands of the sequence. Use this when the caller does not know which enzymes to look for. Coordinates, the two records reported per palindromic site, overlapping-site reporting, the `Linear` cut-feasibility rule and the `Circular` origin-spanning search are all exactly as documented for [find_restriction_sites](find_restriction_sites.md).

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L330](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L330)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to scan (non-empty). |
| `topology` | string | No | `Linear` (default) or `Circular`. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sites` | array | All restriction sites found across every enzyme. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: `AAAGAATTCAAA` → includes the EcoRI site at position 3 (forward cut 4, plus reverse strand).

### Example 2: `AAAAAAAA` → `sites: []` (no recognition sequence matches).

### Example 3: `TTCAAAAAAAGAA` + `topology: Circular` → includes the origin-spanning EcoRI site at position 10 (forward cut 11), which the linear scan does not report.

## See Also

- [find_restriction_sites](find_restriction_sites.md), [restriction_map](restriction_map.md)
