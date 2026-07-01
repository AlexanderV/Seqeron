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

Performs a comprehensive restriction scan: every enzyme in the built-in database is searched against both strands of the sequence. Use this when the caller does not know which enzymes to look for. Palindromic enzymes report both strands at the same position.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L221](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L221)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to scan (non-empty). |

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

## See Also

- [find_restriction_sites](find_restriction_sites.md), [restriction_map](restriction_map.md)
