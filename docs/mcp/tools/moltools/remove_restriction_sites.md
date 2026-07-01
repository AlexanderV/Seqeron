# remove_restriction_sites

Synonymously eliminate restriction sites from a coding sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `remove_restriction_sites` |
| **Method ID** | `CodonOptimizer.RemoveRestrictionSites` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Rewrites codons to remove the listed restriction recognition sequences while preserving the encoded protein (only synonymous swaps are used). Site strings may be DNA or RNA (converted to RNA internally); the output is RNA-alphabet. Sites for which no synonymous alternative eliminates the match are left in place.

## Core Documentation Reference

- Source: [CodonOptimizer.cs#L531](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs#L531)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `coding_sequence` | string | Yes | Coding sequence (DNA or RNA), non-empty. |
| `restriction_sites` | string[] | Yes | Recognition sequences to remove (≥1). |
| `target_organism` | object | Yes | Preset id or inline custom table. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `optimizedSequence` | string | RNA-alphabet sequence with sites removed where possible. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Coding sequence cannot be null or empty |
| 1002 | At least one restriction site is required |

## Examples

### Example 1: Remove EcoRI

`GAATTC` (Glu-Phe) with site `GAATTC` → e.g. `GAGUUC` (same protein, no `GAAUUC`).

### Example 2: No site present → sequence returned unchanged (as RNA).

## See Also

- [optimize_codons](optimize_codons.md), [find_restriction_sites](find_restriction_sites.md)
