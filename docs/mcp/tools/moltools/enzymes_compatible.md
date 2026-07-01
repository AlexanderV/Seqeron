# enzymes_compatible

Check whether two named enzymes produce ligatable ends.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `enzymes_compatible` |
| **Method ID** | `RestrictionAnalyzer.AreCompatible` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns `true` when the two named enzymes produce ligatable ends: both blunt, or both leaving the same overhang type (5′/3′) and overhang sequence. An unknown enzyme name resolves to `false` (no error). Names are case-insensitive.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L524](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L524)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `enzyme1_name` | string | Yes | First enzyme name (non-blank). |
| `enzyme2_name` | string | Yes | Second enzyme name (non-blank). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `compatible` | boolean | True if the two enzymes' ends can be ligated. |

## Errors

| Code | Message |
|------|---------|
| 1001 | First enzyme name cannot be null or blank |
| 1002 | Second enzyme name cannot be null or blank |

## Examples

### Example 1: `EcoRV` + `SmaI` (both blunt) → `true`.

### Example 2: `EcoRI` (5′ AATT) + `BamHI` (5′ GATC) → `false` (different overhang sequences).

## See Also

- [compatible_enzymes](compatible_enzymes.md) - Enumerate all compatible pairs.
