# compatible_enzymes

Enumerate all pairs of built-in enzymes that produce ligatable ends.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `compatible_enzymes` |
| **Method ID** | `RestrictionAnalyzer.FindCompatibleEnzymes` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns every unordered pair of built-in restriction enzymes whose ends can be ligated together. Two ends are compatible when **both are blunt**, or **both produce the same overhang type (5′/3′) and the same overhang sequence**. Each result carries the shared end (`"blunt"` or the overhang string).

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L502](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L502)

## Input Schema

_No parameters._

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `pairs` | array | Compatible enzyme pairs, each `{enzyme1, enzyme2, compatibleEnd}`. |

## Errors

_None._

## Examples

### Example 1: Sticky-end compatibility

BamHI (`GGATCC`, cut 1/5) and BglII (`AGATCT`, cut 1/5) both leave a 5′ `GATC` overhang, so they appear as a pair with `compatibleEnd = "GATC"`.

### Example 2: Blunt-end compatibility

All 10 blunt cutters are pairwise compatible, contributing C(10,2) = **45** pairs with `compatibleEnd = "blunt"`.

## Performance

- **Time Complexity:** O(n²) over the enzyme database (n = 39).
- **Space Complexity:** O(number of compatible pairs).

## See Also

- [enzymes_compatible](enzymes_compatible.md) - Test a single named pair.
- [blunt_cutters](blunt_cutters.md), [sticky_cutters](sticky_cutters.md)
