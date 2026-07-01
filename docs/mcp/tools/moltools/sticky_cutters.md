# sticky_cutters

List all built-in restriction enzymes that produce sticky (cohesive) ends.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `sticky_cutters` |
| **Method ID** | `RestrictionAnalyzer.GetStickyCutters` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns every enzyme in the built-in database whose cut is staggered (`CutPositionForward != CutPositionReverse`), leaving a 5′ or 3′ single-stranded overhang. These cohesive ends anneal to complementary overhangs and are used for directional / sticky-end cloning.

The database has **29** sticky cutters (all 39 enzymes minus the 10 blunt cutters).

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L114](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L114)

## Input Schema

_No parameters._

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `enzymes` | array | Sticky-cutting enzymes. |

## Errors

_None._

## Examples

### Example 1: EcoRI (`GAATTC`, cut 1/5) is a sticky cutter leaving a 5′ `AATT` overhang.

### Example 2: The list contains 29 enzymes; blunt cutters (EcoRV, SmaI, …) are excluded.

## See Also

- [blunt_cutters](blunt_cutters.md), [compatible_enzymes](compatible_enzymes.md)
