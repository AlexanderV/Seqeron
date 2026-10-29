# sticky_cutters

List the restriction enzymes that produce sticky (cohesive) ends (common set by default, or all 754 REBASE enzymes).

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `sticky_cutters` |
| **Method ID** | `RestrictionAnalyzer.GetStickyCutters` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns every enzyme of `enzyme_set` — `"common"` (default, the 39 common laboratory enzymes) or `"rebase"` (all 754 REBASE enzymes with a defined cut) — whose cut is staggered (`CutPositionForward != CutPositionReverse`), leaving a 5′ or 3′ single-stranded overhang. These cohesive ends anneal to complementary overhangs and are used for directional / sticky-end cloning. The complement of `blunt_cutters`.

The common set has **29** sticky cutters (all 39 enzymes minus the 10 blunt cutters): **21** leave a
5′ overhang and **8** a 3′ overhang. Each returned record carries `overhangType`
(`FivePrime` / `ThreePrime`); for the overhang *sequence* use `get_enzyme`, and to test whether two
of them can be ligated use `enzymes_compatible`. All 39 entries were checked against Biopython
`is_5overhang()` / `is_3overhang()` on 2026-10-14: 39/39 agreement.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L175](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L175)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `enzyme_set` | string | No | `"common"` (default, the 39 common laboratory enzymes) or `"rebase"` (all 754 REBASE enzymes with a defined cut). Case-insensitive. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `enzymes` | array | Sticky-cutting enzymes. |

## Errors

| Code | Message |
|------|---------|
| 1010 | Unknown enzyme set '&lt;value&gt;'; expected "common" or "rebase" (`enzyme_set`) |

## Examples

### Example 1: EcoRI (`GAATTC`, cut 1/5) is a sticky cutter leaving a 5′ `AATT` overhang (`overhangType: "FivePrime"`).

### Example 2: PstI (`CTGCAG`, cut 5/1) leaves a 3′ `TGCA` overhang (`overhangType: "ThreePrime"`).

### Example 3: The list contains 29 enzymes; blunt cutters (EcoRV, SmaI, HincII, …) are excluded.

### Example: the REBASE set

`sticky_cutters(enzyme_set: "rebase")` → **617** enzymes — **388** `FivePrime`, **229** `ThreePrime`
(Biopython 1.88 `is_5overhang()` / `is_3overhang()` over the 754) — starting `AarI` (`CACCTGC`,
11/15, Type IIS), `AasI` (`GACNNNNNNGTC`, 7/5), `AatII`, … The default call still returns the 29
common sticky cutters.

### Example: unknown enzyme set

`sticky_cutters(enzyme_set: "all")` → error 1010.

## See Also

- [blunt_cutters](blunt_cutters.md), [compatible_enzymes](compatible_enzymes.md), [enzymes_compatible](enzymes_compatible.md), [get_enzyme](get_enzyme.md)
