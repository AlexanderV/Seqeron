# enzymes_by_cut_length

List built-in restriction enzymes by recognition-sequence length.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `enzymes_by_cut_length` |
| **Method ID** | `RestrictionAnalyzer.GetEnzymesByCutLength` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns every enzyme in the built-in database whose recognition sequence is exactly `length` base pairs. Short recognition sites (4 bp) cut frequently, 6 bp are typical cloning enzymes, and 8 bp are rare cutters. Lengths with no matching enzyme return an empty list.

The built-in database contains **9** 4-cutters, **24** 6-cutters, **5** 8-cutters (plus one 13-bp interrupted palindrome, SfiI).

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L84](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L84)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `length` | integer | Yes | Recognition-sequence length in bp (positive). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `enzymes` | array | Enzymes with recognition length == `length`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Recognition-sequence length must be positive |

## Examples

### Example 1: `length = 6` → 24 six-cutters (EcoRI, BamHI, HindIII, …).

### Example 2: `length = 8` → 5 rare cutters (NotI, PacI, AscI, FseI, SwaI).

## See Also

- [get_enzyme](get_enzyme.md), [blunt_cutters](blunt_cutters.md), [sticky_cutters](sticky_cutters.md)
