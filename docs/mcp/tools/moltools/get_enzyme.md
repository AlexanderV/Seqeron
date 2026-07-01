# get_enzyme

Look up a built-in restriction enzyme by name.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `get_enzyme` |
| **Method ID** | `RestrictionAnalyzer.GetEnzyme` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Case-insensitive lookup of a single enzyme in the built-in Type II restriction-enzyme database. Returns the enzyme's recognition sequence, forward/reverse cut positions and source organism, or `enzyme = null` when the name is not present.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L76](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L76)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `name` | string | Yes | Enzyme name (non-blank). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `enzyme` | object \| null | The enzyme record, or `null` if not found. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Enzyme name cannot be null or blank |

## Examples

### Example 1: EcoRI

`get_enzyme("EcoRI")` → `{ name: "EcoRI", recognitionSequence: "GAATTC", cutPositionForward: 1, cutPositionReverse: 5, organism: "Escherichia coli" }`. Lookup is case-insensitive (`ecori` works too).

### Example 2: Unknown enzyme → `{ enzyme: null }`.

## See Also

- [enzymes_by_cut_length](enzymes_by_cut_length.md), [find_restriction_sites](find_restriction_sites.md)
