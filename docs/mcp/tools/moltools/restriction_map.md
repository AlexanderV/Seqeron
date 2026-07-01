# restriction_map

Build a restriction map of a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `restriction_map` |
| **Method ID** | `RestrictionAnalyzer.CreateMap` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Maps a DNA sequence: reports all forward+reverse restriction sites, the sites grouped by enzyme (positions), the total number of forward-strand sites, the **unique cutters** (enzymes with exactly one forward-strand site — ideal for cloning), and the **non-cutters** from the queried set. An empty `enzyme_names` list considers every built-in enzyme.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L463](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L463)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to map (non-empty). |
| `enzyme_names` | string[] | Yes | Enzymes to consider; empty = all built-in enzymes. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sequenceLength` | integer | Sequence length. |
| `sites` | array | All sites (both strands). |
| `sitesByEnzyme` | object | Enzyme → sorted positions. |
| `totalSites` | integer | Forward-strand site count. |
| `uniqueCutters` | string[] | Enzymes with one forward-strand site. |
| `nonCutters` | string[] | Queried enzymes with no site. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: `AAAGAATTCAAA` with `[EcoRI, BamHI]`

EcoRI cuts once (unique cutter, positions `[3, 3]` for both strands); BamHI does not cut → non-cutter. `totalSites = 1`.

## See Also

- [find_restriction_sites](find_restriction_sites.md), [restriction_digest](restriction_digest.md)
