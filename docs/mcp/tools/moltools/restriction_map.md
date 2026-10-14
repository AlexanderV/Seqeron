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

Maps a DNA sequence: reports all forward+reverse restriction sites, the sites grouped by enzyme (positions), the total number of forward-strand sites, the **unique cutters** (enzymes with exactly one forward-strand site — ideal for cloning), and the **non-cutters**.

The *considered set* is the enzyme set the map is built over: the names in `enzyme_names`, resolved case-insensitively and de-duplicated, or — when `enzyme_names` is empty — every built-in enzyme. `sitesByEnzyme` and `nonCutters` always **partition** that set, and both are keyed by the table's own canonical enzyme name whatever casing the caller used. Site coordinates, the two records per palindromic site and the `Linear`/`Circular` rules are as documented for [find_restriction_sites](find_restriction_sites.md).

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L803](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L803)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to map (non-empty). |
| `enzyme_names` | string[] | Yes | Enzymes to consider (case-insensitive, de-duplicated); empty = all built-in enzymes. |
| `topology` | string | No | `Linear` (default) or `Circular`. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sequenceLength` | integer | Sequence length. |
| `sites` | array | All sites (both strands). |
| `sitesByEnzyme` | object | Enzyme → sorted positions. |
| `totalSites` | integer | Forward-strand site count. |
| `uniqueCutters` | string[] | Enzymes with one forward-strand site. |
| `nonCutters` | string[] | Considered enzymes with no site, under their canonical name. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: `AAAGAATTCAAA` with `[EcoRI, BamHI]`

EcoRI cuts once (unique cutter, positions `[3, 3]` for both strands); BamHI does not cut → non-cutter. `totalSites = 1`. The same holds for `[ecori, bamhi]`: `sitesByEnzyme` is keyed `EcoRI` and `nonCutters` is `["BamHI"]`.

### Example 2: `AAAGAATTCAAA` with `enzyme_names: []`

The considered set is the whole built-in table, so `sitesByEnzyme` is `{ "EcoRI": [3, 3] }` and `nonCutters` lists the other 38 built-in enzymes.

## See Also

- [find_restriction_sites](find_restriction_sites.md), [restriction_digest](restriction_digest.md)
