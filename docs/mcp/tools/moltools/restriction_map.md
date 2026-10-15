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

The *considered set* is the enzyme set the map is built over: the names in `enzyme_names`, resolved case-insensitively and de-duplicated, or — when `enzyme_names` is empty — every enzyme of `enzyme_set`: `"common"` (default, the 39 common laboratory enzymes) or `"rebase"` (all 754 REBASE enzymes with a defined cut). Named enzymes may be any of the 754; `enzyme_set` is ignored when names are given. `sitesByEnzyme` and `nonCutters` always **partition** that set, and both are keyed by the table's own canonical enzyme name whatever casing the caller used. Site coordinates, the two records per palindromic site and the `Linear`/`Circular` rules are as documented for [find_restriction_sites](find_restriction_sites.md).

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L1344](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L1344)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to map (non-empty). |
| `enzyme_names` | string[] | Yes | Enzymes to consider (case-insensitive, de-duplicated; any of the 754 REBASE enzymes with a defined cut); empty = every enzyme of `enzyme_set`. |
| `topology` | string | No | `Linear` (default) or `Circular`. |
| `enzyme_set` | string | No | Table used when `enzyme_names` is empty: `"common"` (default, the 39) or `"rebase"` (all 754). Case-insensitive; ignored when names are given. |

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
| 1010 | Unknown enzyme set '&lt;value&gt;'; expected "common" or "rebase" (`enzyme_set`) |

## Examples

### Example 1: `AAAGAATTCAAA` with `[EcoRI, BamHI]`

EcoRI cuts once (unique cutter, positions `[3, 3]` for both strands); BamHI does not cut → non-cutter. `totalSites = 1`. The same holds for `[ecori, bamhi]`: `sitesByEnzyme` is keyed `EcoRI` and `nonCutters` is `["BamHI"]`.

### Example 2: `AAAGAATTCAAA` with `enzyme_names: []`

The considered set is the whole common table, so `sitesByEnzyme` is `{ "EcoRI": [3, 3] }` and `nonCutters` lists the other 38 common enzymes.

### Example: all REBASE enzymes

`restriction_map(sequence: "AAGGTCTCAAAAAAAAAAAAAAAGAGACCAAAAAAAAAAA", enzyme_names: [], enzyme_set: "rebase")`
→ `sitesByEnzyme` has 14 enzymes (e.g. `BsaI: [2, 23]` — its forward site and its reverse-strand
`GAGACC` site), `nonCutters` the other **740** (14 + 740 = 754), `totalSites: 41`,
`uniqueCutters: ["SetI"]`. The same call with the default set maps 0 sites and 39 non-cutters.

### Example: unknown enzyme set

`restriction_map(sequence: "GAATTC", enzyme_names: [], enzyme_set: "all")` → error 1010.

## See Also

- [find_restriction_sites](find_restriction_sites.md), [restriction_digest](restriction_digest.md)
