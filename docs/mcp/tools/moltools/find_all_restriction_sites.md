# find_all_restriction_sites

Find sites for every enzyme of a set (the 39 common enzymes by default, or all 754 REBASE enzymes).

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `find_all_restriction_sites` |
| **Method ID** | `RestrictionAnalyzer.FindAllSites` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Performs a comprehensive restriction scan: every enzyme of `enzyme_set` — `"common"` (default, the 39 common laboratory enzymes) or `"rebase"` (all 754 REBASE enzymes with a defined cut) — is searched against both strands of the sequence. Use this when the caller does not know which enzymes to look for. Coordinates, the two records reported per palindromic site, overlapping-site reporting, the `Linear` cut-feasibility rule and the `Circular` origin-spanning search are all exactly as documented for [find_restriction_sites](find_restriction_sites.md).

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L427](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L427)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to scan (non-empty). |
| `topology` | string | No | `Linear` (default) or `Circular`. |
| `enzyme_set` | string | No | `"common"` (default, the 39 common laboratory enzymes) or `"rebase"` (all 754 REBASE enzymes with a defined cut). Case-insensitive. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sites` | array | All restriction sites found across every enzyme. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1010 | Unknown enzyme set '&lt;value&gt;'; expected "common" or "rebase" (`enzyme_set`) |

## Examples

### Example 1: `AAAGAATTCAAA` → includes the EcoRI site at position 3 (forward cut 4, plus reverse strand).

### Example 2: `AAAAAAAA` → `sites: []` (no recognition sequence matches).

### Example 3: `TTCAAAAAAAGAA` + `topology: Circular` → includes the origin-spanning EcoRI site at position 10 (forward cut 11), which the linear scan does not report.

### Example: the REBASE set

`find_all_restriction_sites(sequence: "AAGGTCTCAAAAAAAAAAAAAAAGAGACCAAAAAAAAAAA", enzyme_set: "rebase")`
→ 42 site records from 14 enzymes (`AbaSI`, `Alw26I`, `BcoDI`, `BmeDI`, `BsaI`, `BsmAI`, `Bso31I`,
`BspTNI`, `BstMAI`, `Eco31I`, `EcoO157SI`, `MspJI`, `SetI`, `YkrI`); BsaI's two records are
`{ position: 2, isForwardStrand: true, cutPosition: 9, recognizedSequence: "GGTCTC" }` and
`{ position: 23, isForwardStrand: false, cutPosition: 22, recognizedSequence: "GGTCTC" }`. The default
(common) set finds **no** site in this sequence. The set of enzymes with a site equals Biopython
1.88 `RestrictionBatch(754).search(...)` (locked in `EnzymeSetToolsTests` on a second sequence, linear
and circular).

### Example: unknown enzyme set

`find_all_restriction_sites(sequence: "GAATTC", enzyme_set: "all")` → error 1010.

## See Also

- [find_restriction_sites](find_restriction_sites.md), [restriction_map](restriction_map.md)
