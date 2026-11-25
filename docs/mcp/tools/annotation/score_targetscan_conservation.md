# score_targetscan_conservation

TargetScan 7.0 branch-length score, PCT and conservation call per site of a multi-species 3′UTR alignment.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `score_targetscan_conservation` |
| **Method ID** | `MiRnaAnalyzer.ScoreTargetScanConservation` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

End-to-end TargetScan 7.0 conservation of a miRNA family's sites: site groups (`targetscan_70.pl`) → UTR bin 1–10 from the median per-position branch length (`targetscan_70_BL_bins.pl`) → per-site branch length on the bin tree, PCT with the bundled TargetScan 7.0 sigmoid rows and the Conserved call (`targetscan_70_BL_PCT.pl`; Friedman et al. 2009). `pct` is null for TargetScan's `NA` (seed not one of the 111 PCT families, or ORF-overlapping site). `utrBin` replaces the computed bin (`medianBranchLength` is then null); `referenceSpecies` defaults to 9606; `thresholds` is `TargetScanHuman7` (8mer 1.8, 7mer-m8 2.8, 7mer-A1 3.6; default) or `TargetScanMouse7` (0.6 / 1.8 / 2.5).

## Core Documentation Reference

- Source: [MiRnaAnalyzer.TargetScanConservation.cs#L612](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.TargetScanConservation.cs#L612)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `miRnaSeedOrMature` | string | Yes | 7-nt seed region (miRNA nt 2–8) or a mature miRNA of ≥ 8 nt |
| `alignedUtrs` | object | Yes | Species identifier → aligned 3′UTR (gaps '-'; lower case = ORF overlap) |
| `speciesWithMiRna` | array | No | Species in which the miRNA family is annotated (only these are searched); default all |
| `referenceSpecies` | string | No | Reference species for the UTR bin (default 9606) |
| `utrBin` | integer | No | TargetScan bin 1–10 to use instead of the computed one |
| `thresholds` | string | No | TargetScanHuman7 (default) or TargetScanMouse7 |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `utrBin` | integer | The UTR's bin (1–10) |
| `medianBranchLength` | number | null | Median branch length behind the bin; null when `utrBin` was supplied |
| `groups` | array | Site groups (see find_targetscan_site_groups) |
| `sites` | array | Per site: `site`, `groupNumber`, `speciesWithSite`, `branchLength` (4 decimals), `pct` (number | null), `conserved` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Expected the 7-nt seed region or a mature miRNA of >= 8 nt |
| 1001 | alignedUtrs cannot be null or empty |
| 1001 | referenceSpecies cannot be null or empty |
| 1001 | utrBin must be 1-10 |
| 1001 | Unknown thresholds (valid values listed) |

## Examples

### Example 1: Kit gene LIN28B, let-7/98

Bin 1 (median 1.12468545); the 7mer-m8 sites of the five-species groups have branch length 1.6379 and PCT 0.0058, the 8mer group of three species 1.1469 / 0.0, as in targetscan_70_output.BL_PCT.txt.

**Request:**
```json
{"miRnaSeedOrMature": "GAGGUAG", "alignedUtrs": {"9606": "<aligned UTR of the kit sample>", "10090": "<…>"}, "speciesWithMiRna": ["10090", "10116", "9031", "9606", "9615"]}
```

**Response (abridged):**
```json
{"utrBin": 1, "medianBranchLength": 1.12468545, "sites": [{"groupNumber": 1, "branchLength": 1.6379, "pct": 0.0058, "conserved": false}]}
```

## Performance

- **Time Complexity:** O(S · L) plus tree walks per site
- **Space Complexity:** O(S · L)

## See Also

- [find_targetscan_site_groups](find_targetscan_site_groups.md)
- [score_target_site_context_plus_plus](score_target_site_context_plus_plus.md)
