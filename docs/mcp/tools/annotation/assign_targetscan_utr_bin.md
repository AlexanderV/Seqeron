# assign_targetscan_utr_bin

TargetScan 7.0 3'UTR conservation bin 1-10 from a multi-species alignment (targetscan_70_BL_bins.pl).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `assign_targetscan_utr_bin` |
| **Method ID** | `MiRnaAnalyzer.AssignTargetScanUtrBin` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

For every alignment column where the reference species has a nucleotide, the branch length (total length of the minimal subtree on the generic tree connecting the species sharing that nucleotide, case-insensitive) is computed; the median over columns selects the bin (largest k whose threshold the median exceeds, 1 if none). Optional `newick` replaces the generic tree (the bundled thresholds then no longer apply); the string is passed unchanged to `PhylogeneticAnalyzer.ParseNewick` (single-line recommended, no normalisation by the tool).

## Core Documentation Reference

- Source: [MiRnaAnalyzer.TargetScanConservation.cs#L459](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.TargetScanConservation.cs#L459)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `alignedUtrs` | object | Yes | Species identifier → aligned 3'UTR ('-' gaps) |
| `referenceSpecies` | string | No | Reference species, default `9606`; must be a leaf of the tree |
| `newick` | string | No | Optional Newick tree; default bundled generic tree |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `bin` | integer | Bin 1-10 |
| `medianBranchLength` | number | null | Median per-position branch length (null when the reference row has no nucleotide) |
| `binThresholds` | number[] | The 10 bin lower bounds (0, 1.21207417, …, 5.78196252) |

## Errors

| Code | Message |
|------|---------|
| 1001 | alignedUtrs cannot be null or empty |
| 1001 | referenceSpecies cannot be null or empty |
| 1001 | The alignment has no row for the reference species |
| 1001 | Invalid Newick tree |

## Examples

### Example 1

Constructed alignment (11 rows, species 99999 not in the tree) gives median 2.55195311, bin 3 (targetscan_70_BL_bins.pl).

**Request:**
```json
{"alignedUtrs": {"9606": "ACGUACGUAC-GUACGUACGU", "9598": "NCG-ACNUA--GNAC-UANGU", "9544": "NCGUNCGUNC-GNACGNACGN", "10090": "NCGUA-GUAC-GUAC-UACGN", "10116": "ncguacnuac-gnacguangu", "9913": "NCGUACG-AC-GUANGUACGU", "9615": "NCGUACGUNC-GUACGNACGU", "9031": "NCGUACGUA--GUACGUANGU", "8364": "NCGUACGUAC-GUACGUACGN", "7897": "NCGUACGUAC--", "99999": "NCGUACGUAC-GNACGUACGU"}}
```

**Response (abridged):**
```json
{"bin": 3, "medianBranchLength": 2.55195311, "binThresholds": "[10 values]"}
```

## Performance

- **Time Complexity:** linear in the input size (alignment columns × species for the UTR bin)
- **Space Complexity:** O(n)

## See Also

- [score_targetscan_conservation](score_targetscan_conservation.md)
- [score_target_site_context_plus_plus](score_target_site_context_plus_plus.md)
