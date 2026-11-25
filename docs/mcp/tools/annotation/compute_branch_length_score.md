# compute_branch_length_score

Friedman 2009 branch-length score (Bls) of a site on a Newick or bundled TargetScan bin tree.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `compute_branch_length_score` |
| **Method ID** | `MiRnaAnalyzer.ComputeBranchLengthScore` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Total branch length of the minimal subtree connecting `speciesWithSite`, rounded to 4 decimals (targetscan_70_BL_PCT.pl getBranchLength); one species, none, or species that are not leaves give 0. Provide exactly one tree: `newick` (passed unchanged to `PhylogeneticAnalyzer.ParseNewick`; single-line recommended, the tool does not normalise whitespace) or `targetScanBin` 1-10 (bundled TargetScan 7.0 `Tree.bin_NN`, 84 species as NCBI taxonomy ids).

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L2153](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L2153)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `speciesWithSite` | string[] | Yes | Species (leaf names) in which the site is conserved |
| `newick` | string | No | Newick tree; exclusive with targetScanBin |
| `targetScanBin` | integer | No | Bundled bin tree 1-10; exclusive with newick |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `branchLengthScore` | number | Bls, 4 decimals |
| `treeSource` | string | `newick` or `targetScanBin N` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Provide exactly one of newick or targetScanBin |
| 1001 | targetScanBin must be 1-10 |
| 1001 | Invalid Newick tree |
| 1001 | speciesWithSite cannot be null |

## Examples

### Example 1

Kit row (CDC2L6, let-7): species 9031 9606 9615 10090 10116 on bin tree 1 give 1.6379.

**Request:**
```json
{"speciesWithSite": ["9031", "9606", "9615", "10090", "10116"], "targetScanBin": 1}
```

**Response (abridged):**
```json
{"branchLengthScore": 1.6379, "treeSource": "targetScanBin 1"}
```

## Performance

- **Time Complexity:** linear in the input size (alignment columns × species for the UTR bin)
- **Space Complexity:** O(n)

## See Also

- [score_targetscan_conservation](score_targetscan_conservation.md)
- [score_target_site_context_plus_plus](score_target_site_context_plus_plus.md)
