# find_targetscan_site_groups

TargetScan 7.0 seed-match sites of a miRNA family in a multi-species 3′UTR alignment, grouped across species.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_targetscan_site_groups` |
| **Method ID** | `MiRnaAnalyzer.FindTargetScanSiteGroups` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Port of `targetscan_70.pl`: finds 8mer, 7mer-m8, 7mer-A1 and 6mer seed matches in every species of a multi-species 3′UTR alignment (`alignedUtrs`, species → aligned UTR with `-` gaps; lower-case letters flag ORF overlap; the first 14 nucleotides of every UTR are masked as TargetScan does) and groups overlapping sites across species. Coordinates are 1-based (`msaStart`/`msaEnd` alignment columns, `utrStart`/`utrEnd` ungapped UTR). Each group gives its sites, `groupType` (e.g. `7mer-m8+8mer-1a`), species, and the species per site type used by the branch-length / PCT step. `speciesWithMiRna` restricts the searched species (TargetScan 7: those annotated for the family); omit to search all.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.TargetScanConservation.cs#L140](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.TargetScanConservation.cs#L140)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `miRnaSeedOrMature` | string | Yes | 7-nt seed region (miRNA nt 2–8) or a mature miRNA of ≥ 8 nt |
| `alignedUtrs` | object | Yes | Species identifier → aligned 3′UTR (gaps '-'; lower case = ORF overlap) |
| `speciesWithMiRna` | array | No | Species in which the miRNA family is annotated (only these are searched); default all |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `groups` | array | Groups in TargetScan order (numbered from 1): `groupNumber`, `groupType`, `siteTypes`, `species`, `sites` (`species`, `msaStart`, `msaEnd`, `utrStart`, `utrEnd`, `siteType`, `orfOverlap`), `speciesBySiteType` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Expected the 7-nt seed region or a mature miRNA of >= 8 nt |
| 1001 | Seed region contains a character other than A/C/G/U/T |
| 1001 | alignedUtrs cannot be null or empty |
| 1001 | A species key is empty |
| 1001 | The aligned UTR of species 'X' is null |

## Examples

### Example 1: Kit gene LIN28B, let-7/98 (5 species)

Five groups reproduce targetscan_70_output.txt rows 55–59 of the kit sample (group 2 is 7mer-m8+8mer-1a; the 8mer is aligned in 9031, 9606 and 9615).

**Request:**
```json
{"miRnaSeedOrMature": "GAGGUAG", "alignedUtrs": {"9606": "<aligned UTR of the kit sample>", "10090": "<…>"}, "speciesWithMiRna": ["10090", "10116", "9031", "9606", "9615"]}
```

**Response (abridged):**
```json
{"groups": [{"groupNumber": 2, "groupType": "7mer-m8+8mer-1a"}]}
```

## Performance

- **Time Complexity:** O(S · L) for S species of aligned length L (plus pairwise grouping of the sites)
- **Space Complexity:** O(S · L)

## See Also

- [score_targetscan_conservation](score_targetscan_conservation.md)
- [find_mirna_target_sites](find_mirna_target_sites.md)
