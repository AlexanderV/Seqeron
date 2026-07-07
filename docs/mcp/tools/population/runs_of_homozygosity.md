# runs_of_homozygosity

Identify runs of homozygosity (ROH) from per-SNP genotype calls.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Population |
| **Tool Name** | `runs_of_homozygosity` |
| **Method ID** | `PopulationGeneticsAnalyzer.FindROH` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans SNPs in ascending position order (input need not be pre-sorted) with the window-free consecutive-runs method (Marras et al. 2015). A candidate run extends while it holds at most `maxHeterozygotes` opposite (heterozygous) genotypes. A run that exceeds the tolerance is closed at its last homozygous SNP and a new run starts at the breaking SNP. A closed run is reported only if it contains at least `minSnps` SNPs and spans at least `minLength` bp (PLINK 1.9 `--homozyg-snp` / `--homozyg-kb`). Genotype 0 (hom-ref) and 2 (hom-alt) are both homozygous; 1 is heterozygous.

## Core Documentation Reference

- Source: [PopulationGeneticsAnalyzer.cs#L1485](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs#L1485)
- Algorithm doc: [Runs_Of_Homozygosity.md](../../../algorithms/Population_Genetics/Runs_Of_Homozygosity.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genotypes` | object[] | Yes | Per-SNP `{position, genotype}` (0/1/2) |
| `minSnps` | integer | No | Minimum SNPs in a run (default 50) |
| `minLength` | integer | No | Minimum run length in bp (default 1,000,000) |
| `maxHeterozygotes` | integer | No | Maximum tolerated heterozygous calls per run (default 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | object[] | Reported runs |
| `items[].start` | integer | Position of the first SNP in the run |
| `items[].end` | integer | Position of the last SNP in the run |
| `items[].snpCount` | integer | Number of SNPs in the run |

## Errors

| Code | Message |
|------|---------|
| 1001 | `minSnps` < 1, or `minLength` / `maxHeterozygotes` negative. |

## Examples

### Example 1: One uninterrupted run of 100 homozygous SNPs

**User Prompt:**
> Find ROH in 100 homozygous SNPs spaced 20 kb apart.

**Expected Tool Call:**
```json
{
  "tool": "runs_of_homozygosity",
  "arguments": {
    "genotypes": [ { "position": 0, "genotype": 0 }, "…", { "position": 1980000, "genotype": 0 } ],
    "minSnps": 50,
    "minLength": 1000000
  }
}
```

**Response:**
```json
{ "items": [{ "start": 0, "end": 1980000, "snpCount": 100 }] }
```

100 SNPs span 1.98 Mb ≥ 1 Mb and ≥ 50 SNPs → one run.

### Example 2: Run shorter than minSnps → no runs

**User Prompt:**
> 20 homozygous SNPs with minSnps = 100.

**Expected Tool Call:**
```json
{
  "tool": "runs_of_homozygosity",
  "arguments": {
    "genotypes": [ { "position": 0, "genotype": 0 }, "…", { "position": 380000, "genotype": 0 } ],
    "minSnps": 100
  }
}
```

**Response:**
```json
{ "items": [] }
```

## Performance

- **Time Complexity:** O(n log n) (positions sorted internally)
- **Space Complexity:** O(n)

## References

- Marras, G. et al. (2015). *Anim. Genet.* 46(2):110–121 (consecutive-runs method).
- Chang, C. C. et al. (2015). *GigaScience* 4:7 (PLINK 1.9 `--homozyg` defaults).

## See Also

- [inbreeding_from_roh](inbreeding_from_roh.md)
