# inbreeding_from_roh

Estimate the genomic inbreeding coefficient F_ROH from runs of homozygosity.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Population |
| **Tool Name** | `inbreeding_from_roh` |
| **Method ID** | `PopulationGeneticsAnalyzer.CalculateInbreedingFromROH` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes `F_ROH = (Σ L_ROH) / L_AUTO` — the total length of an individual's runs of homozygosity divided by the assayed autosomal genome length (McQuillan et al. 2008). Each segment `[start, end)` contributes length `end − start`. A non-positive `genomeLength` has no defined denominator and returns 0.

## Core Documentation Reference

- Source: [PopulationGeneticsAnalyzer.cs#L1441](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs#L1441)
- Algorithm doc: [Runs_Of_Homozygosity.md](../../../algorithms/Population_Genetics/Runs_Of_Homozygosity.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `rohSegments` | object[] | Yes | Half-open `{start, end}` ROH segments |
| `genomeLength` | integer | Yes | Autosomal genome length L_AUTO (positive) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `inbreedingCoefficient` | number | F_ROH = ΣL_ROH / L_AUTO, in [0, 1] |

## Errors

None (non-positive `genomeLength` returns 0).

## Examples

### Example 1: Two 10 Mb segments over a 100 Mb genome

**User Prompt:**
> F_ROH for two 10 Mb ROH over a 100 Mb genome?

**Expected Tool Call:**
```json
{
  "tool": "inbreeding_from_roh",
  "arguments": {
    "rohSegments": [{ "start": 0, "end": 10000000 }, { "start": 50000000, "end": 60000000 }],
    "genomeLength": 100000000
  }
}
```

**Response:**
```json
{ "inbreedingCoefficient": 0.2 }
```

F_ROH = (10M + 10M)/100M = 0.20.

### Example 2: Whole-genome ROH → F_ROH = 1

**User Prompt:**
> A single ROH covering the whole genome.

**Expected Tool Call:**
```json
{
  "tool": "inbreeding_from_roh",
  "arguments": {
    "rohSegments": [{ "start": 0, "end": 2673768 }],
    "genomeLength": 2673768
  }
}
```

**Response:**
```json
{ "inbreedingCoefficient": 1.0 }
```

## Performance

- **Time Complexity:** O(segments)
- **Space Complexity:** O(1)

## References

- McQuillan, R. et al. (2008). Runs of homozygosity in European populations. *Am. J. Hum. Genet.* 83(3):359–372.

## See Also

- [runs_of_homozygosity](runs_of_homozygosity.md)
