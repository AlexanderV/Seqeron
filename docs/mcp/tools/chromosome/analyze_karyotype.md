# analyze_karyotype

Analyze a karyotype from chromosome descriptors.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `analyze_karyotype` |
| **Method ID** | `ChromosomeAnalyzer.AnalyzeKaryotype` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Groups autosomes by base name (stripping a trailing `_N` numeric copy suffix, e.g. `chr1_1` and
`chr1_2` both group under `chr1`), counts copies per group, and flags any group whose copy count
differs from `expectedPloidyLevel` as an aneuploidy using ISCN-style nomenclature
(Nullisomy/Monosomy/Disomy/Trisomy/Tetrasomy/Pentasomy/Polysomy). Sex chromosomes are reported
separately and are not grouped for aneuploidy. Reports total counts, genome size, mean chromosome
length, and the list of abnormalities.

## Core Documentation Reference

- Source: [ChromosomeAnalyzer.cs#L279](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L279)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `chromosomes` | array | Yes | Chromosome descriptors `{ name, length, isSexChromosome }`. Empty list → empty karyotype. |
| `expectedPloidyLevel` | integer | No | Expected copies per autosome (default 2, must be > 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `totalChromosomes` | integer | Number of input chromosomes. |
| `autosomeCount` | integer | Number of non-sex chromosomes. |
| `sexChromosomes` | string[] | Names of sex chromosomes. |
| `totalGenomeSize` | integer | Sum of all lengths (bp). |
| `meanChromosomeLength` | number | Mean length across all chromosomes. |
| `ploidyLevel` | integer | Echoes `expectedPloidyLevel`. |
| `hasAneuploidy` | boolean | True if any autosome group deviates from expected ploidy. |
| `abnormalities` | string[] | e.g. `"Trisomy chr21"`, `"Monosomy chr1"`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Chromosomes cannot be null |
| 1002 | Expected ploidy level must be positive |

## Examples

### Example 1: Normal diploid set

Two autosome pairs (`chr1`, `chr2`) plus `chrX`/`chrY`.

```json
{
  "totalChromosomes": 6,
  "autosomeCount": 4,
  "sexChromosomes": ["chrX", "chrY"],
  "ploidyLevel": 2,
  "hasAneuploidy": false,
  "abnormalities": []
}
```

### Example 2: Trisomy 21

Three copies of `chr21` where diploid is expected.

```json
{
  "autosomeCount": 3,
  "hasAneuploidy": true,
  "abnormalities": ["Trisomy chr21"]
}
```

## References

- [ChromosomeAnalyzer.AnalyzeKaryotype](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L279)
- ISCN aneuploidy nomenclature (Wikipedia: Aneuploidy).
