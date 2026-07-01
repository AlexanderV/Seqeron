# compare_assemblies

Compare two assemblies by shared k-mer content.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `compare_assemblies` |
| **Method ID** | `GenomeAssemblyAnalyzer.CompareAssemblies` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Extracts k-mer sets from each assembly and reports the fraction of each assembly's k-mers shared with
the other (`alignedFraction1`, `alignedFraction2`) plus an identity proxy = mean of the two fractions.
Structural counts (`breakpoints`, `inversions`, `translocations`) are reported as 0 by this k-mer
comparison.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L842](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L842)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `assembly1` | array | Yes | Assembly 1 sequences `{ id, sequence }`. |
| `assembly2` | array | Yes | Assembly 2 sequences `{ id, sequence }`. |
| `name1` | string | No | Name for assembly 1 (default Assembly1). |
| `name2` | string | No | Name for assembly 2 (default Assembly2). |
| `kmerSize` | integer | No | K-mer size (default 21, > 0). |

## Output Schema

Fields: `assembly1Name`, `assembly2Name`, `alignedFraction1`, `alignedFraction2`, `breakpoints`,
`inversions`, `translocations`, `sequenceIdentity`.

## Errors

| Code | Message |
|------|---------|
| 1001 | Assembly 1 cannot be null |
| 1002 | Assembly 2 cannot be null |
| 1003 | K-mer size must be positive |

## Example

Two identical assemblies → `alignedFraction1 = alignedFraction2 = sequenceIdentity = 1.0`.

## References

- [GenomeAssemblyAnalyzer.CompareAssemblies](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L842)
- [find_syntenic_blocks_assemblies](find_syntenic_blocks_assemblies.md)
