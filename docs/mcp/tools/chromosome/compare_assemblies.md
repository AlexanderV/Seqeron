# compare_assemblies

Compare two assemblies by exact canonical k-mer sets (Mash model).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `compare_assemblies` |
| **Method ID** | `GenomeAssemblyAnalyzer.CompareAssemblies` |
| **Version** | 1.1.0 |
| **Stability** | Stable |

## Description

Alignment-free whole-assembly comparison by exact canonical k-mer sets (Mash model, Ondov et al. 2016).
All contigs of an assembly are pooled; k-mers are canonical (min of k-mer and reverse complement) and
windows with any non-ACGT symbol are skipped. Reports `distinctKmers1/2`, `sharedKmers`, containment
`|K1 ∩ K2| / |Ki|` (`containment1/2`), the Jaccard index `J = |K1 ∩ K2| / |K1 ∪ K2|` and the Mash distance
`D = −ln(2J/(1+J))/k` (identical sets → 0, no shared k-mer → 1). It reports no aligned bases, identity
percentage or breakpoints (the former `alignedFraction*`, `breakpoints`, `inversions`, `translocations`
(always 0) and `sequenceIdentity` (mean of the fractions, unsourced) fields were removed); use
[find_syntenic_blocks_assemblies](find_syntenic_blocks_assemblies.md) for structure.

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

Fields: `assembly1Name`, `assembly2Name`, `kmerSize`, `distinctKmers1`, `distinctKmers2`, `sharedKmers`,
`containment1`, `containment2`, `jaccardIndex`, `mashDistance`.

## Errors

| Code | Message |
|------|---------|
| 1001 | Assembly 1 cannot be null |
| 1002 | Assembly 2 cannot be null |
| 1003 | K-mer size must be positive |

## Example

Two identical assemblies → `containment1 = containment2 = jaccardIndex = 1.0`, `mashDistance = 0`. `ACGTACGT` vs `ACGTAAAA` (k = 4) → shared 2, J = 1/3, D = ln 2 / 4 = 0.1733.

## References

- [GenomeAssemblyAnalyzer.CompareAssemblies](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L842)
- [find_syntenic_blocks_assemblies](find_syntenic_blocks_assemblies.md)
