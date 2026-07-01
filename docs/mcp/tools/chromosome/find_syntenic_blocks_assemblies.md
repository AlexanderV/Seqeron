# find_syntenic_blocks_assemblies

Find syntenic blocks between two assemblies (k-mer based).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `find_syntenic_blocks_assemblies` |
| **Method ID** | `GenomeAssemblyAnalyzer.FindSyntenicBlocks` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Anchors shared k-mers between two assemblies and clusters collinear anchors into syntenic blocks of at
least `minBlockSize`, flagging blocks whose orientation is reversed as inverted.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L895](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L895)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `assembly1` | array | Yes | Assembly 1 sequences `{ id, sequence }`. |
| `assembly2` | array | Yes | Assembly 2 sequences `{ id, sequence }`. |
| `minBlockSize` | integer | No | Minimum block size (default 1000, > 0). |
| `kmerSize` | integer | No | K-mer size (default 21, > 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | array | `{ seq1, start1, end1, seq2, start2, end2, isInverted }`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Assembly 1 cannot be null |
| 1002 | Assembly 2 cannot be null |
| 1003 | Minimum block size must be positive |
| 1004 | K-mer size must be positive |

## Example

Two identical 5 kb sequences → one forward block starting at `start1 = 0` with `isInverted = false`.

## References

- [GenomeAssemblyAnalyzer.FindSyntenicBlocks](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L895)
- [compare_assemblies](compare_assemblies.md)
