# find_syntenic_blocks_assemblies

Find syntenic blocks between two assemblies (nucmer model: MUMs + mgaps clustering).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `find_syntenic_blocks_assemblies` |
| **Method ID** | `GenomeAssemblyAnalyzer.FindSyntenicBlocks` |
| **Version** | 1.1.0 |
| **Stability** | Stable |

## Description

Finds blocks as MUMmer's `nucmer` does before extension (Delcher et al. 2002; Marçais et al. 2018). Assembly 1
is the reference, assembly 2 the query. For every (reference sequence, query sequence) pair and both query
strands, maximal unique matches of length ≥ `kmerSize` (unique in the reference, nucmer `--mumreference`; both,
`--mum`, with `uniqueInBoth`) are found with the library suffix tree; non-ACGT symbols never match. The anchors
are clustered with the `mgaps` algorithm of MUMmer 4 (repeat filtering, union-find over anchors whose gap ≤
`maxGap` and diagonal difference ≤ max(`diagonalDifference`, `diagonalFactor`·gap), repeated extraction of the
highest-scoring colinear chain); a chain with summed anchor length ≥ `minBlockSize` is a block. Coordinates are
0-based inclusive on the forward strand of each sequence; `isInverted` marks the reverse-complement strand.
No Smith–Waterman extension between anchors is done (a block spans first anchor start to last anchor end).

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L895](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L895)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `assembly1` | array | Yes | Assembly 1 sequences `{ id, sequence }`. |
| `assembly2` | array | Yes | Assembly 2 sequences `{ id, sequence }`. |
| `minBlockSize` | integer | No | Minimum summed anchor length of a block (nucmer `-c`; default 1000, > 0). |
| `kmerSize` | integer | No | Minimum MUM length (nucmer `-l`; default 21, > 0). |
| `maxGap` | integer | No | Maximum gap between clustered anchors (nucmer `-g`; default 90, ≥ 0). |
| `diagonalFactor` | number | No | Diagonal separation factor (nucmer `-d`; default 0.12, ≥ 0). |
| `diagonalDifference` | integer | No | Fixed diagonal difference (nucmer `-D`; default 5, ≥ 0). |
| `uniqueInBoth` | boolean | No | MUMs unique in both assemblies (nucmer `--mum`); default false (`--mumreference`). |

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
| 1005 | Maximum gap / diagonal difference / diagonal factor must be non-negative |

## Example

Two identical 5 kb sequences → one forward block [0, 4999] on both with `isInverted = false`. Query `A + RC(B) + C` (150-bp segments) against `A + B + C` → three blocks, the middle one `isInverted = true`.

## References

- [GenomeAssemblyAnalyzer.FindSyntenicBlocks](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L895)
- [compare_assemblies](compare_assemblies.md)
