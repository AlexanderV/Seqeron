# find_synteny_blocks

Identify collinear synteny blocks from ortholog gene pairs.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `find_synteny_blocks` |
| **Method ID** | `ChromosomeAnalyzer.FindSyntenyBlocks` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Groups ortholog pairs by the pair of chromosomes they connect, sorts them by position in the first
genome, and emits each collinear run of at least `minGenes` genes (respecting `maxGap`, in Mb) as one
block with its strand (`+`/`-`) and gene count.

## Core Documentation Reference

- Source: [ChromosomeAnalyzer.cs#L1348](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L1348)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `orthologPairs` | array | Yes | `{ chr1, start1, end1, gene1, chr2, start2, end2, gene2 }`. |
| `minGenes` | integer | No | Minimum genes per block (default 3, > 0). |
| `maxGap` | integer | No | Max gap between consecutive genes in Mb (default 10). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | array | Blocks with `species1Chromosome`, `species2Chromosome`, coordinates, `strand`, `geneCount`, `sequenceIdentity`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Ortholog pairs cannot be null |
| 1002 | Minimum genes must be positive |

## Example

Five same-order `chr1`↔`chrA` orthologs → one forward block with `geneCount = 5`, `strand = "+"`.

## References

- [ChromosomeAnalyzer.FindSyntenyBlocks](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L1348)
- [detect_rearrangements](detect_rearrangements.md)
