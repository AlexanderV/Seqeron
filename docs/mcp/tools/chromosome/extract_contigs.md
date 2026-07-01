# extract_contigs

Extract contigs (gap-free runs) from scaffolds.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `extract_contigs` |
| **Method ID** | `GenomeAssemblyAnalyzer.ExtractContigs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Splits each scaffold at `N`/`n` runs and emits every gap-free run whose length ≥ `minContigLength` as
a sequence named `{id}_contig{n}`.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L473](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L473)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `scaffolds` | array | Yes | Scaffolds `{ id, sequence }`. |
| `minContigLength` | integer | No | Minimum contig length (default 200, > 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | array | `{ id, sequence }` contigs. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Scaffolds cannot be null |
| 1002 | Minimum contig length must be positive |

## Example

`300×A` + `20×N` + `250×C` with `minContigLength = 200` → `s_contig1` (300 bp), `s_contig2` (250 bp).

## References

- [GenomeAssemblyAnalyzer.ExtractContigs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L473)
- [analyze_scaffolds](analyze_scaffolds.md)
