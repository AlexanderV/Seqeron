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

QUAST "broken scaffolds" (`qutils.split_by_ns` + the `--min-contig` filter): splits each scaffold at
every maximal `N`/`n` run of length ≥ `minGapLength` (default 10 = QUAST `Ns_break_threshold` and the
NCBI assembly-gap rule; 1 cuts at every N). Shorter N-runs stay inside the contig (AGP v2.1). Every
piece of length ≥ `minContigLength` is emitted, case preserved, as `{id}_contig{n}` (n counts the
emitted contigs). The contigs are those of `analyze_scaffolds` with the same `minGapLength`.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L473](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L473)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `scaffolds` | array | Yes | Scaffolds `{ id, sequence }`. |
| `minContigLength` | integer | No | Minimum contig length (default 200, > 0). |
| `minGapLength` | integer | No | Minimum N-run length that splits contigs (default 10, > 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | array | `{ id, sequence }` contigs. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Scaffolds cannot be null |
| 1002 | Minimum contig length must be positive |
| 1003 | Minimum gap length must be positive |

## Example

`300×A` + `20×N` + `250×C` with `minContigLength = 200` → `s_contig1` (300 bp), `s_contig2` (250 bp).

`300×A` + `9×N` + `300×C` (default `minGapLength = 10`) → one 609-bp contig; with `minGapLength = 1`
→ two 300-bp contigs.

## References

- QUAST `qutils.split_by_ns` / `broke_scaffolds`, NCBI Gapped Format for Genome Submissions, AGP v2.1 — see `docs/Evidence/ASSEMBLY-GAPS-001-Evidence.md`.

- [GenomeAssemblyAnalyzer.ExtractContigs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L473)
- [analyze_scaffolds](analyze_scaffolds.md)
