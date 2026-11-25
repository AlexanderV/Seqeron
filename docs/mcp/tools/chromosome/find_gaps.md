# find_gaps

Find gaps (N-runs) in assembled sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `find_gaps` |
| **Method ID** | `GenomeAssemblyAnalyzer.FindGaps` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans each sequence for maximal runs of gap symbols (`seqtk gap`; QUAST `correct_seq` + `split_by_ns`). A run of length ≥ `minGapLength` is
reported as a gap with 0-based `start`, **inclusive** `end` (`length = end − start + 1`; AGP v2.1
1-based = `[start+1, end+1]`, BED = `[start, end+1)`), its length, and a Seqeron length class
(`< 10` Short, `< 100` Medium, `< 1000` Long, else Scaffold) — not an AGP column-7b gap type, which
records linkage evidence a sequence does not carry. By default (`nonAcgtAsGap = true`) every symbol that is not A/C/G/T — N, IUPAC ambiguity codes, U, `-`, ... — is a gap symbol, as
in QUAST and seqtk; with `nonAcgtAsGap = false` only `N`/`n` is a gap and IUPAC codes are bases (AGP v2.1). The default 1 reports every N-run (the gaps counted by `assembly_statistics`); NCBI counts
runs of ≥ 10 N as assembly gaps, so pass `minGapLength = 10` for that convention.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L341](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L341)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array | Yes | Assembled sequences `{ id, sequence }`. |
| `minGapLength` | integer | No | Minimum gap-run length to report (default 1, > 0). |
| `nonAcgtAsGap` | boolean | No | QUAST `correct_seq`: treat every non-ACGT symbol (IUPAC codes, U, ...) as N, i.e. a gap symbol (default true); false = only N/n. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | array | `{ sequenceId, start, end, length, gapType }`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequences cannot be null |
| 1002 | Minimum gap length must be positive |

## Example

`"AAANNNNNGGG"` → one gap `[3, 7]`, length 5, type `Short`.

## References

- QUAST `qutils.split_by_ns`, seqtk `stk_gap`, NCBI Gapped Format for Genome Submissions, AGP v2.1 — see `docs/Evidence/ASSEMBLY-GAPS-001-Evidence.md`.

- [GenomeAssemblyAnalyzer.FindGaps](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L341)
- [gap_distribution](gap_distribution.md)
