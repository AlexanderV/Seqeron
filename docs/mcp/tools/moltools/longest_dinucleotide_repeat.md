# longest_dinucleotide_repeat

Longest dinucleotide tandem-repeat unit count.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `longest_dinucleotide_repeat` |
| **Method ID** | `PrimerDesigner.FindLongestDinucleotideRepeat` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds the longest tandem repeat of a 2-nt unit (e.g. `ATATAT` = 3 units of `AT`) and returns the number of repeat units, case-insensitive. Sequences shorter than 4 nt return 0.

Any 2-mer counts, including homo-dinucleotides (`AAAA` = 2 units of `AA`); no base is a wildcard.

**Library screen (unsourced heuristic).** No authoritative published definition of this screen, or of the
library default limit `MaxDinucleotideRepeats = 4` (`PrimerDesigner.DefaultParameters`), was found
(audit round 3, A3-8). Primer3 has no dinucleotide-repeat setting (its manual lists only
`PRIMER_MAX_POLY_X`, the mononucleotide run, and repeat-library mispriming), so
`PrimerDesigner.Primer3DefaultParameters` disables the limit. Primer3's sourced alternatives:
[longest_homopolymer](longest_homopolymer.md) (`PRIMER_MAX_POLY_X`) and the `mispriming_library` argument of
[design_primers](design_primers.md) (`PRIMER_MISPRIMING_LIBRARY`).

## Core Documentation Reference

- Source: [PrimerDesigner.cs#L1731](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs#L1731)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence (non-empty). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `repeats` | integer | Longest dinucleotide repeat unit count. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: `ATATAT` → `3` (AT × 3).

### Example 2: `AT` → `0` (shorter than 4 nt).

## See Also

- [longest_homopolymer](longest_homopolymer.md), [design_primers](design_primers.md)
