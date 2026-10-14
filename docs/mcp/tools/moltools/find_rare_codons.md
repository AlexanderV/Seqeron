# find_rare_codons

Report codons rarer than a threshold in the target organism.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `find_rare_codons` |
| **Method ID** | `CodonOptimizer.FindRareCodons` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Splits the coding sequence into frame-0 codons (T→U), and reports each codon whose frequency in the target organism's codon-usage table is below `threshold` (default 0.15; strict <, or ≤ with `inclusive_threshold` = true — the Sherlocc "equal or lower than the threshold" comparison, `sherlocc.pl` `$avrg <= $avrg_lim`). Only the 64 unambiguous codons are screened: a triplet containing an ambiguity code (N, R, Y, …) is skipped without shifting the frame, and a trailing partial triplet is ignored; a valid codon absent from an inline table has frequency 0. Each result carries the codon's 0-based nucleotide position, RNA codon, amino acid (translated with the NCBI genetic code `genetic_code`, default 1 = Standard; `*` for stop — UGA is `*` under table 1 and `W` under table 2), and frequency. The selection itself is a per-codon table lookup and does not depend on the code. Rare codons slow translation and are candidates for optimization.

## Core Documentation Reference

- Source: [CodonOptimizer.cs#L1114](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs#L1114)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `coding_sequence` | string | Yes | Coding sequence (DNA or RNA), non-empty. |
| `target_organism` | object | Yes | Preset id (`EColiK12` \| `Yeast` \| `Human`) or inline custom table. |
| `threshold` | number | No | Frequency threshold (default 0.15). |
| `inclusive_threshold` | boolean | No | `false` (default): frequency < threshold; `true`: frequency ≤ threshold (Sherlocc). |
| `genetic_code` | integer | No | NCBI translation table number (default `1` = Standard, identical to omitting it; e.g. `2` vertebrate mitochondrial, `11` bacterial/plastid). Only names the reported amino acid (e.g. UGA is `*` under table 1, `W` under table 2); the rare-codon selection is a per-codon table lookup and does not depend on the code. One of 1–6, 9–16, 21–33. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `rareCodons` | array | Rare codons, each `{position, codon, aminoAcid, frequency}`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Coding sequence cannot be null or empty |
| 1002 | Unknown genetic code table: N. Supported NCBI tables: 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33. |

## Examples

### Example 1: Table `{UUU:0.8, UUC:0.1}`, threshold 0.15

`TTTTTC` → `UUC` at position 3 (frequency 0.1 < 0.15) is the only rare codon.

## See Also

- [optimize_codons](optimize_codons.md), [cai_from_organism_table](cai_from_organism_table.md)
