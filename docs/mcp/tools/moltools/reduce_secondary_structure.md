# reduce_secondary_structure

Reduce mRNA secondary structure by synonymous codon swaps.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `reduce_secondary_structure` |
| **Method ID** | `CodonOptimizer.ReduceSecondaryStructure` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Greedily swaps codons for synonymous alternatives that lower the structure score of a sliding window (default 40 nt), reducing mRNA secondary structure while preserving the protein. Two scoring models:

- `PairCount` (default) — fraction of position pairs (i, j), j ≥ i + 4, that can form a canonical pair, using the canonical `RnaSecondaryStructure.CanPair` rule (Watson-Crick plus the G·U wobble, ViennaRNA pair set). A pairing-potential heuristic: no stacking, loop penalties or pair nesting. A window is rewritten when the fraction exceeds 0.5.
- `Turner2004Mfe` — the window's minimum free energy from the canonical Turner 2004 nearest-neighbour engine (`RnaSecondaryStructure.CalculateMinimumFreeEnergy`, the energy model of ViennaRNA RNAfold), maximised (made less negative). A window is rewritten when its MFE is below 0 kcal/mol and a swap is accepted only when it raises the window MFE. Window folds are memoised and capped, but this mode is ~150× slower than `PairCount` (~9 s per kilobase vs ~0.06 s, measured on a 999-nt CDS).

Output is upper-case RNA trimmed to whole codons; a sequence shorter than `window_size` keeps all its codons.

## Core Documentation Reference

- Source: [CodonOptimizer.cs#L990](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs#L990)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `coding_sequence` | string | Yes | Coding sequence (DNA or RNA), non-empty. |
| `target_organism` | object | Yes | Preset id or inline custom table. |
| `window_size` | integer | No | Sliding-window size in nt (default 40, positive). |
| `scoring` | string | No | `PairCount` (default) or `Turner2004Mfe`. |
| `genetic_code` | integer | No | NCBI translation table number (default `1` = Standard, identical to omitting it; e.g. `2` vertebrate mitochondrial, `11` bacterial/plastid). Synonymous swaps stay inside this code's families. One of 1–6, 9–16, 21–33. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `optimizedSequence` | string | Sequence with reduced local self-complementarity. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Coding sequence cannot be null or empty |
| 1002 | Window size must be positive |
| 1003 | Unknown secondary-structure scoring |
| 1004 | Unknown genetic code table: N. Supported NCBI tables: 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33. |

## Examples

### Example 1: Below window → `ATGATG` → `AUGAUG` (codons unchanged, RNA-normalised).

### Example 2: A 45-nt sequence (> 40-nt window) is processed with synonymous swaps that preserve its length and protein.

### Example 3: `scoring = Turner2004Mfe` on a 99-nt CDS (`MKLADLAQLAFREDLRQALDRFMAQLDKALAE*`) rewrites 19 codons and raises the summed 45-nt window MFE from -252.88 to -96.45 kcal/mol (worst window -17.00 -> -7.18), protein unchanged.

## See Also

- [optimize_codons](optimize_codons.md) (strategy `MinimizeSecondary` runs this pass)
- [minimum_free_energy](../analysis/minimum_free_energy.md)
