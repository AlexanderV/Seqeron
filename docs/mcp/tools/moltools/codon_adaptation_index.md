# codon_adaptation_index

Codon Adaptation Index from a reference RSCU table.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `codon_adaptation_index` |
| **Method ID** | `CodonUsageAnalyzer.CalculateCai` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes CAI (Sharp & Li 1987) for a coding DNA sequence using a caller-supplied **reference RSCU** table (DNA-alphabet codons). Relative adaptiveness `w = RSCU(codon) / max RSCU(synonymous)`; CAI is the geometric mean of `w` over scored codons. Single-codon amino acids (Met, Trp) and stop codons are excluded from the mean; a codon with `w < 0.0001` (absent from the reference) is scored with `w = 0.01` (CodonW `cai_out`; Bulmer 1988) — it lowers CAI rather than being dropped. DNA or RNA input, case-insensitive; non-nucleotide triplets are skipped without shifting the frame.

This differs from [cai_from_organism_table](cai_from_organism_table.md), which takes a codon-usage **frequency** table.

## Core Documentation Reference

- Source: [CodonUsageAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs) `CalculateCai`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Coding DNA sequence (frame 0), non-empty. |
| `reference_rscu` | object | Yes | Reference RSCU table: DNA codon → RSCU value. |
| `genetic_code` | integer | No | NCBI translation table number (default `1` = Standard; e.g. `2` vertebrate mitochondrial, `11` bacterial/plastid). Defines the synonymous families and stops (CodonW `-code`). One of 1–6, 9–16, 21–33. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `cai` | number | Codon Adaptation Index in 0..1. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Reference RSCU table cannot be null |
| 1003 | Unknown genetic code table: N. Supported NCBI tables: 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33. |

## Examples

### Example 1: Non-optimal codon

RSCU `{TTT:1.0, TTC:0.5}` → w[TTC] = 0.5. Sequence `TTC` → CAI = 0.5.

### Example 2: Two codons

Sequence `TTTTTC` → w = {1.0, 0.5}; CAI = √0.5 ≈ 0.7071.

### Example 3: Non-standard genetic code (`genetic_code = 2`)

Reference `{ATG: 1.0, ATA: 0.5}`, sequence `ATAATG`. Under NCBI table 2 Met = {ATA, ATG}, so w[ATA] = 0.5, w[ATG] = 1 and CAI = √0.5 ≈ 0.7071; under the default table 1 ATG (Met) is a single-codon family and is excluded, ATA is Ile with family maximum 0.5, so CAI = 1.0. With the Sharp & Li E. coli reference (`CodonUsageAnalyzer.EColiOptimalCodons`) CodonW 1.4.4 gives `ATAATG` 0.055 (`-code 1`, NCBI 2) vs 0.003 (NCBI 1).

## Performance

- **Time Complexity:** O(n) in sequence length.
- **Space Complexity:** O(1) beyond the reference table.

## See Also

- [cai_from_organism_table](cai_from_organism_table.md) - CAI from a frequency table.
- [rscu](rscu.md) - Compute RSCU for a sequence.
