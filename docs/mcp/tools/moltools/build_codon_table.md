# build_codon_table

Derive a per-organism codon-usage table from a reference coding sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `build_codon_table` |
| **Method ID** | `CodonOptimizer.CreateCodonTableFromSequence` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Splits the reference into frame-0 codons, converts `T`→`U`, groups codons by the amino acid they encode under the NCBI genetic code `genetic_code` (default 1 = Standard; stop codons of that code form their own family), and computes each codon's relative frequency **within its amino-acid family**. A codon that does not occur in the reference set is counted as **0.5** (Sharp & Li 1987, as implemented by Biopython `CodonAdaptationIndex`), so all 64 codons are present and the relative adaptiveness derived from the table matches Biopython exactly. The resulting table can be fed to the codon-optimizer tools as a custom organism.

## Core Documentation Reference

- Source: [CodonOptimizer.cs#L1794](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs#L1794)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reference_sequence` | string | Yes | Reference coding sequence (DNA or RNA), non-empty. |
| `organism_name` | string | Yes | Organism name to attach to the table (non-blank). |
| `genetic_code` | integer | No | NCBI translation table number (default `1` = Standard, identical to omitting it; e.g. `2` vertebrate mitochondrial, `11` bacterial/plastid). Defines the synonymous families the frequencies are normalised within and the returned codon → amino-acid map. One of 1–6, 9–16, 21–33. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `organismName` | string | The supplied organism name. |
| `codonFrequencies` | object | RNA codon → relative frequency within its amino-acid group (0-1). |
| `codonToAminoAcid` | object | RNA codon → single-letter amino acid (`*` = stop) under `genetic_code` (default Standard). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Reference sequence cannot be null or empty |
| 1002 | Organism name cannot be null or blank |
| 1003 | Unknown genetic code table: N. Supported NCBI tables: 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33. |

## Examples

### Example 1: Two phenylalanine codons

**User Prompt:**
> Build a codon table named MyOrg from "TTTTTTTTC".

**Expected Tool Call:**
```json
{ "tool": "build_codon_table", "arguments": { "reference_sequence": "TTTTTTTTC", "organism_name": "MyOrg" } }
```

Codons `TTT,TTT,TTC` → RNA `UUU,UUU,UUC`; both encode Phe, so `UUU` = 2/3 ≈ 0.6667 and `UUC` = 1/3 ≈ 0.3333.

**Response (abridged):**
```json
{ "organismName": "MyOrg", "codonFrequencies": { "UUU": 0.6667, "UUC": 0.3333 }, "codonToAminoAcid": { "UUU": "F", "UUC": "F" } }
```

## Performance

- **Time Complexity:** O(n) in sequence length.
- **Space Complexity:** O(k) in distinct codons.

## See Also

- [cai_from_organism_table](cai_from_organism_table.md) - Uses a frequency table to score CAI.
- [optimize_codons](optimize_codons.md) - Optimizes against a codon-usage table.
