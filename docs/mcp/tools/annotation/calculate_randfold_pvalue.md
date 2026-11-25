# calculate_randfold_pvalue

randfold randomisation test (Bonnet et al. 2004): P-value of the MFE of a sequence against dinucleotide-shuffled sequences, deterministic for a fixed seed.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `calculate_randfold_pvalue` |
| **Method ID** | `MiRnaAnalyzer.CalculateRandfoldPValue` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

randfold dinucleotide mode (`randfold -d`; Bonnet, Wuyts, Rouzé & Van de Peer 2004, Bioinformatics 20:2911): the input is folded with the canonical Turner 2004 MFE engine, `iterations` Altschul–Erickson dinucleotide shuffles (uShuffle-compatible Eulerian-trail shuffle) are folded, and `P = (1 + #{shuffle ΔG° ≤ input ΔG°}) / (iterations + 1)` (|Δ| ≤ 1e-7 counts as equal, as `randfold.c`). Pre-miRNAs give small P.

**Determinism.** The shuffles are drawn from `System.Random(seed)`: the same `sequence` + `iterations` + `seed` (+ `minLoopSize`) always return the same P-value (default seed 42). The P-value is a Monte-Carlo estimate; vary `seed` to see its spread. `iterations` is capped at 10000 (cost: iterations + 1 MFE folds). `found = false` (P = 1) for an empty sequence.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L3707](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L3707)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Candidate nucleotide sequence (RNA or DNA; T read as U) |
| `iterations` | integer | No | Number of dinucleotide shuffles N (randfold default 999; 0 gives P = 1; max 10000) (default 999) |
| `seed` | integer | No | Seed of the System.Random driving the shuffles (default 42) (default 42) |
| `minLoopSize` | integer | No | Minimum hairpin loop size of the folder (NNDB minimum 3) (default 3) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `found` | boolean | False for an empty sequence |
| `freeEnergy` | number | ΔG° of the input's MFE fold (kcal/mol) |
| `pValue` | number | (1 + #lower-or-equal shuffles)/(N + 1) |
| `iterations` | integer | Shuffles folded |
| `lowerOrEqualCount` | integer | Shuffles with ΔG° ≤ the input's |
| `seed` | integer | Seed used |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | iterations must be in [0, 10000] |
| 1001 | minLoopSize must be >= 0 |

## Examples

### Example 1: hsa-mir-21, 99 shuffles, seed 7

**Request:**
```json
{"sequence": "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA", "iterations": 99, "seed": 7}
```

**Response (abridged):**
```json
{"found": true, "freeEnergy": -35.13, "pValue": 0.01, "iterations": 99, "lowerOrEqualCount": 0, "seed": 7}
```

## Performance

- **Time Complexity:** O(iterations · n³) (one MFE fold per shuffle)

## See Also

- [classify_pre_mirna](classify_pre_mirna.md)
- [find_pre_mirna_hairpins](find_pre_mirna_hairpins.md)
