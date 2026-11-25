# fickett_score

Fickett (1982) TESTCODE coding statistic of a whole DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `fickett_score` |
| **Method ID** | `GenomeAnnotator.CalculateFickettScore` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

For each base the position parameter `max/(min+1)` of its counts in the three codon phases and its composition
fraction are mapped through Fickett's look-up tables to coding probabilities, weighted and summed. Call: `< 0.74`
NonCoding, `>= 0.95` Coding, otherwise NoOpinion. `tableSet = "fickett1982"` (default) uses the published tables
(EMBOSS `Etcode.dat`), composition over A/C/G/T, exact 4-decimal value; `"cpat"` reproduces CPAT `fickett_value`
(content table without the 0.19 threshold, fractions over the full length, length < 2 → 0), the feature of CPAT's models.

## Core Documentation Reference

- Source: [GenomeAnnotator.cs#L1230](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs#L1230) — `CalculateFickettScore`
- Algorithm: [Coding_Potential_Calculation.md §9](../../../algorithms/Extended_Annotation/Coding_Potential_Calculation.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence |
| `tableSet` | string | No | `fickett1982` (default) or `cpat` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `score` | number | TESTCODE value |
| `call` | string | `NonCoding`, `NoOpinion` or `Coding` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | tableSet must be 'fickett1982' or 'cpat' |

## Example

`TTATATTTTCCACCCAACATCGTGCACATTAGGCTCACCGCGCCCTACGG` → `{ "score": 0.8787, "call": "NoOpinion" }`;
with `tableSet: "cpat"` → `0.8847` (G fraction 0.16 is read one bin early by CPAT's table).

## See Also

- [fickett_scan](fickett_scan.md), [cpat_coding_probability](cpat_coding_probability.md), [coding_potential](coding_potential.md)
