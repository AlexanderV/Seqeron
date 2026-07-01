# crispr_specificity_score

Aggregate a guide RNA's off-targets into a single specificity score.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `crispr_specificity_score` |
| **Method ID** | `CrisprDesigner.CalculateSpecificityScore` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Runs the naïve off-target scan ([find_off_targets](find_off_targets.md)) with `maxMismatches = 4`, then reduces it to a single number: `100` if there are no off-targets, otherwise `max(0, 100 − Σ off-target penalty)`. Each off-target's penalty weights seed-region (PAM-proximal 12 bp) mismatches more heavily (5 vs 2 per mismatch). The on-target (0 mismatches) is not counted. The guide length must equal the system's guide length.

## Core Documentation Reference

- Source: [CrisprDesigner.cs#L378](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs#L378)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `guide_sequence` | string | Yes | Guide RNA (length must match the system). |
| `genome` | string | Yes | Genome / reference to scan. |
| `system_type` | enum | No | CRISPR system (default `SpCas9`). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `specificity` | number | Specificity score 0..100 (100 = no off-targets). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Guide sequence cannot be null or empty |
| 1002 | Genome cannot be null or empty |
| 4001 | Guide length does not match the system's guide length |

## Examples

### Example 1: Unique guide → 100

A 20-nt guide whose only genomic occurrence is its exact on-target has no off-targets, so the score is `100.0`.

## Performance

- **Time Complexity:** O(genome × guide).
- **Space Complexity:** O(number of off-targets).

## See Also

- [find_off_targets](find_off_targets.md), [design_guide_rnas](design_guide_rnas.md)
