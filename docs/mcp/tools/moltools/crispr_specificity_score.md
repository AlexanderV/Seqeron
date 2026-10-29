# crispr_specificity_score

This library's own heuristic aggregate of a guide RNA's off-targets (NOT a published model).

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `crispr_specificity_score` |
| **Method ID** | `CrisprDesigner.CalculateHeuristicSpecificityScore` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

> **This is this library's own heuristic, not a published model.** Hsu et al. 2013 motivates only the *direction* of the seed weighting (PAM-proximal mismatches are less tolerated). The weights `5` / `2`, the 12-nt seed width, the linear sum and the subtraction from 100 are **unsourced** — no source was found for them and none is claimed. For a publication-grade number use [crispr_mit_specificity_score](crispr_mit_specificity_score.md) (MIT/Hsu 2013) or [crispr_cfd_specificity_score](crispr_cfd_specificity_score.md) (Doench 2016 CFD); both reproduce CRISPOR numerically.

Runs the naïve off-target scan ([find_off_targets](find_off_targets.md)) with `maxMismatches = 4`, then reduces it to a single number: `100` if there are no off-targets, otherwise `max(0, 100 − Σ off-target penalty)`. Each off-target's penalty weights seed-region (PAM-proximal 12 bp) mismatches more heavily (5 vs 2 per mismatch). The on-target (0 mismatches) is not counted. The guide length must equal the system's guide length.

Known consequences of the heuristic: it is **linear** in the mismatch count where the published models are multiplicative, it ignores which bases are involved and the off-target's own PAM, and it **saturates at 0**, so two guides with very different real off-target profiles can share a score of 0. Treat it as coarse triage only.

The C# method is named `CalculateHeuristicSpecificityScore` so that the declaration is structural and not doc-only; the shorter historical `CalculateSpecificityScore` is kept as a backward-compatible alias that forwards to it. The **tool name `crispr_specificity_score` is unchanged** (renaming it would be a breaking catalog change).

## Core Documentation Reference

- Source: [CrisprDesigner.cs#L1086](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs#L1086)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `guide_sequence` | string | Yes | Guide RNA (length must match the system). |
| `genome` | string | Yes | Genome / reference to scan. |
| `system_type` | enum | No | CRISPR system (default `SpCas9`). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `specificity` | number | Heuristic specificity score 0..100 (100 = no off-targets). |

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

- [crispr_mit_specificity_score](crispr_mit_specificity_score.md), [crispr_cfd_specificity_score](crispr_cfd_specificity_score.md), [find_off_targets](find_off_targets.md), [design_guide_rnas](design_guide_rnas.md)
