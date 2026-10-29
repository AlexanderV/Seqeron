# find_off_targets

Naïve genome scan for a guide RNA's off-target sites.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `find_off_targets` |
| **Method ID** | `CrisprDesigner.FindOffTargets` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Enumerates every PAM site in the genome for the chosen system and reports those whose target differs from the guide by **1..max_mismatches** mismatches (the exact on-target with 0 mismatches is excluded). Each hit reports its position, sequence, mismatch count and positions, strand, and an `offTargetScore` that is **this library's own heuristic, not a published model** (5 per mismatch inside the PAM-proximal 12-nt seed + 2 per other mismatch; Hsu 2013 motivates only the direction of the weighting, the weights 5/2 and the 12-nt width are unsourced). For published models over the same hits use [crispr_mit_specificity_score](crispr_mit_specificity_score.md) (MIT/Hsu 2013) or [crispr_cfd_specificity_score](crispr_cfd_specificity_score.md) (Doench 2016 CFD). The guide length must equal the system's guide length. Complexity is O(genome × guide); keep genomes ≲ 1 Mb.

## Core Documentation Reference

- Source: [CrisprDesigner.cs#L674](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs#L674)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `guide_sequence` | string | Yes | Guide RNA (length must match the system). |
| `genome` | string | Yes | Genome / reference to scan. |
| `max_mismatches` | integer | No | 0..5 (default 3). |
| `system_type` | enum | No | CRISPR system (default `SpCas9`). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `offTargets` | array | Off-target hits, each `{position, sequence, mismatches, mismatchPositions, isForwardStrand, offTargetScore}`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Guide sequence cannot be null or empty |
| 1002 | Genome cannot be null or empty |
| 1003 | Maximum mismatches must be in the range 0..5 |
| 1004 | Invalid nucleotide '{char}' at position {i}. (the genome must be A/C/G/T; the guide likewise, U read as T) |
| 4001 | Guide length does not match the system's guide length |

## Examples

### Example 1: 1-mismatch off-target

Guide `ATATATATATATATATATAT` vs genome `GTATATATATATATATATATTGG` → one hit, `GTATATATATATATATATAT`, 1 mismatch at position 0.

### Example 2: On-target only → `offTargets: []`.

## See Also

- [crispr_specificity_score](crispr_specificity_score.md), [design_guide_rnas](design_guide_rnas.md)
