# crispr_mit_specificity_score

Published MIT / Hsu 2013 aggregate guide specificity score.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `crispr_mit_specificity_score` |
| **Method ID** | `CrisprDesigner.CalculateMitSpecificityScore` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The **published** MIT / Hsu 2013 aggregate guide specificity score — the "Guide score" of crispr.mit.edu, implemented in CRISPOR as `calcMitGuideScore`. Scans the genome for off-targets with **1..max_mismatches** mismatches ([find_off_targets](find_off_targets.md)), scores each hit with the Hsu position-weight model, and returns `100 / (100 + Σ hitScores) × 100` in `0..100` (100 = no off-targets). The value is **unrounded** — the reference rounds to an integer.

The per-hit model is chosen by system, exactly as the reference's `annotateOfftargets` does: a **SaCas9** PAM (`NNGRRT`/`NNNRRT`) with a 21-nt guide is scored with the SaCas9-specific trained model of **Tycko et al. 2018** (`calcSaHitScore`), every other PAM-after-target system with the Hsu model. A non-SaCas9 guide longer than 20 nt falls back to `calcHitScore`'s own last-20 fudge; a 19-nt guide to its `A`-pad.

Cas12a/Cpf1-class and CasX systems (PAM **before** target) are **refused**: the reference reports no MIT off-target score for them at all (`mitScore = 0.0`, `guideScore = -1`). Prefer this tool over [crispr_specificity_score](crispr_specificity_score.md), which is this library's own unsourced heuristic.

## Core Documentation Reference

- Source: [CrisprDesigner.cs#L1005](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs#L1005)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `guide_sequence` | string | Yes | Guide RNA (length must match the system). |
| `genome` | string | Yes | Genome / reference to scan. |
| `max_mismatches` | integer | No | 0..5 (default 4, the reference's default). |
| `system_type` | enum | No | CRISPR system (default `SpCas9`); must be PAM-after-target. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `specificity` | number | MIT/Hsu aggregate specificity 0..100 (100 = no off-targets). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Guide sequence cannot be null or empty |
| 1002 | Genome cannot be null or empty |
| 1003 | Maximum mismatches must be in the range 0..5 |
| 4001 | Guide length does not match the system's guide length |
| 4002 | The MIT/Hsu 2013 off-target model is not defined for PAM-before-target (Cas12a-class / CasX) systems |

## Examples

### Example 1: one PAM-distal off-target → 50

Guide `GACGCATAAAGATGAGACGC` vs genome `TACGCATAAAGATGAGACGCAGG`: one off-target differing at the PAM-distal position 1 (Hsu weight 0) → single-hit score 100 → `100/(100+100)×100 = 50.0`.

### Example 2: SaCas9 21-nt guide

Guide `ACGTACGTACGTACGTACGTA` vs genome `ACGTACGTACGTACGTACGTCAAGAAT` with `system_type: SaCas9`: one off-target differing at the PAM-proximal position → Tycko 2018 `calcSaHitScore` = 4.2613097364785695 → `95.91285612347565`.

### Example 3: on-target only → `specificity: 100.0`.

## Performance

- **Time Complexity:** O(genome × guide).
- **Space Complexity:** O(number of off-targets).

## See Also

- [crispr_cfd_specificity_score](crispr_cfd_specificity_score.md), [find_off_targets](find_off_targets.md), [crispr_specificity_score](crispr_specificity_score.md)
