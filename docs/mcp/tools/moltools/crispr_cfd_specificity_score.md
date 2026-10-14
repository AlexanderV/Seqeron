# crispr_cfd_specificity_score

Published Doench 2016 CFD aggregate guide specificity score.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `crispr_cfd_specificity_score` |
| **Method ID** | `CrisprDesigner.CalculateCfdSpecificityScore` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The **published** aggregate CFD (Cutting Frequency Determination, Doench et al. 2016) guide specificity score, exactly as CRISPOR's `annotateOfftargets` computes it. Scans the genome for off-targets with **1..max_mismatches** mismatches ([find_off_targets](find_off_targets.md)), scores each hit with the CFD mismatch/PAM matrices using **that site's own PAM** (read on the protospacer strand, so reverse-strand hits are handled correctly), and aggregates the **raw** per-hit CFD values in `[0, 1]`:

| `aggregation` | Formula | Notes |
|---|---|---|
| `MitStyle` (default) | `100 / (100 + Σ CFD) × 100` | CRISPOR's default, `calcMitGuideScore(sum(cfdScores))`. 100 = no off-targets. |
| `ParkinsonNormalized` | `100 / Σ CFD` | CRISPOR's `doCfdFix` (`fixCfd=1`) variant. Not bounded by 100; **errors when Σ CFD = 0**. |

`ParkinsonNormalized` is written in the reference as `calcCfdGuideScore(sum) = 100 × 100 / sum`, but the same `doCfdFix` flag also rescales every per-hit CFD by ×100 inside `calcCfdScore`; over the raw per-hit values this API takes, the two factors cancel to `100 / Σ raw`.

CFD is defined only for a **20-nt protospacer**, so SaCas9 (21 nt), Cas12a-class and CasX are **refused**. The reference does not score them either: it sets the per-hit sentinel `cfdScore = -1` for a SaCas9 PAM and `cfdScore = 0.0` for a Cpf1 PAM, and the guide-level `guideCfdScore = -1` for a Cpf1 PAM. For SaCas9 its guide-level value is not a sentinel but a vacuous `calcMitGuideScore(sum([])) = 100` — every hit having been dropped by the `if cfdScore != -1` filter — which is exactly the misleading perfect score this tool refuses to return. Prefer this tool over [crispr_specificity_score](crispr_specificity_score.md), which is this library's own unsourced heuristic.

## Core Documentation Reference

- Source: [CrisprDesigner.cs#L1539](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs#L1539)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `guide_sequence` | string | Yes | Guide RNA (must be 20 nt for the supported systems). |
| `genome` | string | Yes | Genome / reference to scan. |
| `max_mismatches` | integer | No | 0..5 (default 4, the reference's default). |
| `system_type` | enum | No | CRISPR system (default `SpCas9`); must be a 20-nt PAM-after-target system. |
| `aggregation` | enum | No | `MitStyle` (default) or `ParkinsonNormalized`. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `specificity` | number | Aggregate CFD specificity (0..100 for `MitStyle`; unbounded for `ParkinsonNormalized`). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Guide sequence cannot be null or empty |
| 1002 | Genome cannot be null or empty |
| 1003 | Maximum mismatches must be in the range 0..5 |
| 4001 | Guide length does not match the system's guide length |
| 4002 | The CFD model is defined only for a 20-nt protospacer (SaCas9 / Cas12a-class / CasX are refused) |
| 4003 | The ParkinsonNormalized CFD aggregate is undefined for a CFD sum of 0 |

## Examples

### Example 1: reverse-strand off-target with a NAG PAM

Guide `GACGCATAAAGATGAGACGC` vs genome `CTTGCGTCTCATCTTTATGCGTA` with `system_type: SpCas9_NAG`: one reverse-strand off-target `TACGCATAAAGATGAGACGC` with PAM `AAG`. Per-hit CFD = `rG:dA@1 (1.0) × PAM "AG" (0.259259259)` → `MitStyle = 99.74141115651946`, `ParkinsonNormalized = 385.7142861`.

### Example 2: on-target only → `specificity: 100.0` (`MitStyle`); `ParkinsonNormalized` errors (Σ CFD = 0).

## Performance

- **Time Complexity:** O(genome × guide).
- **Space Complexity:** O(number of off-targets).

## See Also

- [crispr_mit_specificity_score](crispr_mit_specificity_score.md), [find_off_targets](find_off_targets.md), [crispr_specificity_score](crispr_specificity_score.md)
