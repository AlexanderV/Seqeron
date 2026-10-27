# evaluate_guide_rna

Score a single guide RNA against on-target quality heuristics.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `evaluate_guide_rna` |
| **Method ID** | `CrisprDesigner.EvaluateGuideRna` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Evaluates a single guide RNA and returns a candidate record: overall GC%, seed-region GC% (last 10 nt for Cas9-like, first 10 nt for Cas12a-like), polyT (Pol III terminator) presence, self-complementarity score, a 0..100 quality score, and an issues list. The composite score is CHOPCHOP's own penalty sum shown as `clamp(100 − penalty, 0, 100)`: 1 point per CHOPCHOP self-complementary 4-bp stem (`selfComplementaryStems`, `SCORE['FOLDING']`) plus CHOPCHOP's flat `SCORE['CRISPR_BAD_GC'] = 300` (score 0) for a GC outside `minGcContent`..`maxGcContent` (default 40–70; `scoreGcContent = false` reproduces CHOPCHOP's default of no GC term). TTTT, seed GC and common cloning sites are reported in `issues` only (no score effect, as in CRISPOR). `selfComplementarityScore` is `selfComplementaryStems` divided by the `length − 4` windows CHOPCHOP examines (range [0, 1]; it replaced the former unsourced pair fraction) and is reported only. The guide must contain only A/C/G/T (`U` is read as `T`, case is folded ASCII-only); anything else, IUPAC codes included, is rejected. For ad-hoc evaluation `position` is `-1` and `isForwardStrand` is `true`.

## Core Documentation Reference

- Source: [CrisprDesigner.cs#L325](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs#L325)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `guide_sequence` | string | Yes | Guide RNA sequence; its length must be the system's guide length (SpCas9 20 nt, SaCas9 21, Cas12a/AsCas12a/LbCas12a 23, CasX 20). |
| `system_type` | enum | No | CRISPR system (default `SpCas9`). |
| `parameters` | object | No | Optional guide-design parameters. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `gcContent` / `seedGcContent` | number | Overall and seed-region GC%. |
| `hasPolyT` | boolean | TTTT terminator present. |
| `selfComplementarityScore` / `score` | number | `selfComplementaryStems / (length − 4)` (CHOPCHOP stem count over its own window count, [0, 1]) and the 0..100 score. |
| `issues` | string[] | Quality warnings. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Guide sequence cannot be null or empty |
| 1002 | Expected A/C/G/T (U is read as T) but found '{char}' at index {i}. |
| 1003 | Guide length ({n}) does not match expected length ({L}) for {system}. (SpCas9 20 nt, SaCas9 21, Cas12a/AsCas12a/LbCas12a 23, CasX 20.) |
| 1004 | `parameters` invalid: GC window not 0 <= minGcContent <= maxGcContent <= 100, NaN `minScore`, negative `maxSelfComplementaryStems` (null = no filter; CHOPCHOP's -1 is refused), or all-zero parameters. |

## Examples

### Example 1: `ATGCATGCATGCATGCATGC` → GC 50%, seed GC 60%, `hasPolyT = false`, `position = -1`.

### Example 2: a guide containing `TTTT` → `hasPolyT = true`.

## See Also

- [design_guide_rnas](design_guide_rnas.md), [crispr_specificity_score](crispr_specificity_score.md)
