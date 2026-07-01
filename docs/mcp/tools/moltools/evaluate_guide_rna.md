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

Evaluates a single guide RNA and returns a candidate record: overall GC%, seed-region GC% (last 10 nt for Cas9-like, first 10 nt for Cas12a-like), polyT (Pol III terminator) presence, self-complementarity score, a 0..100 quality score, and an issues list. For ad-hoc evaluation `position` is `-1` and `isForwardStrand` is `true`.

## Core Documentation Reference

- Source: [CrisprDesigner.cs#L200](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs#L200)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `guide_sequence` | string | Yes | Guide RNA sequence (non-empty). |
| `system_type` | enum | No | CRISPR system (default `SpCas9`). |
| `parameters` | object | No | Optional guide-design parameters. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `gcContent` / `seedGcContent` | number | Overall and seed-region GC%. |
| `hasPolyT` | boolean | TTTT terminator present. |
| `selfComplementarityScore` / `score` | number | Self-complementarity and 0..100 score. |
| `issues` | string[] | Quality warnings. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Guide sequence cannot be null or empty |

## Examples

### Example 1: `ATGCATGCATGCATGCATGC` → GC 50%, seed GC 60%, `hasPolyT = false`, `position = -1`.

### Example 2: a guide containing `TTTT` → `hasPolyT = true`.

## See Also

- [design_guide_rnas](design_guide_rnas.md), [crispr_specificity_score](crispr_specificity_score.md)
