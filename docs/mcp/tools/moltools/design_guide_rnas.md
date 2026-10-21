# design_guide_rnas

Design and score CRISPR guide-RNA candidates targeting a region.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `design_guide_rnas` |
| **Method ID** | `CrisprDesigner.DesignGuideRnas` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds every PAM site of the chosen CRISPR system whose target/cut window falls inside `[region_start, region_end]` (0-based, `region_end` inclusive), scores each candidate guide (GC%, seed-region GC%, polyT terminator, self-complementarity, common restriction sites → 0–100), and returns those scoring at or above `parameters.MinScore`. The composite 0–100 score is this library's own heuristic with unsourced weights, except that each CHOPCHOP self-complementary 4-bp stem (`selfComplementaryStems`) deducts 1 point, as in CHOPCHOP's `SCORE['FOLDING']`; the legacy `selfComplementarityScore` is an unsourced fraction that is reported only.

## Core Documentation Reference

- Source: [CrisprDesigner.cs#L188](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs#L188)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence containing the target region, non-empty. |
| `region_start` | integer | Yes | 0-based start (`0 ≤ region_start < sequence.Length`). |
| `region_end` | integer | Yes | 0-based inclusive end (`region_start ≤ region_end < sequence.Length`). |
| `system_type` | enum | No | CRISPR system (default `SpCas9`). |
| `parameters` | object | No | Optional `GuideRnaParameters` (`minGcContent`, `maxGcContent`, `minScore`, `avoidPolyT`, `checkSelfComplementarity`, `ranking`, `maxSelfComplementaryStems`, `selfComplementarityBackboneRegions`); defaults when null. `maxSelfComplementaryStems` (default null = no filter, CHOPCHOP `filterSelfCompMax = -1`) drops guides with more CHOPCHOP 4-bp self-complementary stems; `selfComplementarityBackboneRegions` (default null = guide-internal stems only) lists sgRNA scaffold regions on the guide's strand that the stems may also pair against (CHOPCHOP `--backbone`). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `guides` | array | `GuideRnaCandidate` records (sequence, position, strand, GC%, seed GC%, polyT, self-complementarity, score, issues, system). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Region start must be within the sequence |
| 1003 | Region end must satisfy region_start <= region_end < sequence.Length |

## Examples

### Example 1: Single SpCas9 guide with an NGG PAM

Sequence `ACGTACGT…ACGTAGG` (46 nt), region `[20, 45]`, SpCas9. The trailing `AGG` PAM yields one forward-strand guide `ACGTACGTACGTACGTACGT` at position 24 with a perfect score of 100 (50% GC, no polyT).

**Input:** `{ "region_start": 20, "region_end": 45, "system_type": "SpCas9" }`

**Response (abridged):**
```json
{ "guides": [ { "position": 24, "score": 100.0, "isForwardStrand": true, "gcContent": 50.0, "sequence": "ACGTACGTACGTACGTACGT" } ] }
```

### Example 2: Invalid region

`design_guide_rnas("ACGT…", 30, 10)` (end before start) throws `ArgumentException`.

## See Also

- [find_pam_sites](find_pam_sites.md), [evaluate_guide_rna](evaluate_guide_rna.md), [crispr_specificity_score](crispr_specificity_score.md)
