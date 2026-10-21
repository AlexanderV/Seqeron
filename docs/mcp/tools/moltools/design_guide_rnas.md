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

Finds every PAM site of the chosen CRISPR system whose target/cut window falls inside `[region_start, region_end]` (0-based, `region_end` inclusive), scores each candidate guide (GC%, seed-region GC%, polyT terminator, self-complementarity, common restriction sites → 0–100), and returns those scoring at or above `parameters.MinScore`. The composite 0–100 score is CHOPCHOP's own penalty sum shown as `clamp(100 − penalty, 0, 100)`: 1 point per CHOPCHOP self-complementary 4-bp stem (`selfComplementaryStems`, `SCORE['FOLDING']`) plus CHOPCHOP's flat `SCORE['CRISPR_BAD_GC'] = 300` (so score 0) for a GC outside `minGcContent`..`maxGcContent` (default 40–70, edges inclusive). `scoreGcContent` (default `true`) = `false` reproduces CHOPCHOP's own default, which leaves the GC term off. TTTT, seed-region GC and common cloning sites are listed in `issues` only, with no score effect (CRISPOR warns the same way). `minScore` is not from any tool (CHOPCHOP and CRISPOR have no score threshold; their filters are the GC window and `maxSelfComplementaryStems`); for a 20-nt guide any value from 1 to 84 keeps exactly the guides inside the GC window. `selfComplementarityScore` is `selfComplementaryStems` divided by the `length − 4` windows CHOPCHOP examines (range [0, 1]; it replaced the former unsourced pair fraction) and is reported only.

## Core Documentation Reference

- Source: [CrisprDesigner.cs#L188](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs#L188)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence containing the target region, non-empty. |
| `region_start` | integer | Yes | 0-based start (`0 ≤ region_start < sequence.Length`). |
| `region_end` | integer | Yes | 0-based inclusive end (`region_start ≤ region_end < sequence.Length`). |
| `system_type` | enum | No | CRISPR system (default `SpCas9`). |
| `parameters` | object | No | Optional `GuideRnaParameters` (`minGcContent`, `maxGcContent`, `minScore`, `avoidPolyT` (report TTTT), `checkSelfComplementarity`, `scoreGcContent`, `ranking`, `maxSelfComplementaryStems`, `selfComplementarityBackboneRegions`); defaults when null. `maxSelfComplementaryStems` (default null = no filter, CHOPCHOP `filterSelfCompMax = -1`) drops guides with more CHOPCHOP 4-bp self-complementary stems; `selfComplementarityBackboneRegions` (default null = guide-internal stems only) lists sgRNA scaffold regions on the guide's strand that the stems may also pair against (CHOPCHOP `--backbone`). |

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
