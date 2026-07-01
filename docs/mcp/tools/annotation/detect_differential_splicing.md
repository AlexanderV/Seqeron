# detect_differential_splicing

Detect differential splicing between two conditions.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `detect_differential_splicing` |
| **Method ID** | `TranscriptomeAnalyzer.DetectDifferentialSplicing` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Reports splicing events whose PSI change between two conditions (`deltaPSI = psi2 − psi1`) is at least
`deltaPsiThreshold` in absolute value. Positive changes are labelled **IncreasedInclusion**, negative
changes **IncreasedSkipping**; each event also carries the condition-2 PSI as its inclusion level.

## Core Documentation Reference

- Source: [TranscriptomeAnalyzer.cs#L974](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs#L974)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `splicingData` | array | Yes | Per-event `{ geneId, start, end, psiCondition1, psiCondition2 }` |
| `deltaPsiThreshold` | number | No | Minimum \|PSI2 − PSI1\| (default 0.1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `events[].eventType` | string | IncreasedInclusion / IncreasedSkipping |
| `events[].inclusionLevel` | number | PSI in condition 2 |
| `events[].deltaPsi` | number | PSI2 − PSI1 |

## Errors

| Code | Message |
|------|---------|
| 1001 | splicingData cannot be null |

## Examples

### Example 1: Increased inclusion and skipping

G1 (0.2→0.8) is IncreasedInclusion; G3 (0.8→0.3) is IncreasedSkipping.

### Example 2: Below threshold

A 0.02 PSI change is not reported:

**Response:**
```json
{ "events": [] }
```

## Performance

- **Time Complexity:** O(n) for n events
- **Space Complexity:** O(n)

## See Also

- [find_skipped_exon_events](find_skipped_exon_events.md) — single-sample PSI
