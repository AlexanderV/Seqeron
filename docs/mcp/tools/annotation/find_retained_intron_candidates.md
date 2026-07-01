# find_retained_intron_candidates

Find short, moderately-scored intron candidates likely to be retained.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_retained_intron_candidates` |
| **Method ID** | `SpliceSitePredictor.FindRetainedIntronCandidates` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Predicts introns (with length bounds 60–500 nt) and returns those that are short (`Length < 500`) and have
only moderate splice scores (`Score < 0.8`). Such weak, short introns are the most likely to be retained in
some transcripts (intron retention).

## Core Documentation Reference

- Source: [SpliceSitePredictor.cs#L905](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs#L905)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | RNA/DNA sequence to analyze |
| `minScore` | number | No | Minimum combined intron score (default 0.5) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `introns` | array | Introns with `Length < 500` and `Score < 0.8` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Short moderate-score intron

A two-exon sequence with an ~83 nt GT-AG intron scoring below 0.8 is returned as a candidate.

### Example 2: No introns present

A sequence with no GU donor yields no candidates:

**Response:**
```json
{ "introns": [] }
```

## Performance

- **Time Complexity:** O(D·A) for donor/acceptor pairing
- **Space Complexity:** O(k) for k candidates

## See Also

- [predict_introns](predict_introns.md) — all intron candidates
- [detect_alternative_splicing](detect_alternative_splicing.md) — alternative splicing events
