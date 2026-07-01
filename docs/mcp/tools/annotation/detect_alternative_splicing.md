# detect_alternative_splicing

Detect candidate alternative splicing patterns.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `detect_alternative_splicing` |
| **Method ID** | `SpliceSitePredictor.DetectAlternativeSplicing` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Heuristically flags candidate alternative-splicing events from predicted donor/acceptor sites: **ExonSkipping**
(a donor with more than one downstream acceptor >60 nt away), **Alt5SS** (multiple donors clustered in a
50-nt window) and **Alt3SS** (multiple acceptors clustered in a 50-nt window). Each event carries its type,
a representative position, and a human-readable description.

## Core Documentation Reference

- Source: [SpliceSitePredictor.cs#L864](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs#L864)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | RNA/DNA sequence to analyze |
| `minScore` | number | No | Minimum splice-site score (default 0.4) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `events[].type` | string | `ExonSkipping`, `Alt5SS`, or `Alt3SS` |
| `events[].position` | integer | Representative position |
| `events[].description` | string | Human-readable description |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Alternative 5' splice sites

Two clustered GU donors produce an `Alt5SS` event at the first donor position:

**Response:**
```json
{ "events": [ { "type": "Alt5SS", "position": 3 } ] }
```

### Example 2: Single donor

A single donor yields no alternative-splicing event:

**Response:**
```json
{ "events": [] }
```

## Performance

- **Time Complexity:** O(D·A) for D donors and A acceptors
- **Space Complexity:** O(k) for k events

## See Also

- [predict_gene_structure](predict_gene_structure.md) — full exon/intron structure
- [find_retained_intron_candidates](find_retained_intron_candidates.md) — retained-intron candidates
