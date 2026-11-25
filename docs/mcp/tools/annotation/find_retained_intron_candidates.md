# find_retained_intron_candidates

Find retained-intron candidates by sequence-level features (open reading frame; optional length / weak-site / GC filters).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_retained_intron_candidates` |
| **Method ID** | `SpliceSitePredictor.FindRetainedIntronCandidates` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Predicts introns (length 60 and up) and returns those compatible with retention by sequence features:
length a multiple of 3 and no stop codon in some phase of the surrounding frame (the ORF stays open;
Jacob & Smith 2017). Braunschweig et al. 2014 report retained introns are shorter, more G+C-rich and have
weaker splice sites but give no numeric cut-offs, so `maxIntronLength`, `maxIntronScore` and `minGcFraction`
are optional filters, off by default. Retention itself needs transcript (read-coverage) evidence.

## Core Documentation Reference

- Source: [SpliceSitePredictor.cs#L905](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs#L905)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | RNA/DNA sequence to analyze |
| `minScore` | number | No | Minimum combined intron score (default 0.5) |
| `maxIntronLength` | integer | No | Maximum intron length (nt); omitted = 100000 |
| `maxIntronScore` | number | No | Keep only introns with combined score below this (weak sites); omitted = off |
| `minGcFraction` | number | No | Minimum G+C fraction (0-1); omitted = off |
| `requireOpenReadingFrame` | boolean | No | No in-frame stop in some phase (default true) |
| `requireLengthMultipleOf3` | boolean | No | Length divisible by 3 (default true) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `introns` | array | Introns satisfying the criteria |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Frame-compatible intron

A two-exon sequence with an intron of length divisible by 3 and no in-frame stop is returned (set `requireOpenReadingFrame` / `requireLengthMultipleOf3` to false to return every predicted intron).

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
