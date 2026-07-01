# find_repetitive_elements

Find tandem and inverted repeats in a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_repetitive_elements` |
| **Method ID** | `GenomeAnnotator.FindRepetitiveElements` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans a DNA sequence for two classes of repeats: **tandem repeats** (head-to-tail adjacent copies of
a primitive motif — each maximal array reported once by its shortest period) and **inverted repeats**
(a sequence followed by its reverse complement, form W·G·Wᴿ, allowing a bracketed gap/loop). Each hit
reports a 0-based inclusive `start`, an exclusive `end`, a `type` of `tandem_repeat` or
`inverted_repeat`, and the repeat `sequence`. Tandem arrays require at least `minCopies` (≥2) copies
and a total span of at least `minRepeatLength`.

## Core Documentation Reference

- Source: [GenomeAnnotator.cs#L875](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs#L875)

## Input Schema

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `dnaSequence` | string | Yes | — | DNA sequence to scan |
| `minRepeatLength` | integer | No | 10 | Minimum repeat length (nt) |
| `minCopies` | integer | No | 2 | Minimum number of tandem copies (≥2) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `repeats` | array | `{ start, end, type, sequence }` per repeat element |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Tandem repeat

`ATTCGATTCGATTCG` (ATTCG × 3), `minRepeatLength=5` → one `tandem_repeat` spanning [0, 15).

**Response:**
```json
{ "repeats": [ { "start": 0, "end": 15, "type": "tandem_repeat", "sequence": "ATTCGATTCGATTCG" } ] }
```

### Example 2: Inverted repeat

`GAATTC` (the EcoRI palindrome) → one `inverted_repeat` spanning [0, 6).

## Performance

- **Time Complexity:** O(n²) inverted-repeat scan (bounded arm length)
- **Space Complexity:** O(k)

## See Also

- [find_promoter_motifs](find_promoter_motifs.md) - Promoter motifs
- [codon_usage](codon_usage.md) - Codon usage counts
