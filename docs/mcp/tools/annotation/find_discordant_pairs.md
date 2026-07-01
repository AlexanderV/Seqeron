# find_discordant_pairs

Identify discordant read pairs (SV signatures).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_discordant_pairs` |
| **Method ID** | `StructuralVariantAnalyzer.FindDiscordantPairs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Flags read pairs whose mapping is inconsistent with the expected library: mates on different chromosomes
(translocation signature), insert size outside `mean ± 3·sd` (deletion/insertion signature), abnormal
orientation (non-FR — inversion or duplication signature), or insert size above `maxInsertSize`. Concordant
FR pairs with an in-window span are excluded (Medvedev et al. 2009; BreakDancer/DELLY conventions).

## Core Documentation Reference

- Source: [StructuralVariantAnalyzer.cs#L154](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs#L154)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `readPairs` | array | Yes | Read pairs `{ readId, chr1, pos1, strand1, chr2, pos2, strand2, insertSize }` |
| `expectedInsertSize` | integer | No | Expected insert size (default 400) |
| `insertSizeStdDev` | integer | No | Insert-size standard deviation (default 50) |
| `maxInsertSize` | integer | No | Hard maximum insert size (default 10000) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `pairs[]` | object | Discordant read-pair signatures (each with `isDiscordant = true`) |

## Errors

| Code | Message |
|------|---------|
| 1001 | readPairs cannot be null |

## Examples

### Example 1: Multiple SV signatures

Given a concordant FR pair, an interchromosomal pair, a large-span pair and an FF (same-strand) pair, the
last three are flagged:

**Response:**
```json
{ "pairs": [ { "readId": "interchrom" }, { "readId": "largespan" }, { "readId": "inversion" } ] }
```

### Example 2: All concordant

```json
{ "pairs": [] }
```

## Performance

- **Time Complexity:** O(n) for n read pairs
- **Space Complexity:** O(k) for k discordant pairs

## See Also

- [cluster_discordant_pairs](cluster_discordant_pairs.md) — cluster pairs into SV candidates
- [find_split_reads](find_split_reads.md) — split-read SV evidence
