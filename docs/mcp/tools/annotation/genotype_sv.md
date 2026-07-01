# genotype_sv

Genotype a structural variant from supporting-read counts.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `genotype_sv` |
| **Method ID** | `StructuralVariantAnalyzer.GenotypeSV` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Assigns a diploid genotype from the alternate-allele read fraction (`altReads / totalReads`):
`altFraction < 0.1` → **0/0** (quality `refReads·3`); `altFraction > 0.9` → **1/1** (quality `altReads·3`);
`0.3 ≤ altFraction ≤ 0.7` → **0/1** (quality `(ref+alt)·2`); otherwise **0/1** (quality `(ref+alt)·1.5`).
Quality is capped at 99. Zero total reads yields the missing genotype **./.** with quality 0.

## Core Documentation Reference

- Source: [StructuralVariantAnalyzer.cs#L1206](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs#L1206)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sv` | object | Yes | Structural variant to genotype |
| `refReads` | integer | Yes | Reference-supporting reads (≥ 0) |
| `altReads` | integer | Yes | Alternate-supporting reads (≥ 0) |
| `totalReads` | integer | Yes | Total reads spanning the locus (≥ 0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `genotype` | string | `0/0`, `0/1`, `1/1`, or `./.` |
| `quality` | number | Genotype quality (capped at 99) |

## Errors

| Code | Message |
|------|---------|
| 1001 | sv cannot be null |
| 1001 | refReads/altReads/totalReads must be non-negative (ArgumentOutOfRangeException) |

## Examples

### Example 1: Heterozygous

10 ref / 10 alt of 20 total (altFraction 0.5) → `0/1`, quality 40.

**Response:**
```json
{ "genotype": "0/1", "quality": 40 }
```

### Example 2: Homozygous reference

20 ref / 0 alt → `0/0`, quality 60:

**Response:**
```json
{ "genotype": "0/0", "quality": 60 }
```

## Performance

- **Time Complexity:** O(1)
- **Space Complexity:** O(1)

## See Also

- [filter_svs](filter_svs.md) — filter SVs before genotyping
- [annotate_svs](annotate_svs.md) — annotate SVs with genes
