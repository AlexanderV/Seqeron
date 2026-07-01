# haplotype_blocks

Detect haplotype blocks from adjacent-variant linkage disequilibrium.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Population |
| **Tool Name** | `haplotype_blocks` |
| **Method ID** | `PopulationGeneticsAnalyzer.FindHaplotypeBlocks` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Sorts variants by position and walks them left to right, extending the current block while the r² between adjacent variants is ≥ `ldThreshold`. When r² falls below the threshold the current block is closed (only if it holds ≥ 2 variants) and a new one starts. A block's `start`/`end` are the positions of its first and last variants. The `haplotypes` list is present but left empty in this simplified Gabriel et al. (2002) variant.

## Core Documentation Reference

- Source: [PopulationGeneticsAnalyzer.cs#L787](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs#L787)
- Algorithm doc: [Linkage_Disequilibrium.md](../../../algorithms/Population_Genetics/Linkage_Disequilibrium.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `variants` | object[] | Yes | Variants with `variantId`, `position`, and per-individual `genotypes` (0/1/2) |
| `ldThreshold` | number | No | Minimum adjacent-pair r² to extend a block (default 0.7) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | object[] | Detected blocks |
| `items[].start` | integer | Position of the first variant in the block |
| `items[].end` | integer | Position of the last variant in the block |
| `items[].variants` | string[] | Variant ids in the block |
| `items[].haplotypes` | object[] | Haplotype frequencies (empty in this simplified variant) |

## Errors

None (fewer than 2 variants yields no blocks).

## Examples

### Example 1: Three fully-correlated variants form one block

**User Prompt:**
> Find haplotype blocks for three variants with identical genotypes at positions 100, 200, 300.

**Expected Tool Call:**
```json
{
  "tool": "haplotype_blocks",
  "arguments": {
    "variants": [
      { "variantId": "V1", "position": 100, "genotypes": [0, 0, 1, 1, 2, 2] },
      { "variantId": "V2", "position": 200, "genotypes": [0, 0, 1, 1, 2, 2] },
      { "variantId": "V3", "position": 300, "genotypes": [0, 0, 1, 1, 2, 2] }
    ]
  }
}
```

**Response:** one block spanning 100–300 with variants V1, V2, V3 (identical genotypes → r² = 1.0 ≥ 0.7).

### Example 2: Two blocks separated by an LD break

**User Prompt:**
> Detect blocks where V3 is uncorrelated with its neighbours.

**Expected Tool Call:**
```json
{
  "tool": "haplotype_blocks",
  "arguments": {
    "variants": [
      { "variantId": "V1", "position": 100, "genotypes": [0, 0, 0, 1, 1, 1, 2, 2, 2] },
      { "variantId": "V2", "position": 200, "genotypes": [0, 0, 0, 1, 1, 1, 2, 2, 2] },
      { "variantId": "V3", "position": 300, "genotypes": [0, 1, 2, 0, 1, 2, 0, 1, 2] },
      { "variantId": "V4", "position": 400, "genotypes": [0, 0, 0, 1, 1, 1, 2, 2, 2] },
      { "variantId": "V5", "position": 500, "genotypes": [0, 0, 0, 1, 1, 1, 2, 2, 2] }
    ]
  }
}
```

**Response:** two blocks — V1–V2 (100–200) and V4–V5 (400–500); the balanced V3 breaks LD (r² = 0).

## Performance

- **Time Complexity:** O(v · n) — one adjacent-pair r² per consecutive variant
- **Space Complexity:** O(v)

## References

- Gabriel, S. B. et al. (2002). The structure of haplotype blocks in the human genome. *Science* 296:2225–2229.

## See Also

- [linkage_disequilibrium](linkage_disequilibrium.md)
