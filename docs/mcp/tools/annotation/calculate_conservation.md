# calculate_conservation

Compute PhyloP-, PhastCons- and GERP-like conservation scores from multi-species alleles.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `calculate_conservation` |
| **Method ID** | `VariantAnnotator.CalculateConservation` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

For each aligned position, treats the first supplied allele as the reference and measures the
conservation fraction = (species matching the reference) / (total species). From that fraction it
derives three scores:

- **phyloP** = `(fraction - 0.5) * 12` (approx. −6…+6),
- **phastCons** = `fraction` (0…1),
- **gerp** = a rejected-substitution-style score `((conserved - total*0.25) / max(1, total - total*0.25)) * 6`,
  clamped to `[-12.36, 6.18]`,

plus **conservedSpeciesCount**, the raw number of species matching the reference. A position with no
alleles yields all-zero scores.

## Core Documentation Reference

- Source: [VariantAnnotator.cs#L1180](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantAnnotator.cs#L1180)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `positions` | array | Yes | Aligned positions `{ chromosome, position, speciesAlleles }` (one char per species; index 0 is the reference) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `scores` | array | Per-position `{ chromosome, position, phyloP, phastCons, gerp, conservedSpeciesCount }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Positions cannot be null or empty |

## Examples

### Example 1: Fully conserved position

10 identical alleles `AAAAAAAAAA` → fraction 1.0.

**Response:**
```json
{ "scores": [ { "chromosome": "chr1", "position": 100, "phyloP": 6.0, "phastCons": 1.0, "gerp": 6.0, "conservedSpeciesCount": 10 } ] }
```

### Example 2: Half-conserved position

`AATT` (reference `A`, 2 of 4 match) → fraction 0.5.

**Response:**
```json
{ "scores": [ { "chromosome": "chr1", "position": 100, "phyloP": 0.0, "phastCons": 0.5, "gerp": 2.0, "conservedSpeciesCount": 2 } ] }
```

## Performance

- **Time Complexity:** O(p · k) — positions times species per position
- **Space Complexity:** O(p)

## See Also

- [find_conserved_elements](find_conserved_elements.md) - Runs of high PhastCons positions
- [predict_pathogenicity](predict_pathogenicity.md) - Uses conservation among ACMG evidence
