# annotate_regulatory_elements

Find regulatory regions overlapping a variant.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `annotate_regulatory_elements` |
| **Method ID** | `VariantAnnotator.AnnotateRegulatoryElements` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the subset of candidate regulatory regions that overlap a variant. A region overlaps
when it is on the same chromosome and the variant span `[position, position + ref.Length − 1]`
intersects the region interval `[start, end]`. Overlapping regions are echoed back with their
feature type, cell type, score, and associated transcription factors.

## Core Documentation Reference

- Source: [VariantAnnotator.cs#L1287](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantAnnotator.cs#L1287)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `variant` | object | Yes | Variant `{ chromosome, position, reference, alternate, type }` |
| `regulatoryRegions` | array | Yes | Candidate `{ chromosome, start, end, type, cellType?, score?, transcriptionFactors }` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `annotations` | array | Overlapping regions `{ chromosome, start, end, featureType, cellType?, score?, transcriptionFactors }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | variant cannot be null |
| 1001 | regulatoryRegions cannot be null |

## Examples

### Example 1: Variant overlapping a promoter

**Expected Tool Call:**
```json
{
  "tool": "annotate_regulatory_elements",
  "arguments": {
    "variant": { "chromosome": "chr1", "position": 100, "reference": "A", "alternate": "G", "type": "SNV" },
    "regulatoryRegions": [
      { "chromosome": "chr1", "start": 50, "end": 150, "type": "promoter", "transcriptionFactors": [] },
      { "chromosome": "chr1", "start": 200, "end": 300, "type": "enhancer", "transcriptionFactors": [] }
    ]
  }
}
```

**Response:**
```json
{ "annotations": [ { "chromosome": "chr1", "start": 50, "end": 150, "featureType": "promoter" } ] }
```

### Example 2: No overlap

**Response:**
```json
{ "annotations": [] }
```

## Performance

- **Time Complexity:** O(n) in number of candidate regions
- **Space Complexity:** O(k) in number of overlapping regions

## See Also

- [predict_tf_binding_change](predict_tf_binding_change.md) - TF binding disruption by a SNV
- [annotate_variants](annotate_variants.md) - Full variant consequence annotation
