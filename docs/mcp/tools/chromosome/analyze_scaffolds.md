# analyze_scaffolds

Decompose scaffolds into contigs and gaps.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `analyze_scaffolds` |
| **Method ID** | `GenomeAssemblyAnalyzer.AnalyzeScaffolds` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Splits each scaffold at runs of `N`/`n`. Gap-free runs become contigs named `{id}_contig{n}`
over inclusive `[start, end]` indices. An N-run of length ≥ `minGapLength` is recorded as a gap and
classified by length (`< 10` Short, `< 100` Medium, `< 1000` Long, else Scaffold). Reports per
scaffold the total length, summed contig length, and summed gap length.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L415](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L415)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `scaffolds` | array | Yes | Scaffold sequences `{ id, sequence }`. Empty list → empty result. |
| `minGapLength` | integer | No | Minimum N-run length to record as a gap (default 10, must be > 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[].scaffoldId` | string | Scaffold id. |
| `items[].contigs[]` | array | `{ contigId, start, end }` (end inclusive). |
| `items[].gaps[]` | array | `{ sequenceId, start, end, length, gapType }`. |
| `items[].totalLength` | integer | Scaffold sequence length. |
| `items[].contigLength` | integer | Sum of contig lengths. |
| `items[].gapLength` | integer | Sum of gap lengths. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Scaffolds cannot be null |
| 1002 | Minimum gap length must be positive |

## Example

Scaffold `AAAAA` + 10×`N` + `TTTTT` with `minGapLength = 10`:

```json
{
  "items": [
    {
      "scaffoldId": "scaf1",
      "contigs": [
        { "contigId": "scaf1_contig1", "start": 0, "end": 4 },
        { "contigId": "scaf1_contig2", "start": 15, "end": 19 }
      ],
      "gaps": [ { "start": 5, "end": 14, "length": 10, "gapType": "Medium" } ],
      "totalLength": 20,
      "contigLength": 10,
      "gapLength": 10
    }
  ]
}
```

## References

- [GenomeAssemblyAnalyzer.AnalyzeScaffolds](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L415)
