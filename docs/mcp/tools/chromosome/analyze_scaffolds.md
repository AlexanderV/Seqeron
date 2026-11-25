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

Splits each scaffold at every maximal run of gap symbols (`N`/`n` and, by default (`nonAcgtAsGap = true`, QUAST `correct_seq`), every other non-ACGT symbol) of length ≥ `minGapLength` (QUAST
`qutils.split_by_ns`; default 10 = QUAST `Ns_break_threshold` and the NCBI rule that every run of
≥ 10 N is an assembly gap). Shorter N-runs are ambiguous bases inside a contig (AGP v2.1), so
`contigLength + gapLength = totalLength`. The pieces between gaps become contigs named
`{id}_contig{n}` over 0-based inclusive `[start, end]` indices; each gap is classified by length
(`< 10` Short, `< 100` Medium, `< 1000` Long, else Scaffold — a Seqeron length class, not an AGP gap
type). Reports per scaffold the total length, summed contig length, and summed gap length.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L415](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L415)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `scaffolds` | array | Yes | Scaffold sequences `{ id, sequence }`. Empty list → empty result. |
| `minGapLength` | integer | No | Minimum gap-run length that is a gap and splits contigs (default 10, must be > 0). |
| `nonAcgtAsGap` | boolean | No | QUAST `correct_seq`: treat every non-ACGT symbol (IUPAC codes, U, ...) as N, i.e. a gap symbol (default true); false = only N/n. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[].scaffoldId` | string | Scaffold id. |
| `items[].contigs[]` | array | `{ contigId, start, end }` (end inclusive). |
| `items[].gaps[]` | array | `{ sequenceId, start, end, length, gapType }`. |
| `items[].totalLength` | integer | Scaffold sequence length. |
| `items[].contigLength` | integer | Sum of contig lengths (short N-runs inside contigs included). |
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

Scaffold `AAAAANNNNNTTTTT` with `minGapLength = 10` → one contig `[0, 14]`, no gaps,
`contigLength 15`, `gapLength 0`.

## References

- QUAST `qutils.split_by_ns`, NCBI Gapped Format for Genome Submissions, AGP v2.1 — see `docs/Evidence/ASSEMBLY-GAPS-001-Evidence.md`.

- [GenomeAssemblyAnalyzer.AnalyzeScaffolds](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L415)
