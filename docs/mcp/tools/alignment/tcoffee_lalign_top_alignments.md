# tcoffee_lalign_top_alignments

T-Coffee lalign_id_pair: the top-N non-intersecting local alignments (Huang-Miller SIM).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `tcoffee_lalign_top_alignments` |
| **Method ID** | `SequenceAligner.TCoffeeLalignTopAlignments` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The 'lalign_id_pair' pair-wise method of T-Coffee 13.41: the topAlignments (t_coffee -lalign_n_top, default 10) best NON-INTERSECTING local alignments of two sequences by Huang & Miller's SIM algorithm, with the program's scoring (DNA: match 10 / mismatch 0 style identity matrix, gaps -10 / -1; protein: BLOSUM62), best first. Each alignment has its two aligned rows, 0-based starts and score. Verified against the program's lalign libraries (see multiple_align_tcoffee).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.TCoffeePrimary.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.TCoffeePrimary.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequenceA` | string | Yes | First sequence (letters; gaps rejected). |
| `sequenceB` | string | Yes | Second sequence. |
| `topAlignments` | integer | No | t_coffee -lalign_n_top: number of alignments to report (>= 1; default 10). (default: 10) |
| `protein` | boolean | No | true = protein scoring, false (default) = DNA. |
| `gapOpen` | integer | No | Method GOP override (gap open score, negative). |
| `gapExtend` | integer | No | Method GEP override (gap extend score, negative). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignments` | array | Alignments `{ rowA, rowB, startA, startB, score }`, best first |
| `count` | integer | Number returned |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence is required |
| 1002 | topAlignments must be at least 1 |

## Examples

### Example 1: Top 2

**User Prompt:**
> Top two lalign alignments.

**Expected Tool Call:**
```json
{"tool": "tcoffee_lalign_top_alignments", "arguments": {"sequenceA": "TTTTCGAACTTC", "sequenceB": "TTTTCAACTC", "topAlignments": 2}}
```

**Response:**
```json
{"alignments": [{"rowA": "TTTTCGAACT", "rowB": "TTTTC-AACT", "startA": 0, "startB": 0, "score": 790}, {"rowA": "TTCGAACTTC", "rowB": "TTTTCAACTC", "startA": 2, "startB": 0, "score": 500}], "count": 2}
```

The four best scores of this pair are 790, 500, 400, 370.

## See Also

- [tcoffee_slow_pair_align](tcoffee_slow_pair_align.md)
- [tcoffee_primary_library](tcoffee_primary_library.md)
