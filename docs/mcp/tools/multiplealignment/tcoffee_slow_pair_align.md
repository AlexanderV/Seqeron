# tcoffee_slow_pair_align

T-Coffee slow_pair: Myers-Miller global alignment of two sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MultipleAlignment |
| **Tool Name** | `tcoffee_slow_pair_align` |
| **Method ID** | `SequenceAligner.TCoffeeSlowPairAlign` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The 'slow_pair' pair-wise method of T-Coffee 13.41 (the global alignment behind its primary libraries): Myers-Miller linear-space global alignment with the program's matrices and penalties - DNA: dna_idmat (10 / 0), gap open -10, extend -1; protein: BLOSUM62. gapOpen / gapExtend / matrix are the method's own overrides (t_coffee -method slow_pair@EP@GOP@-30@GEP@-3@MATRIX@blosum50mt); matrix names are the T-Coffee matrix names (blosum62mt, blosum50mt, pam250mt, dna_idmat, ...). Returns the two aligned rows. Verified against the program's slow_pair libraries (see multiple_align_tcoffee).

## Core Documentation Reference

- Source: [Seqeron.Genomics.TCoffee/TCoffeeAligner.Primary.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.TCoffee/TCoffeeAligner.Primary.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequenceA` | string | Yes | First sequence (letters; gaps rejected). |
| `sequenceB` | string | Yes | Second sequence. |
| `protein` | boolean | No | true = protein matrices / penalties, false (default) = DNA. |
| `matrix` | string | No | The method's MATRIX override (T-Coffee matrix name); default per type. |
| `gapOpen` | integer | No | The method's GOP override (gap open, an integer score, negative = penalty). |
| `gapExtend` | integer | No | The method's GEP override (gap extend, an integer score, negative = penalty). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `rowA` | string | Aligned first sequence |
| `rowB` | string | Aligned second sequence |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence is required |

## Examples

### Example 1: Gap placed inside

**User Prompt:**
> slow_pair alignment of TTTCAACTC and TTTTCAACTC.

**Expected Tool Call:**
```json
{"tool": "tcoffee_slow_pair_align", "arguments": {"sequenceA": "TTTCAACTC", "sequenceB": "TTTTCAACTC"}}
```

**Response:**
```json
{"rowA": "T-TTCAACTC", "rowB": "TTTTCAACTC"}
```

The routine places the gap inside, not at the start.

## See Also

- [tcoffee_lalign_top_alignments](tcoffee_lalign_top_alignments.md)
- [multiple_align_tcoffee](multiple_align_tcoffee.md)
