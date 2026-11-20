# merger

EMBOSS merger: merge two overlapping nucleotide sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `merger` |
| **Method ID** | `SequenceAligner.Merger` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 merger: merges two overlapping nucleotide sequences into one. A global alignment with free end gaps (the program's embAlignPathCalc / embAlignWalkNWMatrix: single-precision path matrices, a gap of length k costs gapOpen + (k - 1) * gapExtend, EDNAFULL scores or the named matrix) is walked column by column: identical bases are copied, an N takes the other sequence's base, and a mismatch or gap takes the base of the sequence with the better quality around the position (2 points per A/C/G/T/U, 1 per other IUPAC code, 0 per N and off the end, windows of 5, 10, 20 on both sides, equal = sequence A) and is listed in conflicts (positions 1-based, '.' = gap). The unaligned head / tail of the sequence that sticks out is appended; mergedSequence is lower case except the bases chosen at conflicts and next to an N (upper case), as the program writes it. Verified against merger 6.6.0 on 450 random overlapping read pairs (point errors, indels, N / IUPAC codes, lower case, six gap-penalty settings): merged sequence, conflict lines and score equal on all. Note the program's own quirk: a tail of the longer sequence that overhangs after the alignment is appended twice when both sequences have unaligned residues at the end - unreachable, since the walk always ends on a matrix border. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.Merger.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.Merger.cs)
- Algorithm doc: [Semi_Global_Alignment.md](../../../algorithms/Alignment/Semi_Global_Alignment.md)
- Oracle: `scripts/alignment/generate_merger_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/merger_oracle.tsv`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequenceA` | string | Yes | First sequence (merger -asequence). |
| `sequenceB` | string | Yes | Second sequence (merger -bsequence). |
| `gapOpen` | number | No | Gap opening penalty (-gapopen, 0..100, default 50): a gap of length 1. (default: 50) |
| `gapExtend` | number | No | Gap extension penalty (-gapextend, 0..10, default 5). (default: 5) |
| `matrix` | string | No | Optional substitution matrix: any built-in matrix by EMBOSS or NCBI name - BLOSUM62, PAM250, EBLOSUM30 ... EBLOSUM90, EPAM10 ... EPAM500 (protein) or EDNAFULL, EDNAMAT, EDNASIMPLE, ENUC.4.4 (IUPAC DNA). When given, match/mismatch are ignored, residues are scored by the published table (any protein or IUPAC symbol of the matrix; case-insensitive) and gapOpen/gapExtend may be fractional (the EMBOSS needle/water models: gapOpen = -(G - E), gapExtend = -E for -gapopen G -gapextend E). The result then also carries the real-valued exactScore. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `mergedSequence` | string | The merged sequence (lower case, conflict bases upper case) |
| `score` | number | Alignment score (single precision) |
| `startA / startB` | integer | Residues of A / B before the aligned region (one of them is 0) |
| `alignedA / alignedB` | string | The aligned region with '.' for gaps |
| `conflicts` | array | `{ positionA, baseA, positionB, baseB, using }` (1-based) |
| `warnings` | array | Warning lines of the report (empty: unreachable with the program's walk) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence is required (null) |
| 1001 | EMBOSS merger cannot merge an empty sequence. |
| 1002 | -gapopen is a finite number between 0.0 and 100.0 in EMBOSS merger. |

## Examples

### Example 1: Dovetail overlap with one mismatch

**User Prompt:**
> Merge two reads that overlap by 20 bases.

**Expected Tool Call:**
```json
{"tool": "merger", "arguments": {"sequenceA": "GATTACAGATTACACCGGTTAACGTACGTTAGC", "sequenceB": "CCGGTTCACGTACGTTAGCAAATTTGGGCCCAAA"}}
```

**Response:**
```json
{"mergedSequence": "gattacagattacaccggttAacgtacgttagcaaatttgggcccaaa", "score": 86, "startA": 14, "startB": 0, "alignedA": "ccggttaacgtacgttagc", "alignedB": "ccggttcacgtacgttagc", "conflicts": [{"positionA": 21, "baseA": "A", "positionB": 7, "baseB": "c", "using": "A"}], "warnings": []}
```

merger prints `# Score: 86.0`, the conflict `21 'A' 7 'c' 'A'` and the merged sequence `gattacagattacaccggttAacgtacgttagcaaatttgggcccaaa`.

## See Also

- [semi_global_align](semi_global_align.md) - ends-free alignment (EndGaps.OverlapSuffix1Prefix2)
