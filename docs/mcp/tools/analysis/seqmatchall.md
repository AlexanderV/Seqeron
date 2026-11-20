# seqmatchall

EMBOSS seqmatchall: all-against-all comparison of a set of sequences by exact word matches.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `seqmatchall` |
| **Method ID** | `GenomicAnalyzer.SeqMatchAll` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 seqmatchall: the all-against-all comparison of a set of sequences by exact word matches. For every pair i < j with length(i) > wordSize and length(j) >= wordSize every maximal exact match of at least wordSize symbols is reported (forward strand only; overlapping matches and repeated positions are all listed). Nucleotide or protein sequences (any letters, compared literally and case-insensitively); the program's word table leaves out words around its unknown symbol (N for a nucleotide sequence, X otherwise) with a scan that can keep a word that holds an N after a skipped one - reproduced exactly. Positions are 0-based (the program prints 1-based starts). Matches are listed per pair, by position (the program lists them in hash order: compare as sets). Verified against seqmatchall 6.6.0 on 330 random sets (2-5 sequences, nucleotide / N / lower case / protein, word sizes 2-12, 700+ pairs); on A/C/G/T input it equals find_common_regions mode MaximalExactMatches (MUMmer -maxmatch). Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program). More than 4 194 304 maximal word matches are refused with an ArgumentException (SequenceAligner.WordMatchMaxMatches).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/GenomicAnalyzer.WordMatches.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.WordMatches.cs)
- Algorithm doc: [Common_Region_Detection.md](../../../algorithms/Sequence_Comparison/Common_Region_Detection.md)
- Oracle: `scripts/analysis/generate_seqmatchall_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Analysis/seqmatchall_oracle.tsv`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | The sequences (at least one; nucleotide or protein letters). |
| `wordSize` | integer | No | Word size (seqmatchall -wordsize, default 4, >= 1). (default: 4) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `matches` | array | Matches `{ sequence1, sequence2, position1, position2, length }`: 0-based sequence indices (sequence1 < sequence2), 0-based starts, match length |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least one sequence is required |
| 1001 | Sequences cannot be null or empty |
| 1002 | wordSize must be at least 1 (seqmatchall -wordsize). |

## Examples

### Example 1: Three sequences

**User Prompt:**
> Find all word matches of size 6 between these sequences.

**Expected Tool Call:**
```json
{"tool": "seqmatchall", "arguments": {"sequences": ["ACGTACGTAC", "GGACGTACGTACTT", "TTTT"], "wordSize": 6}}
```

**Response:**
```json
{"matches": [{"sequence1": 0, "sequence2": 1, "position1": 0, "position2": 2, "length": 10}, {"sequence1": 0, "sequence2": 1, "position1": 0, "position2": 6, "length": 6}, {"sequence1": 0, "sequence2": 1, "position1": 4, "position2": 2, "length": 6}]}
```

seqmatchall -wordsize 6 prints `10 a 1..10 b 3..12`, `6 a 1..6 b 7..12`, `6 a 5..10 b 3..8` (1-based); the third sequence is shorter than the word size.

## See Also

- [find_common_regions](find_common_regions.md) - MUMmer-style common regions of two sequences (mode MaximalExactMatches)
