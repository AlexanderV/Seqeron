# tcoffee_primary_library

The primary library T-Coffee builds for a set of sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MultipleAlignment |
| **Tool Name** | `tcoffee_primary_library` |
| **Method ID** | `Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeePrimaryLibraryEntries` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The primary library the T-Coffee 13.41 program builds for a set of sequences with the given pair-wise methods (t_coffee -seq s.fa -method m1,m2 -lib_only -out_lib x.tc_lib): residue pairs (sequence index of the first smaller than that of the second; 0-based positions) with weight 10 x percent identity of the alignment they come from and the number of times the pair was added. methods in command-line order: 'proba_pair' (default; ProbCons pair-HMM posteriors), 'slow_pair', 'lalign_id_pair', 'clustalw_pair', 'biphasic_pair' (two affine gap models), 'lalign_rs_s_pair' (SIM on every pair including each sequence against itself), 'promo_pair' / 'blastr_pair' (slow_pair on the dinucleotide-extended alphabet), 'fast_pair' (FASTA heuristic of util_dp_fasta_nw.c: banded Gotoh DP on the best k-tuple diagonals), 'exon_pair' / 'exon2_pair' / 'exon3_pair' (alignment of the genomic templates given in genomicTemplates, expressed in the residues of the sequences); the program refuses 'clean_slow_pair', 'hash_pair', 'subop1_pair', 'subop2_pair', 'cfast_pair' (TCoffeeAbortedException); with two or more methods the program's own merge loses and duplicates entries - the result is what it writes. Verified on 800+ committed libraries (see multiple_align_tcoffee). Entries are ordered by (seqA, seqB, posA, posB).

## Core Documentation Reference

- Source: [Seqeron.Genomics.TCoffee/TCoffeeAligner.Library.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.TCoffee/TCoffeeAligner.Library.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Sequences (at least one, none empty; letters only). |
| `methods` | array<string> | No | Pair-wise methods in command-line order ('proba_pair' (default), 'slow_pair', 'lalign_id_pair', 'clustalw_pair', 'biphasic_pair', 'lalign_rs_s_pair', 'promo_pair', 'blastr_pair', 'fast_pair', 'exon_pair', 'exon2_pair', 'exon3_pair'). |
| `lalignTopAlignments` | integer | No | t_coffee -lalign_n_top (default 10). |
| `protein` | boolean | No | true = protein, false = DNA, omitted = decided from the letters as the program does. |
| `genomicTemplates` | array<string> | No | The genomic (_G_) template of every sequence for exon_pair / exon2_pair / exon3_pair (one entry per sequence, null / empty = none). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `entries` | array | Library entries `{ seqA, posA, seqB, posB, weight, count }` (0-based) |
| `count` | integer | Number of entries |

## Errors

| Code | Message |
|------|---------|
| 1001 | sequences must contain at least one sequence |
| 1001 | Unknown pair-wise method |

## Examples

### Example 1: slow_pair library of two sequences

**User Prompt:**
> Primary library of two sequences with slow_pair.

**Expected Tool Call:**
```json
{"tool": "tcoffee_primary_library", "arguments": {"sequences": ["TTTTCAACTC", "TTTCAACTC"], "methods": ["slow_pair"]}}
```

**Response:**
```json
{"entries": [{"seqA": 0, "posA": 1, "seqB": 1, "posB": 0, "weight": 1000, "count": 1}, {"seqA": 0, "posA": 2, "seqB": 1, "posB": 1, "weight": 1000, "count": 1}, {"seqA": 0, "posA": 3, "seqB": 1, "posB": 2, "weight": 1000, "count": 1}, {"seqA": 0, "posA": 4, "seqB": 1, "posB": 3, "weight": 1000, "count": 1}, {"seqA": 0, "posA": 5, "seqB": 1, "posB": 4, "weight": 1000, "count": 1}, {"seqA": 0, "posA": 6, "seqB": 1, "posB": 5, "weight": 1000, "count": 1}, {"seqA": 0, "posA": 7, "seqB": 1, "posB": 6, "weight": 1000, "count": 1}, {"seqA": 0, "posA": 8, "seqB": 1, "posB": 7, "weight": 1000, "count": 1}, {"seqA": 0, "posA": 9, "seqB": 1, "posB": 8, "weight": 1000, "count": 1}], "count": 9}
```

## See Also

- [multiple_align_tcoffee](multiple_align_tcoffee.md)
- [tcoffee_guide_tree](tcoffee_guide_tree.md)
