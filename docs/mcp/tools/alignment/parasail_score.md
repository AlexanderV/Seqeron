# parasail_score

The score API of parasail 1.3.4 for one query against one or many references: the optimal score with its end coordinates (endQuery, endReference: 0-based index of the last aligned residue), the statistics of the optimal alignment WITHOUT a traceback (matches, similar = positive pairs, length: parasail nw_stats / sw_stats / sg_stats), the full score tables (table; with statistics also the matches / similar / length tables) and the last row and column of the score matrix (rowcol), for global (function nw), local (sw) and every semi-global variant (sg = all end gaps free, sg_qb / sg_qe / sg_db / sg_de, sg_qx, sg_dx and the combinations sg_qb_db, sg_qb_de, sg_qe_db, sg_qe_de: q = query, d = database = reference; b = beginning, e = end, x = both ends). Gaps are parasail's affine penalties: a gap of length k costs open + (k - 1) * extend (both >= 0; BLAST-style 11 / 1 is open 11, extend 1). The program keeps O(min(length)) cells (the table needs m * n) and the query is prepared once for the whole list of references (profile reuse); parallel runs the references on all cores. Scoring: matrix (name or text; integer scores) or match / mismatch over the letters of the sequences. VERIFICATION: equal to parasail 1.3.4 (nw, sw and the eleven sg variants x plain, _stats, _table, _rowcol, _stats_table, _stats_rowcol) on 520 / 520 recorded pairs.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `parasail_score` |
| **Method ID** | `ParasailAligner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The score API of parasail 1.3.4 for one query against one or many references: the optimal score with its end coordinates (endQuery, endReference: 0-based index of the last aligned residue), the statistics of the optimal alignment WITHOUT a traceback (matches, similar = positive pairs, length: parasail nw_stats / sw_stats / sg_stats), the full score tables (table; with statistics also the matches / similar / length tables) and the last row and column of the score matrix (rowcol), for global (function nw), local (sw) and every semi-global variant (sg = all end gaps free, sg_qb / sg_qe / sg_db / sg_de, sg_qx, sg_dx and the combinations sg_qb_db, sg_qb_de, sg_qe_db, sg_qe_de: q = query, d = database = reference; b = beginning, e = end, x = both ends). Gaps are parasail's affine penalties: a gap of length k costs open + (k - 1) * extend (both >= 0; BLAST-style 11 / 1 is open 11, extend 1). The program keeps O(min(length)) cells (the table needs m * n) and the query is prepared once for the whole list of references (profile reuse); parallel runs the references on all cores. Scoring: matrix (name or text; integer scores) or match / mismatch over the letters of the sequences. VERIFICATION: equal to parasail 1.3.4 (nw, sw and the eleven sg variants x plain, _stats, _table, _rowcol, _stats_table, _stats_rowcol) on 520 / 520 recorded pairs.

## Core Documentation Reference

- Source: [ParasailAligner.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ParasailAligner.cs)
- Algorithm: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `query` | string | Yes | The query (parasail s1). |
| `references` | array<string> | Yes | The references (parasail s2, the database sequences), at least one. |
| `function` | string | No | parasail function: nw (global), sw (local) or a semi-global variant sg, sg_qb, sg_qe, sg_db, sg_de, sg_qx, sg_dx, sg_qb_db, sg_qb_de, sg_qe_db, sg_qe_de (default sw). |
| `open` | integer | No | Penalty of the first gap position (>= 0; default 11). (default: 11) |
| `extend` | integer | No | Penalty of each further gap position (>= 0; default 1). (default: 1) |
| `matrix` | string | No | Optional substitution matrix: any built-in matrix by EMBOSS, NCBI or Biopython name - BLOSUM62, PAM250, EBLOSUM30 ... EBLOSUM90, EPAM10 ... EPAM500, BLASTP (BLOSUM62 with the residues U = selenocysteine, O = pyrrolysine and J = leucine / isoleucine, as blastp scores them), GONNET1992, JONES, DAYHOFF, MDM78, BENNER6 / 22 / 74, JOHNSON, RISLER, GENETIC, FENG, LEVIN, MCLACHLAN, RAO, STR (protein; the fractional ones are scored exactly) or EDNAFULL, EDNAMAT, EDNASIMPLE, ENUC.4.4, BLASTN, MEGABLAST, HOXD70, TRANS (IUPAC DNA); or THE TEXT of a matrix instead of a name: several lines in the NCBI / EMBOSS / Biopython data-file format (a header line of symbols, then one row per symbol: the symbol and its scores; # lines are comments; scores may be decimals). When given, match/mismatch are ignored, residues are scored by the published table (any protein or IUPAC symbol of the matrix; case-insensitive) and gapOpen/gapExtend may be fractional (the EMBOSS needle/water models: gapOpen = -(G - E), gapExtend = -E for -gapopen G -gapextend E). The result then also carries the real-valued exactScore. |
| `match` | integer | No | Score of an identical pair when no matrix is given (default 2). |
| `mismatch` | integer | No | Score of a different pair when no matrix is given (default -3). |
| `statistics` | boolean | No | true: matches, similar and length of the optimal alignment (the *_stats functions). (default: false) |
| `table` | boolean | No | true: the score tables (and, with statistics, the matches / similar / length tables); at most 4 000 000 cells per reference. (default: false) |
| `rowcol` | boolean | No | true: the last row and the last column of the score matrix (the *_rowcol functions). (default: false) |
| `parallel` | boolean | No | true: score the references in parallel (the results keep their order). (default: false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `function` | string | The parasail function that ran |
| `count` | integer | Number of references |
| `results` | array | score, endQuery, endReference, matches, similar, length, the tables and the last row / column per reference |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence is required |
| 1002 | Unknown function, a gap penalty below zero, a letter outside the matrix, an empty reference list or a table above 4 000 000 cells |

## Examples

### Example 1: Example

**Input:** `{"query": "ACGTACGTAC", "references": ["TTACGTACGTAA"], "function": "sw", "open": 5, "extend": 2, "statistics": true}`

**Response:** score 18, endQuery 8, endReference 10, 9 matches, 9 similar, length 9

## See Also

- [pairwise_aligner](pairwise_aligner.md)
- [local_align_affine](local_align_affine.md)
