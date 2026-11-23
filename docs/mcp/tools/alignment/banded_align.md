# banded_align

Banded pairwise alignment: only the cells of the diagonals diagonalLow .. diagonalHigh (j - i) of the dynamic program are computed, O(length1 x band) time and memory instead of O(length1 x length2) (30 000 x 30 000 residues with affine gaps: about 1 s instead of 98 s in the Debug build). mode global (default), local or semi_global (endGaps as in semi_global_align), linear (default) or affine gaps (gapOpen + k x gapExtend), match / mismatch or a matrix as in global_align. The alignment is the one of the full program (the same rows, score and coordinates) whenever the full alignment lies inside the band; otherwise the best alignment inside the band (never better than the full optimum). GLOBAL: bandWidth w gives the diagonals from min(0, n - m) - w to max(0, n - m) + w; exact = true (the default when no bandWidth is given) doubles the band from initialBandWidth until every path that leaves it is proven to score strictly less (Ukkonen's bound; the result is then the full program's alignment) and reports proven. LOCAL / SEMI-GLOBAL: give diagonalLow and diagonalHigh (a semi-global band must contain the diagonal 0). xDrop prunes cells more than X below the best score seen (BLAST gapped extension; a heuristic, equal to the banded result when X is at least the largest fall of the score along the optimal path). VERIFICATION: equal to GlobalAlign / GlobalAlignAffine / LocalAlign[Affine] / SemiGlobalAlign[Affine] field for field with a band that covers the matrix and with the full alignment inside a narrow band (seeded random pairs, linear / affine, scalar / BLOSUM62).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `banded_align` |
| **Method ID** | `SequenceAligner.BandedAlign` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Banded pairwise alignment: only the cells of the diagonals diagonalLow .. diagonalHigh (j - i) of the dynamic program are computed, O(length1 x band) time and memory instead of O(length1 x length2) (30 000 x 30 000 residues with affine gaps: about 1 s instead of 98 s in the Debug build). mode global (default), local or semi_global (endGaps as in semi_global_align), linear (default) or affine gaps (gapOpen + k x gapExtend), match / mismatch or a matrix as in global_align. The alignment is the one of the full program (the same rows, score and coordinates) whenever the full alignment lies inside the band; otherwise the best alignment inside the band (never better than the full optimum). GLOBAL: bandWidth w gives the diagonals from min(0, n - m) - w to max(0, n - m) + w; exact = true (the default when no bandWidth is given) doubles the band from initialBandWidth until every path that leaves it is proven to score strictly less (Ukkonen's bound; the result is then the full program's alignment) and reports proven. LOCAL / SEMI-GLOBAL: give diagonalLow and diagonalHigh (a semi-global band must contain the diagonal 0). xDrop prunes cells more than X below the best score seen (BLAST gapped extension; a heuristic, equal to the banded result when X is at least the largest fall of the score along the optimal path). VERIFICATION: equal to GlobalAlign / GlobalAlignAffine / LocalAlign[Affine] / SemiGlobalAlign[Affine] field for field with a band that covers the matrix and with the full alignment inside a narrow band (seeded random pairs, linear / affine, scalar / BLOSUM62).

## Core Documentation Reference

- Source: [SequenceAligner.Banded.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.Banded.cs)
- Algorithm: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence. |
| `sequence2` | string | Yes | Second sequence. |
| `mode` | string | No | global (default), local or semi_global. |
| `affine` | boolean | No | true: affine gaps (gapOpen + k x gapExtend); false (default): linear gaps (every position costs gapExtend). |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | number | No | Gap-open score (affine only; <= 0). May be fractional when matrix is given. (default: -2) |
| `gapExtend` | number | No | Gap score per gap position (<= 0). May be fractional when matrix is given. (default: -1) |
| `matrix` | string | No | Optional substitution matrix: any built-in matrix by EMBOSS, NCBI or Biopython name - BLOSUM62, PAM250, EBLOSUM30 ... EBLOSUM90, EPAM10 ... EPAM500, BLASTP (BLOSUM62 with the residues U = selenocysteine, O = pyrrolysine and J = leucine / isoleucine, as blastp scores them), GONNET1992, JONES, DAYHOFF, MDM78, BENNER6 / 22 / 74, JOHNSON, RISLER, GENETIC, FENG, LEVIN, MCLACHLAN, RAO, STR (protein; the fractional ones are scored exactly) or EDNAFULL, EDNAMAT, EDNASIMPLE, ENUC.4.4, BLASTN, MEGABLAST, HOXD70, TRANS (IUPAC DNA); or THE TEXT of a matrix instead of a name: several lines in the NCBI / EMBOSS / Biopython data-file format (a header line of symbols, then one row per symbol: the symbol and its scores; # lines are comments; scores may be decimals). When given, match/mismatch are ignored, residues are scored by the published table (any protein or IUPAC symbol of the matrix; case-insensitive) and gapOpen/gapExtend may be fractional (the EMBOSS needle/water models: gapOpen = -(G - E), gapExtend = -E for -gapopen G -gapextend E). The result then also carries the real-valued exactScore. |
| `bandWidth` | integer | No | global: the half-width w of the band beyond the corridor between the diagonal 0 and the end diagonal n - m (omit for the exact mode). |
| `exact` | boolean | No | global: double the band until the optimum is proven (default true when bandWidth is omitted). |
| `initialBandWidth` | integer | No | exact mode: the first band half-width (default 16). |
| `diagonalLow` | integer | No | local / semi_global: the lowest diagonal j - i of the band. |
| `diagonalHigh` | integer | No | local / semi_global: the highest diagonal j - i of the band. |
| `endGaps` | string | No | Which end gaps are free (default 'Fitting' = the classic semi-global: free leading/trailing sequence2 residues). Comma-separated, case-insensitive names of the flags FreeStart1 (leading sequence1 residues free), FreeEnd1 (trailing sequence1), FreeStart2 (leading sequence2), FreeEnd2 (trailing sequence2) or of a preset: None (= global, all end gaps penalised; result tagged Global), Fitting (FreeStart2+FreeEnd2; sequence1 inside sequence2; parasail sg_dx, Biopython end_insertion_score=0), FittingSequence2 (FreeStart1+FreeEnd1; sequence2 inside sequence1; sg_qx), OverlapSuffix1Prefix2 (FreeStart1+FreeEnd2: a suffix of sequence1 overlaps a prefix of sequence2; sg_qb_de), OverlapSuffix2Prefix1 (FreeEnd1+FreeStart2; sg_qe_db), All (all four free; parasail sg, EMBOSS needle's default -endweight N). Example: 'FreeStart1,FreeEnd2'. startPosition1/endPosition1 and startPosition2/endPosition2 are the consumed intervals of each sequence (the free flanks lie outside them; an empty interval has end = start - 1); a free flank is shown against '-' in the aligned rows. |
| `xDrop` | integer | No | X-drop pruning: cells more than xDrop below the best score so far are dropped (heuristic; omit for none). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignment` | object | The alignment (rows, score, coordinates) |
| `diagonalLow`, `diagonalHigh` | integer | The diagonals j - i of the band that was computed |
| `proven` | boolean | Exact mode: the banded optimum is proven to be the full optimum (null otherwise) |
| `rounds` | integer | Number of bandings (exact mode doubles the band) |
| `cells` | integer | Band cells filled |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence is required |
| 1002 | Unknown mode, a band that does not contain the origin / the end cell, a band wider than 64 000 000 cells per matrix, bandWidth negative, or no alignment inside the band |

## Examples

### Example 1: Example

**Input:** `{"sequence1": "ACGTACGTTAGCTAGCTAGG", "sequence2": "ACGTACGTAGCTAGCTAGG", "mode": "global", "match": 2, "mismatch": -3, "gapExtend": -2}`

**Response:** the global alignment with proven = true (the rows of global_align)

## See Also

- [global_align](global_align.md)
- [global_align_affine](global_align_affine.md)
- [local_align](local_align.md)
- [semi_global_align](semi_global_align.md)
