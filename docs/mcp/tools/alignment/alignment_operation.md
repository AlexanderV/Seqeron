# alignment_operation

One operation of Bio.Align's Alignment object on an alignment given by its block coordinates (see parse_coordinate_alignment): shape (rows, columns), rows (the printed rows, alignment[i]), row (alignment[row, columns] with a slice or a list), column (alignment[rows, col]), letter (alignment[row, col]), slice (alignment[rows], alignment[rows, columns] with a slice 'a:b:c' or a list [c1, c2, ...]: a new alignment), aligned (pairwise: the aligned blocks), indices (the index of the letter in every column, -1 for a gap), inverse_indices (the column of every letter, -1 if not aligned), frequencies (the weight of each letter / gap per column; the 'weight' annotation of a record), substitutions (the letter-pair counts), reverse_complement, sort (by identifier or residues; reverse), add (alignment + alignment2), map (alignmentJson.map(alignment2Json): reads on a transcript on a genome), mapall (alignmentJson = the multiple alignment, alignmentsJson = one alignment per row), from_same_reference (alignmentsJson: pairwise alignments with the same first sequence combined), parse_printed (rows: gapped rows -> sequences and coordinates), text (str(alignment): the Biopython pretty print, with an optional matrix for the ':' marks), equal / compare (alignment vs alignment2). Slices are Python slices ('1:3', ':', '::-1', '-2:'); negative indices count from the end. VERIFICATION: equal to Biopython 1.88 on 572 recorded cases (every operation, 12 index expressions per alignment).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `alignment_operation` |
| **Method ID** | `AlignmentOperations` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

One operation of Bio.Align's Alignment object on an alignment given by its block coordinates (see parse_coordinate_alignment): shape (rows, columns), rows (the printed rows, alignment[i]), row (alignment[row, columns] with a slice or a list), column (alignment[rows, col]), letter (alignment[row, col]), slice (alignment[rows], alignment[rows, columns] with a slice 'a:b:c' or a list [c1, c2, ...]: a new alignment), aligned (pairwise: the aligned blocks), indices (the index of the letter in every column, -1 for a gap), inverse_indices (the column of every letter, -1 if not aligned), frequencies (the weight of each letter / gap per column; the 'weight' annotation of a record), substitutions (the letter-pair counts), reverse_complement, sort (by identifier or residues; reverse), add (alignment + alignment2), map (alignmentJson.map(alignment2Json): reads on a transcript on a genome), mapall (alignmentJson = the multiple alignment, alignmentsJson = one alignment per row), from_same_reference (alignmentsJson: pairwise alignments with the same first sequence combined), parse_printed (rows: gapped rows -> sequences and coordinates), text (str(alignment): the Biopython pretty print, with an optional matrix for the ':' marks), equal / compare (alignment vs alignment2). Slices are Python slices ('1:3', ':', '::-1', '-2:'); negative indices count from the end. VERIFICATION: equal to Biopython 1.88 on 572 recorded cases (every operation, 12 index expressions per alignment).

## Core Documentation Reference

- Source: [AlignmentOperations.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/AlignmentOperations.cs)
- Algorithm: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `operation` | string | Yes | The operation: shape, rows, row, column, letter, slice, aligned, indices, inverse_indices, frequencies, substitutions, reverse_complement, sort, add, map, mapall, from_same_reference, parse_printed, text, equal, compare. |
| `alignmentJson` | string | No | An alignment as JSON: {"sequences": ["ACGT...", ...] (or objects {"id": ..., "length": ..., "parts": {"0": "ACGT..."}}), "coordinates": [[target positions], [query positions], ...]} (one coordinate row per sequence; a position that decreases along the alignment is a sequence on the reverse strand), as parse_coordinate_alignment returns it. Not needed for from_same_reference and parse_printed. |
| `alignment2Json` | string | No | The second alignment (add, map, equal, compare). |
| `alignmentsJson` | array<string> | No | A JSON array of alignments (mapall: one per row of alignmentJson; from_same_reference: the pairwise alignments). |
| `rows` | string | No | Row index (row, letter) or row slice 'a:b:c' (column, slice; default ':'). |
| `columns` | string | No | Column index (column, letter), column slice 'a:b:c' or a JSON list [1, 2] (row, slice). |
| `reverse` | boolean | No | sort: sort descending. text: nothing. (default: false) |
| `matrix` | string | No | Optional substitution matrix (name or text) for text: ':' marks positive pairs. |
| `printedRows` | array<string> | No | parse_printed: the gapped rows. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `operation` | string | The operation |
| `resultJson` | string | The value as JSON (an alignment as {sequences, coordinates, attrs}, a list, a string or a number) |
| `text` | string | The pretty text for the operation text |

## Errors

| Code | Message |
|------|---------|
| 1001 | alignmentJson / alignment2Json / alignmentsJson required for the operation |
| 1002 | Unknown operation, inconsistent steps or unequal step sizes in the alignment |

## Examples

### Example 1: Example

**Input:** `{"operation": "slice", "alignmentJson": "{\"sequences\": [\"TTACGTTT\", \"ACGT\"], \"coordinates\": [[0, 2, 6, 8], [0, 0, 4, 4]]}", "columns": "1:5"}`

**Response:** the alignment of columns 1 to 4: coordinates [[1, 2, 5], [0, 0, 3]]

## See Also

- [pairwise_aligner](pairwise_aligner.md)
- [alignment_counts](alignment_counts.md)
- [parse_coordinate_alignment](parse_coordinate_alignment.md)
