# multiple_seq_alignment

The legacy Bio.Align.MultipleSeqAlignment of Biopython 1.88: rows / columns / slices / add / sort / substitutions / alignment / append / extend / delete and the column annotations on an alignment of records.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_seq_alignment` |
| **Method ID** | `MultipleSeqAlignment` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

One operation of the legacy Bio.Align.MultipleSeqAlignment of Biopython 1.88 on an alignment of records (rows of one length with identifier, name, description, annotations, per-letter annotations, cross-references; alignment annotations and per-column annotations): summary (str(alignment): at most 20 rows, 50 letters), length (the number of rows), alignment_length, row (alignment[i], negative from the end), rows (alignment[a:b:c]: a sub-alignment; the column annotations are kept when all rows are), row_columns (alignment[i, a:b:c]: the part of a row), cell (alignment[i, j]), column (alignment[a:b:c, j]: the letters of a column), slice (alignment[rows, columns]: a sub-alignment with sliced column annotations), add (alignment + otherJson: columns appended pairwise, common annotations kept, per-column annotations concatenated), sort (by id, length, sequence, description or gaps; reverse), substitutions (the weighted letter-pair counts: the 'weight' annotations multiply, gaps are skipped, the matrix is symmetric with the counts shared between both sides), alignment (the modern Alignment object: gap-free sequences and block coordinates; Biopython 1.88 empties the letter annotations of the source rows in this call, so the returned rows keep theirs and the source is not returned), append / extend (a record / records of the alignment length), delete_row / delete_rows, column_annotations (the restricted dictionary: set_column_annotation adds one entry of the length the dictionary was created with). Slices are Python slices ('1:3', ':', '::-1', '-2:'). Errors are Biopython's (Sequences must all be the same length; list index out of range; index out of range; the per-column and per-letter annotation messages). Feature objects of records and file output through AlignIO are not part of the tool (use write_alignment for the text formats). VERIFICATION: results and error messages equal to Biopython 1.88 on 7 030 operations of 1 200 recorded alignments (scripts/alignment/generate_multiple_seq_alignment_oracle.py).

## Core Documentation Reference

- Source: [MultipleSeqAlignment.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MultipleSeqAlignment.cs)
- Algorithm: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_multiple_seq_alignment_oracle.py](../../../../scripts/alignment/generate_multiple_seq_alignment_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `operation` | string | Yes | The operation: summary, length, alignment_length, row, rows, row_columns, cell, column, slice, add, sort, substitutions, alignment, append, extend, delete_row, delete_rows, column_annotations, set_column_annotation. |
| `recordsJson` | string | Yes | The rows as a JSON array of records {"id": ..., "name": ..., "description": ..., "seq": "ACGT-", "annotations": {...}, "letter_annotations": {"key": "string or array, one entry per letter"}, "dbxrefs": [...]}; only seq is required (id / name / description default to <unknown id> / <unknown name> / <unknown description>). |
| `annotationsJson` | string | No | Annotations of the whole alignment (JSON object). |
| `columnAnnotationsJson` | string | No | Per-column annotations (JSON object: key -> string or array with one entry per column). |
| `row` | integer | No | Row index (row, row_columns, cell, delete_row); negative counts from the end. (default: 0) |
| `rows` | string | No | Row slice 'a:b:c' (rows, column, slice, delete_rows; default ':'). |
| `columns` | string | No | Column slice 'a:b:c' (row_columns, slice; default ':'). |
| `column` | integer | No | Column index (cell, column); negative counts from the end. (default: 0) |
| `otherJson` | string | No | add: the second alignment as a JSON object {"records": [...], "annotations": {...}, "column_annotations": {...}}. |
| `newRecordsJson` | string | No | append: one record (JSON object); extend: a JSON array of records. |
| `sortKey` | string | No | sort key: id (default), length, sequence, description or gaps. |
| `reverse` | boolean | No | sort: descending. (default: false) |
| `annotationKey` | string | No | set_column_annotation: the key. |
| `annotationValueJson` | string | No | set_column_annotation: the value (JSON: a string or an array). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `operation` | string | The operation |
| `resultJson` | string | The value as JSON: an alignment {records, annotations, column_annotations}, a record, a string, a number or a matrix {letters, values} |
| `text` | string | str(alignment) when the value is an alignment |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequences must all be the same length |
| 1002 | list index out of range / index out of range (a row or a column that does not exist) |
| 1003 | Any per-letter annotation should be a Python sequence (list, tuple or string) of the same length as the biological sequence, here N. |
| 1004 | When adding two alignments they must have the same length (i.e. same number of rows) |

## Examples

### Example 1: A sub-alignment of rows 1-4 and columns 3-5

**User Prompt:**
> Take rows 1 to 4 and columns 3 to 5 of this alignment.

**Expected Tool Call:**
```json
{"tool": "multiple_seq_alignment", "arguments": {"operation": "slice", "recordsJson": "[{\"id\": \"Alpha\", \"seq\": \"AAAACGT\"}, {\"id\": \"Beta\", \"seq\": \"AAA-CGT\"}, {\"id\": \"Gamma\", \"seq\": \"AAAAGGT\"}, {\"id\": \"Delta\", \"seq\": \"AAAACGT\"}, {\"id\": \"Epsilon\", \"seq\": \"AAA-GGT\"}]", "rows": "1:5", "columns": "3:6"}}
```

**Response:**
```json
{"operation": "slice", "text": "Alignment with 4 rows and 3 columns\n-CG Beta\nAGG Gamma\nACG Delta\n-GG Epsilon"}
```

The doctest of MultipleSeqAlignment.__getitem__: align[1:5, 3:6].

## See Also

- [pairwise_aligner](pairwise_aligner.md)
- [alignment_operation](alignment_operation.md)
