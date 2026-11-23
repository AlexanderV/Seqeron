# alignment_summary_info

Bio.Align.AlignInfo.SummaryInfo of Biopython 1.79: dumb / gap consensus, replacement dictionary, position-specific score matrix, information content (ic_vector), get_column, print_info_content.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `alignment_summary_info` |
| **Method ID** | `SummaryInfo` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

One member of Bio.Align.AlignInfo.SummaryInfo of Biopython 1.79 (the last release that has the class; 1.88 keeps only get_column) on an alignment of records: dumb_consensus (per column the unique most common residue when its share of the residues, gaps not counted, reaches threshold, else the ambiguous letter; requireMultiple), gap_consensus (the same with the gap symbols counted), replacement_dictionary (letters: for every ordered pair of letters the summed product of the 'weight' annotations of the two rows that carry them in one column, the earlier row's letter first), pos_specific_score_matrix (the PSSM: per column the weighted count of every letter, with the consensus or axisSeq as the axis; charsToIgnore, the gap symbol is always ignored; the table text is str(pssm)), information_content (sum over the columns start..end-1 of the sum of f ln(f / e) / ln(logBase) with f the weighted frequency of the letter, gaps and charsToIgnore not counted, e from eFreqTableJson, pseudoCount k: f = (count + e k) / (total + k); the per-column values are ic_vector), ic_vector (the vector after information_content), print_info_content (position, letter of row repRecord, value with three decimals; needs eFreqTableJson), get_column. Errors and their messages are Biopython's (information_content without an expected-frequency table fails as it does there: TypeError; a negative pseudo count, a letter missing in the table, a range outside the first row, logBase 1 ...). VERIFICATION: the values of Biopython's own test suite (Tests/test_AlignInfo.py: ic 7.32029999423075, 133.061475107, 7.546, the printed PSSMs) and hand-derived values; the recorded run of 6 199 calls on the real class is the user-run package handoff/R14-13-aligninfo-biopython179 (its README says whether it has been integrated).

## Core Documentation Reference

- Source: [AlignInfo.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/AlignInfo.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md)
- Oracle: [alignment_summary_info_reference.py](../../../../scripts/alignment/alignment_summary_info_reference.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `operation` | string | Yes | The operation: dumb_consensus, gap_consensus, replacement_dictionary, pos_specific_score_matrix, information_content, ic_vector, print_info_content, get_column. |
| `recordsJson` | string | Yes | The rows as a JSON array of records {"id": ..., "seq": "ACGT-", "annotations": {"weight": 0.5}}; only seq is required (id defaults to <unknown id>); the 'weight' annotation of a row weights it (default 1.0). |
| `threshold` | number | No | consensus: share the most common residue needs (default 0.7). |
| `ambiguous` | string | No | consensus: the letter of a column without a clear consensus (default X). |
| `requireMultiple` | boolean | No | consensus: a column with one residue only is ambiguous. (default: false) |
| `letters` | string | No | replacement_dictionary: the letters to count (required). |
| `axisSeq` | string | No | pos_specific_score_matrix: the axis sequence (default the dumb consensus). |
| `charsToIgnoreJson` | string | No | pos_specific_score_matrix / information_content: a JSON array of letters that are not counted, e.g. ["-", "*"]. |
| `start` | integer | No | information_content: first column (0-based, default 0). (default: 0) |
| `end` | integer | No | information_content: end column, exclusive (default the length of the first row). |
| `eFreqTableJson` | string | No | information_content: expected frequency of each letter as a JSON object {"A": 0.25, "C": 0.25, ...} (gap symbols have none). |
| `logBase` | number | No | information_content: logarithm base (default 2: bits). |
| `pseudoCount` | number | No | information_content: pseudo count k (default 0). |
| `col` | integer | No | get_column: the column (negative from the end). (default: 0) |
| `repRecord` | integer | No | print_info_content: the row whose letters label the positions (default 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `operation` | string | The operation |
| `resultJson` | string | The value as JSON: a consensus string, [[a, b, weight], ...], {rows, text} of the PSSM, {value, ic_vector}, a column |
| `text` | string | The printed table of the PSSM or print_info_content |

## Errors

| Code | Message |
|------|---------|
| 1001 | Unknown operation, or an alignment that is not a JSON array of records |
| 1001 | information_content without an expected-frequency table (Biopython's TypeError), a range outside the first row, a negative pseudo count, a letter missing in the table |
| 1001 | A residue that is not among the counted letters (PSSM) |

## Examples

### Example 1: The consensus of three proteins

**User Prompt:**
> Give the dumb consensus of this protein alignment.

**Expected Tool Call:**
```json
{"tool": "alignment_summary_info", "arguments": {"operation": "dumb_consensus", "recordsJson": "[{\"id\": \"ID001\", \"seq\": \"MHQAIFIYQIGYP*LKSGYIQSIRSPEYDNW-\"}, {\"id\": \"ID002\", \"seq\": \"MH--IFIYQIGYAYLKSGYIQSIRSPEY-NW*\"}, {\"id\": \"ID003\", \"seq\": \"MHQAIFIYQIGYPYLKSGYIQSIRSPEYDNW*\"}]"}}
```

**Response:**
```json
{"operation": "dumb_consensus", "resultJson": "\"MHQAIFIYQIGYXXLKSGYIQSIRSPEYDNW*\""}
```

The value of Tests/test_AlignInfo.py of Biopython 1.79.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
