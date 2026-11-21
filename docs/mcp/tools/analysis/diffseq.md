# diffseq

EMBOSS 6.6.0 diffseq: the regions of difference between two nearly identical sequences (336 nucleotide / protein pairs, every -rformat: report and both GFF3 feature files byte for byte). The exact word matches of wordSize symbols (the engine and the minimal non-overlapping set of megamerger) are the regions of identity; the stretches between them are the differences (a difference with length 0 in one sequence is an insertion in the other); with globalDifferences the differences before the first and after the last match are reported too (a start or end that matches by itself is not counted). A nucleotide pair also gets the SNP count and its transitions / transversions. Returns the differences and the program's files: the report in the default format 'diffseq' (7 report tags; reportFormat takes any -rformat or the report qualifiers as typed, e.g. '-rformat table -raccshow'), and the GFF3 files -aoutfeat / -boutfeat. The sequence qualifiers (-sbegin1 ...) and feature-table inputs (GenBank / EMBL with features) are not part of this port.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `diffseq` |
| **Method ID** | `GenomicAnalyzer.Diffseq` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 diffseq: the regions of difference between two nearly identical sequences (336 nucleotide / protein pairs, every -rformat: report and both GFF3 feature files byte for byte). The exact word matches of wordSize symbols (the engine and the minimal non-overlapping set of megamerger) are the regions of identity; the stretches between them are the differences (a difference with length 0 in one sequence is an insertion in the other); with globalDifferences the differences before the first and after the last match are reported too (a start or end that matches by itself is not counted). A nucleotide pair also gets the SNP count and its transitions / transversions. Returns the differences and the program's files: the report in the default format 'diffseq' (7 report tags; reportFormat takes any -rformat or the report qualifiers as typed, e.g. '-rformat table -raccshow'), and the GFF3 files -aoutfeat / -boutfeat. The sequence qualifiers (-sbegin1 ...) and feature-table inputs (GenBank / EMBL with features) are not part of this port.

## Core Documentation Reference

- Algorithm: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence (nucleotide or protein letters). |
| `sequence2` | string | Yes | Second sequence (the same kind). |
| `wordSize` | integer | No | Word size (diffseq -wordsize, default 10, at least 2). (default: 10) |
| `globalDifferences` | boolean | No | Report the differences at the start and the end (diffseq -globaldifferences). (default: false) |
| `name1` | string | No | Name of the first sequence. (default: seqa) |
| `name2` | string | No | Name of the second sequence. (default: seqb) |
| `reportFormat` | string | No | The report format (-rformat name) or the report qualifiers as typed after the files (e.g. '-rformat table -raccshow'); default format diffseq. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `differences` | array | The regions of difference (1-based starts / ends, lengths; a length 0 is a gap). |
| `overlapStart1` | integer | Start of the overlap in the first sequence (-1: no match). |
| `overlapEnd1` | integer | End of the overlap in the first sequence. |
| `overlapStart2` | integer | Start of the overlap in the second sequence. |
| `overlapEnd2` | integer | End of the overlap in the second sequence. |
| `snps` | integer | SNP count (nucleotide pairs). |
| `transitions` | integer | Transitions among the SNPs. |
| `transversions` | integer | Transversions among the SNPs. |
| `report` | string | The report file text (-rformat, default diffseq). |
| `featureFile1` | string | The GFF3 feature file of the first sequence (-aoutfeat). |
| `featureFile2` | string | The GFF3 feature file of the second sequence (-boutfeat). |

## Errors

| Code | Message |
|------|---------|
| InvalidParams | Invalid arguments (an empty or unreadable sequence, a word size below 2, a report format the program refuses for the sequence type). |

## Examples

### Example 1: Compare two nearly identical sequences

**User Prompt:**
> Compare two nearly identical sequences

**Expected Tool Call:**
```json
{"tool": "diffseq", "arguments": {"sequence1": "ACGTACGTACGTAAGGCCTTAACCGGTT", "sequence2": "ACGTACGTACGTAAGGCATTAACCGGTT", "wordSize": 8}}
```

**Response:**
```json
{"differences": [{"start1": 18, "end1": 18, "length1": 1, "start2": 18, "end2": 18, "length2": 1}], "snps": 1, "transitions": 1, "transversions": 0}
```

## See Also

- [megamerger](megamerger.md) - merge by word matches
- [seqmatchall](seqmatchall.md) - all-against-all word matches
