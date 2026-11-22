# wordfinder

EMBOSS 6.6.0 wordfinder: Matches large sequences against one or more other sequences by their words, as wordfinder of EMBOSS 6.6.0.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `wordfinder` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 wordfinder: the target sequences (bsequence) are indexed by their words of wordlen residues, every query sequence (asequence) is looked up in the index, the diagonal with the most word residues fixes the start points and a banded Smith-Waterman alignment of width cells (gap penalties gapopen / gapextend, matrix datafile: an EMBOSS matrix name such as EDNAFULL or EBLOSUM62, or the text of a matrix file) is written in the alignment format alignFormat (every -aformat of EMBOSS, default simple, with the sub-header lines Wordscore and Alignlength) as out.txt; the matches and alignments excluded by limitmatch / lowmatch / limitalign / lowalign and the matches per target are written to err.txt. Equals the installed program on 300 runs (every input format of seqret, ranges, the report and graph options, the messages of the library; standard output, error stream, exit status and written files byte for byte). B09 WP63, F282.

## Core Documentation Reference

- Source: [EmbossSequencePrograms.Wordfinder.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/EmbossSequencePrograms.Wordfinder.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_wp63_oracle.py](../../../../scripts/alignment/generate_emboss_prog_wp63_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The text of the main output of the program (the report, or its standard output). |
| `files` | object | The other files the program wrote (sequences, features, graph data), by name. |
| `messages` | string[] | What the library printed on the error stream (warnings and errors). |

## Errors

| Code | Message |
|------|---------|
| 1001 | The sequence cannot be empty |
| 1001 | The program died: the messages of the library (a value it refuses, an input no reader reads, a sequence of the wrong kind) |
| 1001 | An input on which the reference program never stops or dies by a fault is refused |

## Examples

### Example 1: A run on a small input

**User Prompt:**
> Run wordfinder on this input.

**Expected Tool Call:**
```json
{"tool": "wordfinder", "arguments": {"querySequence": ">q1\nACGTACGTACGGGTTTTAAACCCGGGATATATCGCGCGTAGCTAGCTAGCTAGGGGATCCC\n", "targetSequence": ">t1\nGGGACGTACGTACGGGTTTTAAACCCGGGATATATCGCGCGTAGCTAGCTAGCTAGGGGATCCCAAAAAA\n"}}
```

**Response:**
```json
{"text": "########################################\n# Program: wordfinder\n# Rundate: Tue 14 Nov 2023 22:13:20\n# Commandline: wordfinder\n#    -auto\n#    -asequence asequence.dat\n#    -bsequence bsequence.dat\n#    -outfile out.txt\n#    -errfile err.txt\n# Align_format: simple\n# Report_file: out.txt\n########################################\n\n#=======================================\n#\n# Aligned_sequences: 2\n# 1: q1\n# 2: t1\n# Matrix: EDNAFULL\n# Gap_penalty: 30.0\n# Extend_penalty: 1.5\n#\n# Length: 61\n# Identity:      61/61 (100.0%)\n# Similarity:    61/61 (100.0%)\n# Gaps:           0/61 ( 0.0%)\n# Score: 305.0\n# \n# Wordscore:61\n# Alignlength:61\n#\n#=======================================\n\nq1                 1 ACGTACGTACGGGTTTTAAACCCGGGATATATCGCGCGTAGCTAGCTAGC     50\n                     ||||||||||||||||||||||||||||||||||||||||||||||||||\nt1                 4 ACGTACGTACGGGTTTTAAACCCGGGATATATCGCGCGTAGCTAGCTAGC     53\n\nq1                51 TAGGGGATCCC     61\n                     |||||||||||\nt1                54 TAGGGGATCCC     64\n\n\n#---------------------------------------\n#---------------------------------------\n", "files": {"err.txt": "Target 1 t1 matches 1\n"}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
