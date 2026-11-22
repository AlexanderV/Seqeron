# tcode

EMBOSS 6.6.0 tcode: Identifies protein-coding regions of a nucleotide sequence with the Fickett TESTCODE statistic, as tcode of EMBOSS 6.6.0.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `tcode` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 tcode: the Fickett TESTCODE statistic of the windows of window bases (at least 200) every step bases, from the base compositions and the position asymmetries of the three codon positions with the probability tables of the data file Etcode.dat (or the datafile given); every window is reported with its score and the estimation Non-coding (below 0.74), No opinion or Coding (0.95 and above) as the report out.txt in the report format reportFormat (default table: every report format of EMBOSS, the report qualifiers of reportQualifiers); with plot the values are plotted against the window mid positions (graph: data writes the data file tcode1.dat, bedgraph the track file tcode.bedgraph, none draws nothing; the devices that draw with plplot (ps, png, svg ...) are refused). Equals the installed program on 300 runs (every input format of seqret, ranges, the report and graph options, the messages of the library; standard output, error stream, exit status and written files byte for byte). B09 WP63, F282.

## Core Documentation Reference

- Source: [EmbossSequencePrograms.Tcode.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/EmbossSequencePrograms.Tcode.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_wp63_oracle.py](../../../../scripts/alignment/generate_emboss_prog_wp63_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file (any input format of seqret), or a sequence as letters. |
| `window` | integer | No | The length of the sliding window (-window, default 200: at least 200 and at most the length of the first sequence). |
| `step` | integer | No | The step of the window (-step, default 3). |
| `dataFile` | string | No | The text of a testcode data file in the format of Etcode.dat (-datafile; default the file of EMBOSS). |
| `plot` | boolean | No | Plot the statistic (-plot). (default: false) |
| `graph` | string | No | The graph device of -plot (-graph): data, bedgraph or none. |
| `reportFormat` | string | No | The report format (-rformat, default table): simple, embl, genbank, gff, gff3, pir, swiss, dasgff, debug, listfile, dbmotif, diffseq, draw, restrict, excel, feattable, motif, nametable, regions, seqtable, srs, summary, table, tagseq. |
| `reportQualifiers` | string | No | Report qualifiers as typed on the command line: -raccshow -rdesshow -rusashow -norscoreshow -norstrandshow -rmaxall N -rmaxseq N (-rformat is the parameter reportFormat). |
| `graphQualifiers` | string | No | Graph qualifiers as typed on the command line: -gtitle T -gsubtitle T -gxtitle T -gytitle T -gdesc T -goutfile NAME -gdirectory DIR. |
| `sequenceFormat` | string | No | The input format (-sformat1); default: the program finds the format. |
| `qualifiers` | string | No | Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes. |
| `files` | object | No | More files by name (file name -> text) that the program may read (a file named by a USA, @file). |

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
> Run tcode on this input.

**Expected Tool Call:**
```json
{"tool": "tcode", "arguments": {"sequence": ">s1 test\nATGAAATTTGGGCCCGATGAACTGACCGGCTTAGGTACCATGAAATTTGGGCCCGATGAACTGACCGGCTTAGGTACCATGAAATTTGGGCCCGATGAACTGACCGGCTTAGGTACCATGAAATTTGGGCCCGATGAACTGACCGGCTTAGGTACCATGAAATTTGGGCCCGATGAACTGACCGGCTTAGGTACCATGAAATTTGGGCCCGATGAACTGACCGGCTTAGGTACCATGAAATTTGGGCCCGATGAACTGACCGGCTTAGGTACCATGAAATTTGGGCCCGATGAACTGACCGGCTTAGGTACC\n", "window": 200, "step": 30}}
```

**Response:**
```json
{"text": "########################################\n# Program: tcode\n# Rundate: Tue 14 Nov 2023 22:13:20\n# Commandline: tcode\n#    -auto\n#    -sequence sequence.dat\n#    -window 200\n#    -step 30\n#    -outfile out.txt\n# Report_format: table\n# Report_file: out.txt\n########################################\n\n#=======================================\n#\n# Sequence: s1     from: 1   to: 312\n# HitCount: 4\n#\n# Fickett TESTCODE statistic\n#\n#=======================================\n\n  Start     End  Strand   Score Estimation\n      1     200       +   1.048 Coding    \n     31     230       +   1.034 Coding    \n     61     260       +   1.087 Coding    \n     91     290       +   1.046 Coding    \n\n#---------------------------------------\n#---------------------------------------\n\n#---------------------------------------\n# Total_sequences: 1\n# Total_length: 312\n# Reported_sequences: 1\n# Reported_hitcount: 4\n#---------------------------------------\n", "files": {}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
