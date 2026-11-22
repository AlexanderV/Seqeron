# plotorf

EMBOSS 6.6.0 plotorf: Plots the potential open reading frames of a nucleotide sequence, as plotorf of EMBOSS 6.6.0.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `plotorf` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 plotorf: for each of the six frames the stretches from a start codon (start, default ATG) to the next stop codon of the frame (stop, default TAA,TAG,TGA; a stretch that reaches the end of the sequence ends there) as the filled rectangles of six plots; graph: data writes the six data files plotorf1.dat ... plotorf6.dat (returned in files), bedgraph one track file, none draws nothing (the devices that draw with plplot are refused). Equals the installed program on 60 runs (every input format of seqret, ranges, the report and graph options, the messages of the library; standard output, error stream, exit status and written files byte for byte). B09 WP63, F282.

## Core Documentation Reference

- Source: [EmbossSequencePrograms.Plotorf.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/EmbossSequencePrograms.Plotorf.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_emboss_prog_wp63_oracle.py](../../../../scripts/alignment/generate_emboss_prog_wp63_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file (any input format of seqret), or a sequence as letters. |
| `start` | string | No | The start codons, a comma list (-start, default ATG). |
| `stop` | string | No | The stop codons, a comma list (-stop, default TAA,TAG,TGA). |
| `graph` | string | No | The graph device (-graph): data (default), bedgraph or none. |
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
> Run plotorf on this input.

**Expected Tool Call:**
```json
{"tool": "plotorf", "arguments": {"sequence": ">o1\nATGAAATTTGGGCCCTAGATGCCCAAATTTTAAATGGGGGGGGTGATAGCATGAAACCCTGAGGATAA\n"}}
```

**Response:**
```json
{"text": "Created plotorf1.dat\nCreated plotorf2.dat\nCreated plotorf3.dat\nCreated plotorf4.dat\nCreated plotorf5.dat\nCreated plotorf6.dat\n", "files": {"plotorf1.dat": "##Multi 2D Plot Small\n##Title Potential codons (rectangles)\n##Graphs 6\n##Number 1\n##Points 0\n##XminA 1.000000 XmaxA 68.000000 YminA 0.000000 YmaxA 1.000000\n##Xmin 1.000000 Xmax 68.000000 Ymin 0.000000 Ymax 1.000000\n##ScaleXmin 1.000000 ScaleXmax 68.000000 ScaleYmin 0.000000 ScaleYmax 1.000000\n##Maintitle F1\n##Subtitle \n##Xtitle Sequence\n##Ytitle Orf\n##DataObjects\n##Number 3\nFilled Rectangle x1 16.000000 y1 0.000000 x2 1.000000 y2 1.000000 colour 4\nFilled Rectangle x1 31.000000 y1 0.000000 x2 19.000000 y2 1.000000 colour 4\nFilled Rectangle x1 52.000000 y1 0.000000 x2 34.000000 y2 1.000000 colour 4\n##GraphObjects\n##Number 0\n", "plotorf2.dat": "##Multi 2D Plot Small\n##Title Potential codons (rectangles)\n##Graphs 6\n##Number 2\n##Points 0\n##XminA 1.000000 XmaxA 68.000000 YminA 0.000000 YmaxA 1.000000\n##Xmin 1.000000 Xmax 68.000000 Ymin 0.000000 Ymax 1.000000\n##ScaleXmin 1.000000 ScaleXmax 68.000000 ScaleYmin 0.000000 ScaleYmax 1.000000\n##Maintitle F2\n##Subtitle \n##Xtitle Sequence\n##Ytitle Orf\n##DataObjects\n##Number 0\n##GraphObjects\n##Number 0\n", "plotorf3.dat": "##Multi 2D Plot Small\n##Title Potential codons (rectangles)\n##Graphs 6\n##Number 3\n##Points 0\n##XminA 1.000000 XmaxA 68.000000 YminA 0.000000 YmaxA 1.000000\n##Xmin 1.000000 Xmax 68.000000 Ymin 0.000000 Ymax 1.000000\n##ScaleXmin 1.000000 ScaleXmax 68.000000 ScaleYmin 0.000000 ScaleYmax 1.000000\n##Maintitle F3\n##Subtitle \n##Xtitle Sequence\n##Ytitle Orf\n##DataObjects\n##Number 1\nFilled Rectangle x1 60.000000 y1 0.000000 x2 51.000000 y2 1.000000 colour 4\n##GraphObjects\n##Number 0\n", "plotorf4.dat": "##Multi 2D Plot Small\n##Title Potential codons (rectangles)\n##Graphs 6\n##Number 4\n##Points 0\n##XminA 1.000000 XmaxA 68.000000 YminA 0.000000 YmaxA 1.000000\n##Xmin 1.000000 Xmax 68.000000 Ymin 0.000000 Ymax 1.000000\n##ScaleXmin 1.000000 ScaleXmax 68.000000 ScaleYmin 0.000000 ScaleYmax 1.000000\n##Maintitle R1\n##Subtitle \n##Xtitle Sequence\n##Ytitle Orf\n##DataObjects\n##Number 0\n##GraphObjects\n##Number 0\n", "plotorf5.dat": "##Multi 2D Plot Small\n##Title Potential codons (rectangles)\n##Graphs 6\n##Number 5\n##Points 0\n##XminA 1.000000 XmaxA 68.000000 YminA 0.000000 YmaxA 1.000000\n##Xmin 1.000000 Xmax 68.000000 Ymin 0.000000 Ymax 1.000000\n##ScaleXmin 1.000000 ScaleXmax 68.000000 ScaleYmin 0.000000 ScaleYmax 1.000000\n##Maintitle R2\n##Subtitle \n##Xtitle Sequence\n##Ytitle Orf\n##DataObjects\n##Number 1\nFilled Rectangle x1 1.000000 y1 0.000000 x2 52.000000 y2 1.000000 colour 4\n##GraphObjects\n##Number 0\n", "plotorf6.dat": "##Multi 2D Plot Small\n##Title Potential codons (rectangles)\n##Graphs 6\n##Number 6\n##Points 0\n##XminA 1.000000 XmaxA 68.000000 YminA 0.000000 YmaxA 1.000000\n##Xmin 1.000000 Xmax 68.000000 Ymin 0.000000 Ymax 1.000000\n##ScaleXmin 1.000000 ScaleXmax 68.000000 ScaleYmin 0.000000 ScaleYmax 1.000000\n##Maintitle R3\n##Subtitle \n##Xtitle Sequence\n##Ytitle Orf\n##DataObjects\n##Number 0\n##GraphObjects\n##Number 0\n"}, "messages": []}
```

The output is the text the installed program writes.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
