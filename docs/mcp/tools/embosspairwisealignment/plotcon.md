# plotcon

EMBOSS 6.6.0 plotcon: the similarity of a multiple alignment along its length (the data of the plot, -graph data).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `plotcon` |
| **Method ID** | `EmbossProgramRunner.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 plotcon: plots the similarity of a multiple alignment along its length (the data of the plot: x = the column, y = the similarity; -graph data). Every column scores the sum over all pairs of sequences of the integer matrix score of the two residues (scorefile: EBLOSUM62 for protein, EDNAFULL for nucleic sets, or a matrix file in files; a residue that the matrix does not contain, such as a gap, scores nothing) times the weight of the sequence; the score of the window of winsize columns that ends at column k is the sum of the column scores of every row over the window, put at the middle of the window, divided by numseq (numseq - 1) winsize^2 in single precision as the program computes it (a set of one sequence or winsize 0 gives NaN: the data file prints -000nan as the program does); the first winsize / 2 columns and the last ones have no window and stay 0. sequences is the text of the alignment (any input format of seqret) or the sequences; shorter sequences are padded with gaps. The result text is the data file plotcon1.dat (the other files and messages in files / messages); graph qualifiers (-gtitle, -gxtitle, -goutfile, -gdirectory ...) are accepted as typed on the command line, the title and the y label are set by the program. VERIFICATION: standard output, error stream, exit status and written files equal to the installed plotcon 6.6.0 on 400 recorded runs plus a fresh check of 1 500 runs (every input format, ranges, matrices, window sizes including 0 and negative, the graph devices that write text, the messages of the library).

## Core Documentation Reference

- Source: [EmbossSequencePrograms.Plotcon.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Emboss.PairwiseAlignment/EmbossPairwisePrograms.Plotcon.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md)
- Oracle: [generate_emboss_prog_wp63_oracle.py](../../../../scripts/alignment/generate_emboss_prog_wp63_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | string | Yes | The text of the alignment (any input format of seqret, e.g. aligned FASTA, Clustal, MSF), or sequences as text; at least one sequence. |
| `winsize` | integer | No | Number of columns to average the similarity over (-winsize, default 4). |
| `scorefile` | string | No | The scoring matrix (-scorefile): an EMBOSS matrix name (EBLOSUM62 default for protein, EDNAFULL default for nucleic acids, EBLOSUM30 ... EBLOSUM90, EPAM10 ... EPAM500, EDNAMAT, EDNASIMPLE, ENUC.4.4) or the name of a matrix file given in files. |
| `sequenceFormat` | string | No | The input format (-sformat); default: the program finds the format. |
| `qualifiers` | string | No | Sequence input qualifiers as typed after the sequence on the command line: -sbegin -send -sreverse, -slower -supper, -snucleotide -sprotein, -sid name, -sdbname db, -squick, -iquery query, -ioffset bytes. |
| `graphQualifiers` | string | No | Graph qualifiers as typed on the command line: -gtitle T -gsubtitle T -gxtitle T -gytitle T -gdesc T -goutfile NAME -gdirectory DIR. |
| `files` | object | No | More files by name (file name -> text) that the program may read (a matrix file, a file named by a USA, @file). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The data file plotcon1.dat of -graph data (the title, the extent and x / y of every column) |
| `files` | object | The other files the program wrote, by name |
| `messages` | string[] | What the program wrote to its error stream |

## Errors

| Code | Message |
|------|---------|
| 1001 | The sequences cannot be empty |
| 1001 | The program died: the messages of the library (a matrix it cannot read, a value it refuses, an input no reader reads) |
| 1001 | An input on which the reference program has no defined result is refused with the reason |

## Examples

### Example 1: Three proteins, window of 4 columns

**User Prompt:**
> Plot the similarity of this alignment.

**Expected Tool Call:**
```json
{"tool": "plotcon", "arguments": {"sequences": ">a\nMKTAYIAKQRQISFVKSH\n>b\nMKTAYLAKQRQISFVRSH\n>c\nMKSAYIAK-RQISFVKSH\n"}}
```

**Response:**
```json
{"text": "##2D Plot ... 3.000000\t1.020833 ... (18 points)"}
```

The first window (columns 1-4) has the pair scores 2 x 15 + 2 x 15 + 2 x 7 + 2 x 12 = 98 (BLOSUM62: M/M 5, K/K 5, T/T 5, T/S 1, A/A 4) divided by 3 x 2 x 4 x 4 = 96: 1.020833 at point 3.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
