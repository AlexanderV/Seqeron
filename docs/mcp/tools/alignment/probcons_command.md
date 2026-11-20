# probcons_command

command line (options, MFA / MSF / Clustal input files)

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `probcons_command` |
| **Method ID** | `ProbConsProgram.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

ProbCons 1.12 as a function of its command line and the text of the files it names (ProbConsProgram.Run): -clustalw, -c / --consistency, -ir / --iterative-refinement, -pre / --pre-training, -pairs, -viterbi, -v / --verbose, -annot FILE, -t / --train FILE, -e / --emissions, -p / --paramfile FILE, -a / --alignment-order, -go / --gap-open, -ge / --gap-extension, -co / --cutoff, with the program's number reading (strtol base 0, strtod prefix), ranges and error texts; input files are FASTA (any of '.', '-' and white space inside the sequence are removed; an empty sequence ends the file), MSF or Clustal alignments (ParseMSF with its quirks: letters upper-cased, gaps kept - the program then overflows a buffer, so gapped alignment input is not verified). Output: the alignment as FASTA (60 columns) or Clustal-like text with the conservation line (stdout), the pairwise alignment files of -pairs / -viterbi (header1-header2.fasta / .aln), the -annot column scores (setw(4) per line), the trained parameter file of -t (%.10f) - and the program's standard-error text byte for byte (heading, parameter set at %.10g, 'Loading sequence file', 'Alignment tree: ...', -v progress). Exit code 1 = the program exits with an error (unknown option, bad value, missing file, unknown character, -e without -t); with no input file the Debian build prints 'Zero sequences delivered to function DoAlign()' and exits with 0. VERIFICATION: exit code, stdout, stderr and files equal the installed Debian probcons 1.12-14 on 407 / 407 committed command lines.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/ProbConsProgram.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ProbConsProgram.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md) (section ProbCons 1.12)
- Oracle: `scripts/alignment/generate_probcons_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/probcons_oracle.tsv`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `arguments` | array<string> | Yes | The command-line arguments of probcons, one per element (for example ["-clustalw", "-c", "1", "in.fa"]). |
| `files` | object<string,string> | No | The text of the files the arguments name: file name -> text (sequence files and the -p parameter file). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `exitCode` | integer | 0 = success, 1 = the program exits with an error. |
| `stdout` | string | The alignment (FASTA with 60 columns, or the Clustal-like text with -clustalw). |
| `stderr` | string | Everything the program writes to the error stream: heading, parameter set, "Loading sequence file", pre-training sets, -v progress, the alignment tree, errors. |
| `files` | object<string,string> | The files the program wrote: pairwise alignments (-pairs / -viterbi), the -annot file, the trained parameter file (-t). |

## Errors

| Code | Message |
|------|---------|
| 1001 | arguments is null |

## Examples

### Example 1: Clustal output with ten refinement passes

**Tool Call:**
```json
{
  "tool": "probcons_command",
  "arguments": {
    "arguments": [
      "-clustalw",
      "-ir",
      "10",
      "in.fa"
    ],
    "files": {
      "in.fa": ">a one\nMKVLAAGIVGLLLAQPSTAEEKV\n>b\nMKVLAGIVGLLAQPSTAEEKVLL\n>c\nMKILAAGIAGLLLQPSAEEKV\n>d\nMRVLAAGLVGLLLAQSTAEEK\n"
    }
  }
}
```

**Response:**
```json
{
  "exitCode": 0,
  "stdout": "PROBCONS version 1.12 multiple sequence alignment\n\na    MKVLAAGIVGLLLAQPSTAEEK--V\nb    MKVL-AGIVGLLA-QPSTAEEKVLL\nc    MKILAAGIAGLLL-QPS-AEEK--V\nd    MRVLAAGLVGLLL-AQSTAEEK---\n     *::****:.***    ******  :\n",
  "stderr": "\nPROBCONS version 1.12 - align multiple protein sequences and print to standard output\nWritten by Chuong Do\n\nUsing parameter set:\n    initDistrib[] = { 0.6814756989 8.615339902e-05 8.615339902e-05 0.1591759622 0.1591759622 }\n        gapOpen[] = { 0.0119511066 0.0119511066 0.008008334786 0.008008334786 }\n      gapExtend[] = { 0.3965826333 0.3965826333 0.8988758326 0.8988758326 }\n\nLoading sequence file: in.fa\nAlignment tree: ((a one (c d)) b)\n\n\n",
  "files": {}
}
```

The stderr text is byte for byte the program's (parameter set printed with 10 significant digits).

### Example 2: All-pairs Viterbi alignments (probcons -viterbi in.fa)

**Tool Call:**
```json
{
  "tool": "probcons_command",
  "arguments": {
    "arguments": [
      "-viterbi",
      "in.fa"
    ],
    "files": {
      "in.fa": ">a one\nMKVLAAGIVGLLLAQPSTAEEKV\n>b\nMKVLAGIVGLLAQPSTAEEKVLL\n>c\nMKILAAGIAGLLLQPSAEEKV\n>d\nMRVLAAGLVGLLLAQSTAEEK\n"
    }
  }
}
```

**Response:**
```json
{
  "exitCode": 0,
  "stdout": "",
  "stderr": "\nPROBCONS version 1.12 - align multiple protein sequences and print to standard output\nWritten by Chuong Do\n\nUsing parameter set:\n    initDistrib[] = { 0.6814756989 8.615339902e-05 8.615339902e-05 0.1591759622 0.1591759622 }\n        gapOpen[] = { 0.0119511066 0.0119511066 0.008008334786 0.008008334786 }\n      gapExtend[] = { 0.3965826333 0.3965826333 0.8988758326 0.8988758326 }\n\nLoading sequence file: in.fa\n",
  "files": {
    "a one-b.fasta": ">a one\nMKVLAAGIVGLLLAQPSTAEEKV--\n>b\nMKVL-AGIVGL-LAQPSTAEEKVLL\n",
    "a one-c.fasta": ">a one\nMKVLAAGIVGLLLAQPSTAEEKV\n>c\nMKILAAGIAGLLL-QPS-AEEKV\n",
    "a one-d.fasta": ">a one\nMKVLAAGIVGLLLAQPSTAEEKV\n>d\nMRVLAAGLVGLLLAQ-STAEEK-\n",
    "b-c.fasta": ">b\nMKVL-AGIVGLLAQPSTAEEKVLL\n>c\nMKILAAGIAGLLLQPS-AEEKV--\n",
    "b-d.fasta": ">b\nMKVL-AGIVGLLAQPSTAEEKVLL\n>d\nMRVLAAGLVGLLLAQSTAEEK---\n",
    "c-d.fasta": ">c\nMKILAAGIAGLLLQPSAEEKV\n>d\nMRVLAAGLVGLLLAQSTAEEK\n"
  }
}
```

One file per pair, named header1-header2.fasta.

## See Also

- [multiple_align_probcons](multiple_align_probcons.md)
- [kalign_command](kalign_command.md)
- [clustalw_command](clustalw_command.md)
