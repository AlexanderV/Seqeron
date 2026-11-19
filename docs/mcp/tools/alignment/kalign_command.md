# kalign_command

Kalign 3.3.5 as a function of its command line and the text of the files it names.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `kalign_command` |
| **Method ID** | `KalignProgram.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Kalign 3.3.5 as a function of its command line and the text of the files it names (KalignProgram.Run): the arguments are those of the program (-i / --input / --infile / --in, -o / --output / --outfile / --out, -f / --format fasta, --type rna|dna|internal|protein|divergent, --gpo, --gpe, --tgpe, -n / --nthreads, -q, --set, -h, -v / -V, --showw; getopt_long_only forms: one or two dashes, --name=value, unambiguous abbreviations; non-option arguments are input files, several are merged). The type word is matched like the program's strstr chain (rna, dna, internal, protein, divergent in this order: 'internal' contains 'rna' and therefore means RNA). The alignment is written as FASTA with 60 residues per line to the -o file or to standard output. Exit code 1 = the program stops with an error (unknown / ambiguous option, bad type word or format, fewer than one thread, a missing file, a type that does not fit the sequences); the program segfaults on a single sequence, the port reports exit code 1. --format msf / clu and MSF / Clustal inputs are accepted by the program but not yet written / read here (a NotSupportedException error). VERIFICATION: stdout / output file and exit status equal the installed Debian kalign 1:3.3.5-1 on 83 / 83 committed command lines (scripts/alignment/generate_kalign_command_oracle.py).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/KalignProgram.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/KalignProgram.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md) (section Kalign 3.3.5)
- Oracle: `scripts/alignment/generate_kalign_command_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/kalign_command_oracle.tsv`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `arguments` | array<string> | Yes | The command-line arguments of kalign, one per element (for example ["-i", "in.fa", "-o", "out.fa", "--type", "dna", "-q"]). |
| `files` | object<string,string> | No | The text of the files the arguments name: file name -> FASTA text. |
| `stdin` | string | No | Text of standard input, read when no input file is named. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `exitCode` | integer | 0 = success, 1 = the program stops with an error. |
| `stdout` | string | Standard output (the alignment when no -o file is named, the version / help text). |
| `files` | object<string,string> | The files the program wrote (-o): name -> text. |
| `messages` | array<string> | What the program prints to the error stream (errors). |

## Errors

| Code | Message |
|------|---------|
| 1001 | arguments is null |
| 1002 | --format msf / clu output (accepted by the program, not written here): NotSupportedException |

## Examples

### Example 1: FASTA in, FASTA out, DNA type

**Tool Call:**
```json
{
  "tool": "kalign_command",
  "arguments": {
    "arguments": [
      "-q",
      "-i",
      "in.fa",
      "-o",
      "out.fa",
      "--type",
      "dna"
    ],
    "files": {
      "in.fa": ">s1 first\nACGTACGTTAGCTAGCTAGCATCGATCGATGCATGCATCG\n>s2\nACGTACGTAGCTAGCTAGCATCGATCGATGCATGCCATCG\n>s3\nACGTTCGTTAGCTAGCTTGCATCGATCGATGCATGCATCG\n"
    }
  }
}
```

**Response:**
```json
{
  "exitCode": 0,
  "stdout": "",
  "files": {
    "out.fa": ">s1 first\nACGTACGTTAGCTAGCTAGCATCGATCGATGCATGC-ATCG\n>s2\nACGTACG-TAGCTAGCTAGCATCGATCGATGCATGCCATCG\n>s3\nACGTTCGTTAGCTAGCTTGCATCGATCGATGCATGC-ATCG\n"
  },
  "messages": []
}
```

`kalign -q -i in.fa -o out.fa --type dna` writes exactly this out.fa.

### Example 2: Ambiguous abbreviation: the program stops

**Tool Call:**
```json
{
  "tool": "kalign_command",
  "arguments": {
    "arguments": [
      "-q",
      "-i",
      "in.fa",
      "--gp",
      "3"
    ],
    "files": {
      "in.fa": ">s1 first\nACGTACGTTAGCTAGCTAGCATCGATCGATGCATGCATCG\n>s2\nACGTACGTAGCTAGCTAGCATCGATCGATGCATGCCATCG\n>s3\nACGTTCGTTAGCTAGCTTGCATCGATCGATGCATGCATCG\n"
    }
  }
}
```

**Response:**
```json
{
  "exitCode": 1,
  "stdout": "",
  "files": {},
  "messages": [
    "option '--gp' is ambiguous"
  ]
}
```

The installed program prints `option '--gp' is ambiguous; possibilities: '--gpo' '--gpe'` and exits with status 1.

## See Also

- [multiple_align_kalign](multiple_align_kalign.md)
- [clustalw_command](clustalw_command.md)
- [muscle_msa_command](muscle_msa_command.md)
