# spoa_command

command line (FASTA / FASTQ files)

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `spoa_command` |
| **Method ID** | `SpoaProgram.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

POA 2.0 as a function of its command line and the text of the files it names (SpoaProgram.Run): poa [-read_fasta FILE] [-read_msa FILE] [-read_msa2 FILE] [-read_msa_list FILE] [-subset FILE] [-subset2 FILE] [-remove FILE] [-remove2 FILE] [-do_global] [-do_progressive] [-read_pairscores FILE] [-fuse_all] [-hb] [-hbmin X] [-best] [-preserve_seqorder] [-tolower] [-toupper] [-pir FILE] [-clustal FILE] [-po FILE] [-printmatrix LETTERS] [-v] MATRIXFILE, with the program's argument loop (every argument is tested against every flag; unknown flags are ignored; the first argument that does not start with a dash is the matrix file). Output: the PIR / CLUSTAL / PO files named by -pir / -clustal / -po, standard output (-v gap penalty tables, -printmatrix, the permutation line of -preserve_seqorder). Exit code 1 = the program stops with an error (no matrix or unreadable matrix, no sequence file, no sequences, conflicting flags), 255 = called without arguments (usage). The alignment progress lines of standard error are not reproduced (the program prints uninitialised %ld arguments there); the '...Read N sequences...' / '...Wrote N sequences...' lines are. The MSA / PO inputs (-read_msa, -read_msa2, -read_msa_list; the format is detected from the first line: a PO file, a FASTA-PIR alignment or a CLUSTAL alignment; rows are read in 510-character pieces like the program's fgets) and the filters -subset / -subset2 / -remove / -remove2 (SOURCENAME= lines) are ported; an alignment row name beyond 31 characters, a PO file whose links / sources / rings are inconsistent and a row without residues stop the run (exit 1 / 255) where the program aborts, crashes or loops. VERIFICATION: exit code, standard output, the Read / Wrote lines and the three output files equal the installed Debian poa 2.0+20060928 on 359 / 359 FASTA command lines and 687 / 687 MSA / PO input command lines (scripts/alignment/generate_poa_msa_oracle.py).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SpoaProgram.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SpoaProgram.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md) (section 4.15 SPOA)
- Oracle: `scripts/alignment/generate_spoa_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/spoa_oracle.tsv`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `arguments` | array<string> | Yes | The command-line arguments of spoa, one per element (for example ["-l", "1", "-r", "2", "in.fa"]). |
| `files` | object<string,string> | No | The text of the input files the arguments name: file name -> text (FASTA / FASTQ; bytes are read as Latin-1). |
| `binaryFiles` | object<string,string> | No | Binary input files (for example gzip-compressed FASTA): file name -> base64. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `exitCode` | integer | 0 = success, 1 = the program stops with an error (options, files, penalties), 134 = abort (invalid FASTA / FASTQ, empty input). |
| `stdout` | string | Standard output: consensus / MSA / GFA text, --help, --version. |
| `stderr` | string | Error messages and the `terminate called` text of an abort. |
| `files` | object<string,string> | The files the program wrote: the DOT file of -d. |

## Errors

| Code | Message |
|------|---------|
| 1001 | arguments is null |
| 1002 | (none: a refused command line is an exit status 1 or 134 with the reason in `stderr`) |

## Examples

### Example 1: Global multiple sequence alignment (spoa -l 1 -r 1 u.fa)

**Tool Call:**
```json
{
  "tool": "spoa_command",
  "arguments": {
    "arguments": [
      "-l",
      "1",
      "-r",
      "1",
      "u.fa"
    ],
    "files": {
      "u.fa": ">r1\nACGTACGTAC\n>r2\nACGTTCGTAC\n>r3\nACGACGTAC\n"
    }
  }
}
```

**Response:**
```json
{
  "exitCode": 0,
  "stdout": ">r1\nACGTACGTAC\n>r2\nACGTTCGTAC\n>r3\nACG-ACGTAC\n",
  "stderr": "",
  "files": {}
}
```

### Example 2: Consensus and the DOT graph (spoa -d g.dot u.fa)

**Tool Call:**
```json
{
  "tool": "spoa_command",
  "arguments": {
    "arguments": [
      "-d",
      "g.dot",
      "u.fa"
    ],
    "files": {
      "u.fa": ">r1\nACGTACGTAC\n>r2\nACGTTCGTAC\n>r3\nACGACGTAC\n"
    }
  }
}
```

**Response:**
```json
{
  "exitCode": 0,
  "stdout": ">Consensus LN:i:10\nACGTACGTAC\n",
  "stderr": "",
  "files": {
    "g.dot": "digraph 3 {\n  graph [rankdir = LR]\n  0[label = \"0 - A\", style = filled, fillcolor = goldenrod1]\n  ..."
  }
}
```

## See Also

- [poa_command](poa_command.md)
- [multiple_align_poa](multiple_align_poa.md)
- [kalign_command](kalign_command.md)
