# poa_command

command line (FASTA / matrix / pair-score files)

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `poa_command` |
| **Method ID** | `PoaProgram.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

POA 2.0 as a function of its command line and the text of the files it names (PoaProgram.Run): poa -read_fasta FILE [-do_global] [-do_progressive] [-read_pairscores FILE] [-fuse_all] [-hb] [-hbmin X] [-best] [-preserve_seqorder] [-tolower] [-toupper] [-pir FILE] [-clustal FILE] [-po FILE] [-printmatrix LETTERS] [-v] MATRIXFILE, with the program's argument loop (every argument is tested against every flag; unknown flags are ignored; the first argument that does not start with a dash is the matrix file). Output: the PIR / CLUSTAL / PO files named by -pir / -clustal / -po, standard output (-v gap penalty tables, -printmatrix, the permutation line of -preserve_seqorder). Exit code 1 = the program stops with an error (no matrix or unreadable matrix, no sequence file, no sequences, conflicting flags), 255 = called without arguments (usage). The alignment progress lines of standard error are not reproduced (the program prints uninitialised %ld arguments there); the '...Read N sequences...' / '...Wrote N sequences...' lines are. -read_msa / -read_msa2 / -read_msa_list / -subset / -remove are not ported (NotSupportedException, registered). VERIFICATION: exit code, standard output and the three output files equal the installed Debian poa 2.0+20060928 on 359 / 359 committed command lines.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/PoaProgram.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/PoaProgram.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md) (section POA 2.0)
- Oracle: `scripts/alignment/generate_poa_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/poa_oracle.tsv`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `arguments` | array<string> | Yes | The command-line arguments of poa, one per element (for example ["-read_fasta", "in.fa", "-clustal", "out.aln", "blosum80.mat"]). |
| `files` | object<string,string> | No | The text of the files the arguments name: file name -> text (FASTA file, matrix file, pair-score file). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `exitCode` | integer | 0 = success, 1 = the program stops with an error, 255 = no arguments (usage). |
| `stdout` | string | Standard output: -v gap tables, -printmatrix table, the permutation line of -preserve_seqorder. |
| `stderr` | string | Error messages and the "...Read N sequences..." / "...Wrote N sequences..." lines (the alignment progress lines are not reproduced). |
| `files` | object<string,string> | The files the program wrote: -pir, -clustal, -po. |

## Errors

| Code | Message |
|------|---------|
| 1001 | arguments is null |
| 1002 | -read_msa / -read_msa2 / -read_msa_list (accepted by the program, not ported): NotSupportedException |

## Examples

### Example 1: CLUSTAL output with the gap penalty tables (poa -read_fasta in.fa -clustal out.aln -v m.mat)

**Tool Call:**
```json
{
  "tool": "poa_command",
  "arguments": {
    "arguments": [
      "-read_fasta",
      "in.fa",
      "-clustal",
      "out.aln",
      "-v",
      "m.mat"
    ],
    "files": {
      "in.fa": ">s1 first\nARNDAARNDA\n>s2\nARDDARNDA\n>s3\nARNDNARNDA\n",
      "m.mat": "# demo\nGAP-PENALTIES=12 6 6\n A R N D\nA 7 -3 -3 -3\nR -3 9 -1 -3\nN -3 -1 9 2\nD -3 -3 2 10\n"
    }
  }
}
```

**Response:**
```json
{
  "exitCode": 0,
  "stdout": "X-Gap Penalties (Open, Aff1, Aff2; LTrunc, LDecay): 12 6 6 16 0\nX-Gap Penalties (0, 1, 2, ...): 12 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 ... \nY-Gap Penalties (Open, Aff1, Aff2; LTrunc, LDecay): 12 6 6 16 0\nY-Gap Penalties (0, 1, 2, ...): 12 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 ... \n",
  "stderr": "...Read 3 sequences from sequence file in.fa...\n...Wrote 3 sequences to CLUSTAL file out.aln...\n",
  "files": {
    "out.aln": "CLUSTAL W (1.74) multiple sequence alignment\n\n\ns1                                  ARNDAARNDA\ns2                                  ARDDA-RNDA\ns3                                  ARNDNARNDA\n\n\n"
  }
}
```

### Example 2: PIR output and the score matrix of three symbols

**Tool Call:**
```json
{
  "tool": "poa_command",
  "arguments": {
    "arguments": [
      "-printmatrix",
      "ARN",
      "-read_fasta",
      "in.fa",
      "-pir",
      "out.pir",
      "m.mat"
    ],
    "files": {
      "in.fa": ">s1 first\nARNDAARNDA\n>s2\nARDDARNDA\n>s3\nARNDNARNDA\n",
      "m.mat": "# demo\nGAP-PENALTIES=12 6 6\n A R N D\nA 7 -3 -3 -3\nR -3 9 -1 -3\nN -3 -1 9 2\nD -3 -3 2 10\n"
    }
  }
}
```

**Response:**
```json
{
  "exitCode": 0,
  "stdout": "   A  R  N\nA  7 -3 -3\nR -3  9 -1\nN -3 -1  9\n",
  "stderr": "...Read 3 sequences from sequence file in.fa...\n...Wrote 3 sequences to FASTA-PIR file out.pir...\n",
  "files": {
    "out.pir": ">s1 first\nARNDAARNDA\n>s2 untitled\nARDDA.RNDA\n>s3 untitled\nARNDNARNDA\n"
  }
}
```

## See Also

- [multiple_align_poa](multiple_align_poa.md)
- [probcons_command](probcons_command.md)
- [kalign_command](kalign_command.md)
