# mafft_command_last

`mafft --lastpair` / `--lastmultipair` of MAFFT 7.505 with the programs of LAST 1256 (separately licensed, GPL-3.0-or-later).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Last |
| **Tool Name** | `mafft_command_last` |
| **Method ID** | `Seqeron.Genomics.Alignment.MafftProgram.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The script mafft of MAFFT 7.505 as a function of its command line and the text of the files it names (MafftProgram.Run), with the programs lastdb and lastal of LAST 1256 (ported, GPL) installed: --lastpair aligns every pair of nucleotide sequences with lastal (Last-INS-1). pairlocalalign -R of MAFFT 7.505 overflows its heap when three or more sequences are related (valgrind: invalid writes behind the rows of its table of local homologies; the run then crashes): such runs are refused with an error, the runs with two sequences and the runs without hits are answered. --lastmultipair stops with the messages of the script (it needs --addfragments, and dndpre does not know -Q).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/MafftProgram.cs#L55](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MafftProgram.cs#L55)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `arguments` | array<string> | Yes | The command-line arguments of mafft, one per element, the input file LAST (for example ["--lastpair", "--clustalout", "in.fa"]). |
| `files` | object | No | The text of the files the arguments name: file name -> content (the input FASTA, --addfragments / --add / --merge / --seed files ...). |
| `stdin` | string | No | Text of standard input (read when the input file is '-'). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `exitCode` | integer | Exit status of the script |
| `stdout` | string | Standard output (the alignment) |
| `stderr` | string | Standard error (the messages of the programs and of the script) |
| `files` | object | The files the run wrote (tree and distance files ...) |

## Errors

| Code | Message |
|------|---------|
| 1001 | arguments must contain the command-line arguments of mafft (the input file last). |
| — | `NotSupportedException`: the run is one in which `pairlocalalign -R` of MAFFT 7.505 overflows its heap (three or more related sequences; every such run of the oracle crashed, valgrind: "Invalid write of size 4") |

## Examples

### Example 1: Two related sequences

**Tool Call:**
```json
{
  "tool": "mafft_command_last",
  "arguments": {
    "arguments": ["--lastpair", "in.fa"],
    "files": { "in.fa": ">a\nACGT...\n>b\nACGT...\n" }
  }
}
```

**Response:** `exitCode` 0, the alignment on `stdout`, the messages of `pairlocalalign`, `lastal` and `tbfast` on `stderr`; the text ends with "Pairwise alignments were computed by LAST".

## Worked Example

`--lastpair` with two sequences: `lastdb _dbd _dbd` makes the database of the sequences, for every query k `lastdb _db<k>d` the database of the k sequences before it, `lastal -m max(10, 3k) -e 5000 -f 0 -s 1 -p _scoringmatrixforlast -a 3599 -b 179 _db<k>d _q<k>` the alignments (TAB format); the first alignment of a pair gives its score, the blocks of its alignment become the pair library `hat3`, the distances `hat2`; `tbfast` aligns progressively (Last-INS-1).

## See Also

- [`mafft_command`](../alignment/mafft_command.md) - the same script without LAST (the MIT server): `--lastpair` stops there with the "Install LAST" text of the script.
