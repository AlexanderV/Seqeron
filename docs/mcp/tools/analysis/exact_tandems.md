# exact_tandems

MUMmer exact-tandems: the exact tandem repeats of a sequence (repeat-match -t | sort | tandem-repeat.awk) and the script as installed.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `exact_tandems` |
| **Method ID** | `MummerExactTandems.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MUMmer 3.23 exact-tandems: the exact tandem repeats of the first sequence of a FASTA file with a duplicated unit of at least minLength bases. The script runs repeat-match -t -n N, drops its two header lines, sorts the matches by the first and second column (sort -k1n -k2n) and feeds them to tandem-repeat.awk, which writes the table 'Start Extent UnitLen Copies' with one line per repeat region (the first match of every start whose copies touch; extent = start2 + length - start1, unit = start2 - start1, copies = extent / unit with one decimal). By default (installedScript = false) the result is that intended pipeline. The script as installed by the Debian package is broken - its last command is split over two lines (awk -f $ and, on the next line, $scriptdir/tandem-repeat.awk) - so awk looks for a file named $ and fails with status 2 ('awk: fatal: cannot open source file `$' for reading: No such file or directory'); sh -e then stops the script without a table and without removing its temporary file: installedScript = true reproduces that (text = the two progress lines 'Finding matches' and 'Tandem repeats', exitCode 2, filesLeft = $$.tmp.matches). messages = standard error (the lines of repeat-match, then the awk error with installedScript). A minLength below 1 makes repeat-match fail (the table is just the header). An empty file or a file without a '>' line is refused (see repeat_match). Byte for byte equal to the installed script (installedScript = true) and to the pipeline run by hand on 215 + runs each and 1 200 more random runs each.

## Core Documentation Reference

- Source: [MummerExactTandems.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MummerExactTandems.cs)
- Algorithm doc: [Mummer_Package_Programs.md](../../../algorithms/Alignment/Mummer_Package_Programs.md)
- Oracle: [generate_repeatmatch_oracle.py](../../../../scripts/alignment/generate_repeatmatch_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genome` | string | Yes | Genome FASTA text (only the first record is used). |
| `minLength` | integer | No | Minimum length of the duplicated unit (repeat-match -n; default 20). (default: 20) |
| `installedScript` | boolean | No | Reproduce the script as installed (broken: awk -f $), instead of the intended pipeline (default false). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The table of tandem repeats (or the two progress lines of the installed script). |
| `messages` | string | The standard error of the script. |
| `exitCode` | integer | The exit status (0; 2 for the installed script). |
| `filesLeft` | array | The files the script leaves behind ($$.tmp.matches for the installed script). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | repeat-match: no sequence record: the program works on an uninitialised buffer (undefined; valgrind: use of uninitialised values) |

## Examples

### Example 1: A tandem array and a degenerate one

**User Prompt:**
> Find the exact tandem repeats with a unit of at least 6 bases.

**Expected Tool Call:**
```json
{"tool": "exact_tandems", "arguments": {"genome": ">tandem demo\nacactgctcactccaaccccggcccctgagtccgaggagaggatcacagtctggatcaca\ngtctggatcacagtctggatcacagtctggatcacagtctgggtgcttcagagtatgtat\naccactgggtggatcaggatcaggatcaaggatacggcggagggcacg\n", "minLength": 6}}
```

**Response:**
```json
{"text": "   Start   Extent  UnitLen     Copies\n      41       62       12        5.2\n     131       18        6        3.0\n", "messages": "Genome Length = 168   Used 129 internal nodes\n", "exitCode": 0, "filesLeft": []}
```

exact-tandems t.fa 6 (intended pipeline): a unit of 12 bases repeated 5.2 times from position 41 and a unit of 6 bases from 131.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
