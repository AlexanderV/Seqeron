# minimap2_map

minimap2 2.26 as a program: mapping, base-level alignment (-c / -a), paired reads and the other options.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `minimap2_map` |
| **Method ID** | `Minimap2Program.Run` |
| **Version** | 1.0.0 |
| **Stability** | stable |

## Description

minimap2 2.26 as a program: the whole standard output (PAF or SAM) of `minimap2 [options] ref.fa query.fa [query2.fa]` on the given sequences (every option of the mapping and alignment stages, parsed with the program's getopt rules, `-x` applied first and cumulatively, the option checks of mm_check_opt with the program's exit status and error text). reference = the FASTA text of the target sequences (named ref.fa on the command line), queries = the FASTA / FASTQ text of one query file (query.fa; fragment mode = interleaved reads with the same name, /1 /2 suffixes ignored) or of two (query1.fa, query2.fa: the paired reads of `-x sr`), options = the command-line words before the files, e.g. ["-x", "sr", "-a"] or ["-c", "--cs", "--MD", "-N", "10"]. Base-level alignment (`-c`, `-a`, `--cs`, `--MD`: align.c, ksw2 gap filling between the chained anchors, left / right extension, Z-drop, inversion rescue, spliced alignment with `-x splice` / `-u` / `--junc-bed`, the CIGAR `cg`, `cs` short / long, `MD`, `--eqx`, `-Y`, `-L`, `--secondary-seq`, tags NM ms AS nn ts de zd), the SAM writer (flags, MAPQ, SA, mate fields, @HD / @SQ / @RG / @PG header), the multi-segment paired-read mode (segment ids in the anchors, mm_select_sub_multi, mm_seg_gen, mm_pair with pe_bonus, proper pairs), the index in parts (`-I`) and the `--split-prefix` merge, `--alt` / `--alt-drop`, `--qstrand`, `-X` / `-P`, `--paf-no-hit`, `-o FILE` (the output is returned in files) and the debug-only options (`-t -K --print-seeds --print-chains --print-aln-seq --no-kalloc --cap-kalloc --bucket-bits`, no effect on the output) are ported from the 2.26 source and compared with the installed program as text: 2 700 base-level PAF / SAM cases, 2 100 paired cases, 2 401 option / error cases (exit status and standard error included), all equal; the ksw2 kernels ksw_extz2_sse / ksw_extd2_sse / ksw_exts2_sse / ksw_ll_i16 equal the C functions compiled with SSE4.1 on 2 000 random inputs. The program's own log lines (`[M::...]`) and `-d` (index dump) are not reproduced; -t has no effect; a k above 28 or a window of 256 or more is refused (the C program aborts).

## Core Documentation Reference

- Seqeron.Genomics.Alignment/Minimap2Program.cs (the program: options, readers, pipeline), Minimap2Chaining.Align.cs (align.c), Minimap2Chaining.Pe.cs (pe.c), Minimap2Chaining.Output.cs (PAF / SAM writers)
- Algorithm documentation: docs/algorithms/Alignment/Minimap2_Base_Alignment.md

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reference` | string | Yes | The FASTA text of the target (reference) sequences. |
| `queries` | array<string> | Yes | The FASTA / FASTQ text of one query file, or of two query files (paired reads, `-x sr`). |
| `options` | array<string> | No | The command-line options before the files (e.g. ["-x", "map-ont", "-c", "--cs"]); omit for `minimap2 ref.fa query.fa`. |
| `junctionBed` | string | No | The text of a BED file of annotated junctions for `--junc-bed` (BED12 blocks or BED6 intervals), or null. |
| `altList` | string | No | The text of an ALT contig list for `--alt` (the first word of each line is a target name), or null. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `output` | string | The standard output of the program: PAF or SAM text. |
| `standardError` | string | The error and warning lines of the program (the [M::...] log is not reproduced). |
| `exitCode` | integer | 0, or 1 when the program stops with an error (invalid option, missing file). |
| `files` | object | The file written by `-o FILE` (name -> text); the standard output is then empty. |

## Errors

| Code | Message |
|------|---------|
| 1001 | A reference is required |
| 1001 | One query file (or two for paired reads) is required |
| 1001 | A query file is empty |

## Examples

### Example 1: A query mapped with the base-level alignment (-c)

**User Prompt:**
> Map the query against the reference with the base-level alignment.

**Expected Tool Call:**
```json
{"tool": "minimap2_map", "arguments": {"reference": ">t1\nACGTTGCAGG...\n", "queries": [">q1\nCTAAGCTTCG...\n"], "options": ["-c"]}}
```

**Response:**
```json
{"output": "q1\t97\t0\t97\t+\tt1\t127\t30\t127\t97\t97\t0\t...", "standardError": "", "exitCode": 0}
```

`minimap2 -c ref.fa query.fa` prints the PAF line with `cg:Z:` and the alignment tags (abbreviated here).

## See Also

- [minimap2_chains](minimap2_chains.md) - the seeding and chaining stage on its own
- [ksw2_extend](ksw2_extend.md) - the ksw2 kernels (extz2, extd2, exts2, ll) the base-level alignment uses
