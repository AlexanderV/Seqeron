# repeat_match

MUMmer repeat-match: the maximal exact repeats of a sequence on both strands, in the order of the program's suffix tree, byte for byte.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `repeat_match` |
| **Method ID** | `MummerRepeatMatch.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MUMmer 3.23 repeat-match (Delcher): the maximal exact repeats of the first sequence of a FASTA file that are at least minLength long (-n, default 20), within the forward strand and, unless forwardOnly (-f), between the forward and the reverse-complement strand. Output: 'Long Exact Matches:' and a line per match 'start1 start2[r] length' (a line with the letter r after the second position matches the reverse-complement strand: start2 is read on the original sequence). tandemOnly (-t, implies forwardOnly) lists only the repeats whose two copies touch (start1 + length >= start2). The matches come in the order of the program's own suffix tree (a suffix tree of the sequence and its reverse complement in which a new child is put in front of its siblings), reproduced by the same construction. exhaustive (-E) uses the program's quadratic search instead (a different order); verbose (-V 1) prints the sizes of the program's structures, -V 2 the construction and the traversal of the tree. messages = standard error ('Genome Length = n   Used m internal nodes', the warnings of the FASTA reader). A genome that is an exact reverse-complement palindrome stops the program (exitCode 255, '*** Genome is exact palindrome ***' on standard output). An empty file, a file without a '>' line (the program works on an uninitialised buffer: valgrind use of uninitialised values; the output changes under MALLOC_PERTURB_) and a header whose first word has 500 or more characters (buffer overflow, abort) are refused with an error. Byte for byte equal to the installed repeat-match on 370 + runs and 3 300 more random runs (every option, both strands, N / IUPAC characters, genomes of 30 kb).

## Core Documentation Reference

- Source: [MummerRepeatMatch.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MummerRepeatMatch.cs)
- Algorithm doc: [Mummer_Package_Programs.md](../../../algorithms/Alignment/Mummer_Package_Programs.md)
- Oracle: [generate_repeatmatch_oracle.py](../../../../scripts/alignment/generate_repeatmatch_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genome` | string | Yes | Genome FASTA text (only the first record is used). |
| `minLength` | integer | No | -n: minimum length of an exact match (default 20; at least 1). |
| `forwardOnly` | boolean | No | -f: forward strand only (default false). |
| `tandemOnly` | boolean | No | -t: only tandem repeats; implies forwardOnly (default false). |
| `exhaustive` | boolean | No | -E: exhaustive (quadratic) search instead of the suffix tree (default false). |
| `verbose` | integer | No | -V: verbose level (default 0; 1 prints the structure sizes, 2 the tree construction and traversal). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | What repeat-match writes to standard output. |
| `messages` | string | What the program writes to standard error. |
| `exitCode` | integer | The exit status of the program (0, or 255 when it stops). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | ERROR:  Illegal min match length "0" |
| 1001 | repeat-match: no sequence record: the program works on an uninitialised buffer (undefined; valgrind: use of uninitialised values) |

## Examples

### Example 1: A repeat and its reverse complement

**User Prompt:**
> Find the exact repeats of at least 3 bases of this sequence.

**Expected Tool Call:**
```json
{"tool": "repeat_match", "arguments": {"genome": ">pal\nacgtagcta\n", "minLength": 3}}
```

**Response:**
```json
{"text": "Long Exact Matches:\n   Start1     Start2    Length\n        4          9r        6\n        1          4r        4\n", "messages": "Genome Length = 9   Used 14 internal nodes\n", "exitCode": 0}
```

repeat-match -n 3 pal.fa: two matches between the forward and the reverse-complement strand (letter r): 6 bases from position 4 and 4 bases from position 1.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
