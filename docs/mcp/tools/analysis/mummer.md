# mummer

MUMmer 3.23 mummer (maxmat3): the maximal exact matches between a reference and query files, in the program's order.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `mummer` |
| **Method ID** | `MummerProgram.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MUMmer 3.23 mummer (Kurtz' maxmat3): the maximal exact matches between a reference file and one or more query files. matchMode 'mumreference' (default, -mumreference / -mumcand) reports matches unique in the reference, 'mum' (-mum) matches unique in both, 'maxmatch' (-maxmatch) every maximal match; minLength = -l (default 20); nucleotidesOnly (-n) matches only a, c, g, t (any case) - without it every character is a letter of the alphabet (also protein text), both (-b) adds the reverse-complement strand of the query, reverseOnly (-r) only that strand. Output per query record: a '> id' line ('> id Reverse' for the reverse strand, '  Len = n' appended with showSequenceLengths = -L) and one line per match: reference position, query position, length (%8lu columns; with several reference records or fourColumn = -F the reference id is the first column; for reverse matches the query position is a position in the reverse-complemented query, or with showReversePositions = -c relative to the original query), followed by the matching text with showString (-s). Matches come in the order the suffix tree of the program enumerates them (per query position, the reference occurrences in lexicographic order of their suffixes; -mum by reference position). The reference records are joined by a separator no match can span; text before the first '>' belongs to the first record; the residues are lower-cased. The file read errors of the program (empty file, no header, only empty records) and its option conflicts (-b with -r, -c without -b / -r, a non-positive -l, no query file, more than 32 query files) are errors. Byte for byte equal to the installed mummer on 3 200 random runs (every mode, every option, multi-record files, N / IUPAC / odd characters, CR LF).

## Core Documentation Reference

- Source: [MummerProgram.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MummerProgram.cs)
- Algorithm doc: [Mummer_Package_Programs.md](../../../algorithms/Alignment/Mummer_Package_Programs.md)
- Oracle: [generate_mummer_cli_oracle.py](../../../../scripts/alignment/generate_mummer_cli_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reference` | string | Yes | Reference FASTA text (multi-record). |
| `queries` | array<string> | Yes | Query FASTA texts (at most 32 files); each is matched in turn. |
| `matchMode` | string | No | Matching algorithm: 'mumreference' (default, -mumreference / -mumcand), 'mum' (-mum) or 'maxmatch' (-maxmatch). |
| `nucleotidesOnly` | boolean | No | -n: match only the characters a, c, g, t (default false). |
| `minLength` | integer | No | -l: minimum length of a match (default 20). |
| `both` | boolean | No | -b: compute forward and reverse-complement matches (default false). |
| `reverseOnly` | boolean | No | -r: compute only reverse-complement matches (default false). |
| `showString` | boolean | No | -s: show the matching substrings (default false). |
| `showReversePositions` | boolean | No | -c: report the query position of a reverse match relative to the original query (needs both or reverseOnly; default false). (default: false) |
| `fourColumn` | boolean | No | -F: force the 4 column format also for a single reference record (default false). |
| `showSequenceLengths` | boolean | No | -L: show the length of the query sequence on the header line (default false). |
| `referenceName` | string | No | The reference file name shown in the messages (default reference.fa). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | What mummer writes to standard output. |
| `messages` | string | The deterministic lines of standard error ('# reading input file ...'; the timing lines are left out). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | At least one query file is required |
| 1001 | matchMode must be 'mumreference', 'mum' or 'maxmatch' |
| 1001 | mummer: empty sequence in multiple fasta file |
| 1001 | mummer: option -b and option -r exclude each other |

## Examples

### Example 1: Maximal unique matches on both strands

**User Prompt:**
> Find the maximal unique matches of length >= 8 on both strands.

**Expected Tool Call:**
```json
{"tool": "mummer", "arguments": {"reference": ">r\nacgtacgtacgtaaccggtt\n", "queries": [">q\nacgtacgtaacctacgtacgtacgtaaccggtta\n"], "matchMode": "mum", "minLength": 8, "both": true, "showString": true}}
```

**Response:**
```json
{"text": "> q\n       1        14        20\nacgtacgtacgtaaccggtt\n> q Reverse\n       1        10        13\nacgtacgtacgta\n      12         1         9\ntaaccggtt\n", "messages": "# reading input file \"reference.fa\" of length 20\n# construct suffix tree for sequence of length 20\n# (maximum reference length is 536870908)\n# (maximum query length is 4294967295)\n# reading input file \"query.fa\" of length 34\n# matching query-file \"query.fa\"\n# against subject-file \"reference.fa\"\n"}
```

mummer -mum -l 8 -b -s r.fa q.fa: the whole reference matches the query at 14; the reverse-complement strand matches twice (the positions are in the reverse-complemented query).

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
