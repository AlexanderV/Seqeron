# mummer_annotate

MUMmer mummer-annotate: the edit-distance alignment of every gap of a gaps file, and witherrors.gaps.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `mummer_annotate` |
| **Method ID** | `MummerAnnotate.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MUMmer 3.23 mummer-annotate (annotate): reads the output of gaps (the first sequence is named on its '>' lines - supply it in files by that name - and the word 'reverse' there reverse-complements the query), repeats the gap file and, under every gap line, prints the edit-distance alignment of the stretch between the two matches (10 bases of each flanking match included): '    Errors = n' and the rows 'T:' (first sequence) / 'S:' (query) / mismatch marks, 60 columns per block; a gap longer than 10 000 prints '*** Too long ***'. withErrors is the file witherrors.gaps the program writes: the gap lines with the number of errors appended. A window that reaches before the first base of a sequence reads memory in front of the program's buffer (valgrind: invalid read) and is refused with an error; beyond the end the buffer holds zeros. Byte for byte equal to the installed program on 230 + runs.

## Core Documentation Reference

- Source: [MummerAnnotate.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MummerAnnotate.cs)
- Algorithm doc: [Mummer_Package_Programs.md](../../../algorithms/Alignment/Mummer_Package_Programs.md)
- Oracle: [generate_mummer_gaps_oracle.py](../../../../scripts/alignment/generate_mummer_gaps_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `gaps` | string | Yes | The gap file text (the output of mummer_gaps). |
| `dataFasta` | string | Yes | The data file (the query): its first record is used. |
| `files` | object | Yes | The FASTA texts of the files named on the '>' lines of the gap file, by name. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | What mummer-annotate writes to standard output. |
| `withErrors` | string | The file witherrors.gaps. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Gap file text cannot be null or empty |
| 1001 | Sequence cannot be null or empty |
| 1001 | mummer-annotate: a gap window reaches before the first base of a sequence; ... (undefined behaviour). |

## Examples

### Example 1: Align the gap between two matches

**User Prompt:**
> Align the gap of this gaps file.

**Expected Tool Call:**
```json
{"tool": "mummer_annotate", "arguments": {"gaps": "> ref.fa  Consistent matches\n       6        1     65    none      -      -\n      75       70     66    none      4      4\n> ref.fa  Other matches\n", "dataFasta": ">qry\nAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTGGGGAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA\n", "files": {"ref.fa": ">ref\nGCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA\n"}}}
```

**Response:**
```json
{"text": "> ref.fa  Consistent matches\n       6        1     65    none      -      -\n      75       70     66    none      4      4\n    Errors = 3\nT:  cttaagggttaagtaagtgtgatg\nS:  cttaagggttggggaagtgtgatg\n              ^^ ^          \n> ref.fa  Other matches\n", "withErrors": "> ref.fa  Consistent matches\n       6        1     65    none      -      -\n      75       70     66    none      4      4       3\n> ref.fa  Other matches\n"}
```

mummer-annotate g.gaps qry.fa: the gap between the two matches holds three substitutions.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
