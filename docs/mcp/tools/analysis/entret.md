# entret

EMBOSS 6.6.0 'entret': the full text of the entries of a sequence file (or of the entries a query selects) exactly as the sequence readers store it.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `entret` |
| **Method ID** | `EmbossEntretProgram.Entret` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 entret: the entries of a sequence file (any input format of seqret) are read with the text of every entry stored, and the text of each entry is returned exactly as the file has it (the lines the sequence reader read: the whole alignment for the first sequence of an alignment format, nothing for the later ones). sequenceValue is the USA of the program: a query on the file (file:id, file:acc:value, file:wild*, file{a|b}), a range, a list file or a database defined in definitions (db:id); firstOnly stops after the first entry. The embl, genbank, swiss, gff and the other readers store the lines as the library does (a gff3 entry keeps its version line without the line end). Equals the installed program on 756 runs (every input format, damaged files, queries). A5-8c3h, B09 F257.

## Core Documentation Reference

- Source: [EmbossEntretProgram.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/EmbossEntretProgram.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_entret_oracle.py](../../../../scripts/alignment/generate_entret_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file (any input format of seqret), or a sequence as letters. |
| `sequenceFormat` | string | No | The input format (-sformat); default: the program finds the format. |
| `firstOnly` | boolean | No | Read one entry and stop (-firstonly). (default: false) |
| `qualifiers` | string | No | Sequence input qualifiers as typed after the sequence on the command line: -sformat1 F, -sid1 name, -sdbname1 db, -iquery1 query, -ioffset1 bytes, -sbegin1 -send1 -sreverse1 (the range does not change the text of an entry). |
| `files` | object | No | More files by name (file name -> text) that sequenceValue may name: list files, other sequence files, the data files of the databases of definitions. |
| `sequenceValue` | string | No | The USA of the sequence input (-sequence): [format::]file[begin:end:r], file:id, file:field:value, file{a\|b}, @list, list::list, asis::TEXT, db:id; the text of sequence is the file 'sequence.dat'. |
| `definitions` | string | No | The text of an EMBOSS configuration file (.embossrc) that defines databases for sequenceValue: DBNAME name [ type: N format: fasta method: direct directory: "dir" filename: "*.fa" ] (methods direct, emboss and emblcd), RES / RESOURCE, SERVER, INCLUDE (the included file is an entry of files), IF / IFDEF / ELSE / ENDIF, SET name value; the index files of the methods emboss (dbxfasta / dbxflat: NAME.ent, NAME.xid, NAME.pxid ...) and emblcd (dbifasta / dbiflat: division.lkp, entrynam.idx, acnum.trg ...) are entries of files by path, as latin-1 text. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The text of the output file: the entries as the readers stored them. |
| `entries` | integer | The number of entries written. |
| `messages` | string[] | What the library printed (warnings and errors). |

## Errors

| Code | Message |
|------|---------|
| 1001 | The sequence cannot be empty |
| 1001 | Unable to read sequence ... (no entry matched the query or the file has no sequence) |
| 1001 | An input on which the reference program never stops or dies by a fault is refused |

## Examples

### Example 1: One entry of a two-entry FASTA file

**User Prompt:**
> Retrieve the entry s2 from this FASTA file.

**Expected Tool Call:**
```json
{"tool": "entret", "arguments": {"sequence": ">s1 desc one\nACGTACGTAC\n>s2\nTTGGCCAATT\n", "sequenceValue": "sequence.dat:s2"}}
```

**Response:**
```json
{"text": ">s2\nTTGGCCAATT\n", "entries": 1, "messages": []}
```

The query file:id keeps the entries whose identifier matches; the text is the lines of the entry as read.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
