# showfeat

EMBOSS 6.6.0 'showfeat': the features of the sequences of a file as lines of ASCII graphics under a scale.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `showfeat` |
| **Method ID** | `EmbossShowfeatProgram.Showfeat` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 showfeat: the features of the sequences of a file (any input format of seqret that carries features) as lines of ASCII graphics under a scale: one line per feature (or per part of a join; with joinFeatures the parts of a join on one line, with collapse the features on the same line) with the strand (|---> for forward, <---| for reverse), optionally the strand, source, positions, type and tags with their values. The features are filtered by sourceMatch, typeMatch, tagMatch and valueMatch (wildcards * and ?, several patterns separated by blanks, commas, semicolons or bars; stricttags shows only the tags that match), sorted (sort: start, source, type or nosort; by strand first), the sense is selected by forward / reverse / unknown, and annotation adds ranges with a text ("5-10 foo,20-30 bar") as lines of their own. width is the width of the graphics (below 4 no graphics and no scale), html wraps the output in HTML. Equals the installed program on 1 000 runs (every input format with features, every option, damaged files, the runs on which the program dies by a fault). A5-8c3h, B09 F260.

## Core Documentation Reference

- Source: [EmbossShowfeatProgram.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Emboss.FeatureTools/EmbossShowfeatProgram.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_showfeat_oracle.py](../../../../scripts/alignment/generate_showfeat_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file with features (see the formats above), or a sequence as letters (no features). |
| `sequenceFormat` | string | No | The input format (-sformat); default: the program finds the format. |
| `sourceMatch` | string | No | Sources of the features to show (-sourcematch); default *. (default: *) |
| `typeMatch` | string | No | Types of the features to show (-typematch); default *. (default: *) |
| `tagMatch` | string | No | Tags a feature must have (-tagmatch); default *. (default: *) |
| `valueMatch` | string | No | Values the tag must have (-valuematch); default *. (default: *) |
| `sort` | string | No | Sort order (-sort): start (default), source, type or nosort. |
| `joinFeatures` | boolean | No | Show the parts of a join on one line (-joinfeatures). (default: false) |
| `html` | boolean | No | Output as HTML (-html). (default: false) |
| `id` | boolean | No | Show the name of the sequence (-id); default true. (default: true) |
| `description` | boolean | No | Show the description of the sequence (-description); default true. (default: true) |
| `scale` | boolean | No | Show the scale line (-scale); default true. (default: true) |
| `width` | integer | No | Width of the graphics lines (-width); default 60; below 4 there are no graphics and no scale; negative values are reset to 0 with a warning. (default: 60) |
| `collapse` | boolean | No | Show all features of a type on one line (-collapse). (default: false) |
| `forward` | boolean | No | Show forward features (-forward); default true. (default: true) |
| `reverse` | boolean | No | Show reverse features (-reverse); default true. (default: true) |
| `unknown` | boolean | No | Show features of unknown sense (-unknown); default true. (default: true) |
| `showStrand` | boolean | No | Show the strand of the features (-strand). (default: false) |
| `showSource` | boolean | No | Show the source (-sourceshow). (default: false) |
| `showPosition` | boolean | No | Show the positions (-position). (default: false) |
| `showType` | boolean | No | Show the type (-typeshow); default true. (default: true) |
| `showTags` | boolean | No | Show the tags (-tagshow). (default: false) |
| `showValues` | boolean | No | Show the values of the tags (-valueshow); default true. (default: true) |
| `strictTags` | boolean | No | Show only the tags that match tagMatch and valueMatch (-stricttags). (default: false) |
| `annotation` | string | No | Ranges to show as annotation (-annotation): "5-10 text,20-30 text". (default: ) |
| `qualifiers` | string | No | Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sformat1 F, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes. |
| `files` | object | No | More files by name (file name -> text) that sequenceValue may name. |
| `sequenceValue` | string | No | The USA of the sequence input (-sequence): [format::]file[begin:end:r], file:id, file:field:value, file{a\|b}, @list, list::list, asis::TEXT, db:id; the text of sequence is the file 'sequence.dat'. |
| `definitions` | string | No | The text of an EMBOSS configuration file (.embossrc) that defines databases for sequenceValue: DBNAME definitions (methods direct, emboss and emblcd), RES / RESOURCE, SERVER, INCLUDE (the included file is an entry of files), IF / IFDEF / ELSE / ENDIF and SET blocks; the files of a database, and for the methods emboss (dbxfasta / dbxflat index files NAME.ent, .xid, .pxid ...) and emblcd (dbifasta / dbiflat files division.lkp, entrynam.idx, acnum.trg ...) its index files as latin-1 text, are entries of files by path. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The text of the output file. |
| `sequences` | integer | The number of sequences shown. |
| `messages` | string[] | What the library printed (warnings and errors). |

## Errors

| Code | Message |
|------|---------|
| 1001 | The sequence cannot be empty |
| 1001 | Unable to read sequence ... (no sequence in the file) |
| 1001 | Bad range value ... / Odd integer(s) in range specification (annotation) |
| 1001 | An input on which the reference program never stops or dies by a fault is refused |

## Examples

### Example 1: The features of an EMBL entry, 30 columns wide

**User Prompt:**
> Draw the features of this EMBL entry on 30 columns, joins on one line, with strand and tags.

**Expected Tool Call:**
```json
{"tool": "showfeat", "arguments": {"sequence": "ID   X1; SV 1; linear; genomic DNA; STD; UNC; 60 BP.\nXX\nAC   X1;\nXX\nDE   test\nXX\nFH   Key             Location/Qualifiers\nFH\nFT   CDS             join(3..12,20..29)\nFT                   /gene=\"g1\"\nFT   mRNA            1..40\nSQ   Sequence 60 BP; 15 A; 15 C; 15 G; 15 T; 0 other;\n     acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt        60\n//\n", "width": 30, "joinFeatures": true, "showStrand": true, "showTags": true}}
```

**Response:**
```json
{"text": "X1\ntest\n|============================| 60\n|------------------>           + mRNA \n|---->   |--->                 + CDS  gene=\"g1\"\n\n", "sequences": 1, "messages": []}
```

The scale line shows the width and the length of the sequence; a forward feature starts with | and ends with >, a reverse one starts with < and ends with |.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
