# twofeat

EMBOSS 6.6.0 'twofeat': pairs of features of a sequence in a required relation (overlap, distance, sense) reported.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `twofeat` |
| **Method ID** | `EmbossTwofeatProgram.Twofeat` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 twofeat: the features of the sequences of a file (any input format of seqret that carries features) are tested against the criteria of feature A (aSource, aType, aSense 0 / + / -, aMinScore, aMaxScore, aTag, aValue) and of feature B (bSource ... bValue); wildcards * and ?, several patterns separated by blanks, commas, semicolons or bars; the scores are cut to integers) and every pair of an A and a B feature in the required relation is reported: overlap (A any, O overlap, NO none, NW overlap but not within, AW A within B, BW B within A), distance between minRange and maxRange (ignored unless min is less than max; measured from the nearest, left, right or furthest ends: rangeType N / L / R / F), sense (A any, S same, O opposite: the opposite sense stops the program with a fatal error of the library), as a new feature (type typeOut, source twofeat, covering both, tags startA endA startB endB) or, with twoOut, as the two features. The report is in reportFormat (table default; the names of featreport or the report qualifiers as typed). Equals the installed program on 1 000 runs (every input format with features, every option, every report format). A5-8c3h, B09 F263.

## Core Documentation Reference

- Source: [EmbossTwofeatProgram.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/EmbossTwofeatProgram.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_twofeat_oracle.py](../../../../scripts/alignment/generate_twofeat_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file with features (any input format of seqret), or a sequence as letters (no features). |
| `sequenceFormat` | string | No | The input format (-sformat); default: the program finds the format. |
| `aSource` | string | No | Source of feature A (-asource); default *. (default: *) |
| `aType` | string | No | Type of feature A (-atype); default *. (default: *) |
| `aSense` | string | No | Sense of feature A (-asense): 0 any, + forward ('-' is read as any by the program). (default: 0) |
| `aMinScore` | number | No | Minimum score of feature A (-aminscore). (default: 0) |
| `aMaxScore` | number | No | Maximum score of feature A (-amaxscore). (default: 0) |
| `aTag` | string | No | Tag of feature A (-atag); default *. (default: *) |
| `aValue` | string | No | Value of the tags of feature A (-avalue); default *. (default: *) |
| `bSource` | string | No | Source of feature B (-bsource); default *. (default: *) |
| `bType` | string | No | Type of feature B (-btype); default *. (default: *) |
| `bSense` | string | No | Sense of feature B (-bsense). (default: 0) |
| `bMinScore` | number | No | Minimum score of feature B (-bminscore). (default: 0) |
| `bMaxScore` | number | No | Maximum score of feature B (-bmaxscore). (default: 0) |
| `bTag` | string | No | Tag of feature B (-btag); default *. (default: *) |
| `bValue` | string | No | Value of the tags of feature B (-bvalue); default *. (default: *) |
| `overlap` | string | No | Overlap (-overlap): A any, O overlap, NO none, NW overlap but not within, AW A within B, BW B within A. (default: A) |
| `minRange` | integer | No | Minimum distance between the features (-minrange). (default: 0) |
| `maxRange` | integer | No | Maximum distance between the features (-maxrange). (default: 0) |
| `rangeType` | string | No | Where the distance is measured (-rangetype): N nearest ends, L left ends, R right ends, F furthest ends. (default: N) |
| `sense` | string | No | Sense of the pair (-sense): A any, S same, O opposite (fatal in the program). (default: A) |
| `order` | string | No | Order of the features (-order): accepted, not used by the program. (default: A) |
| `twoOut` | boolean | No | Report the two features of a pair instead of one feature that covers them (-twoout). (default: false) |
| `typeOut` | string | No | The type of the new feature (-typeout); default misc_feature. (default: misc_feature) |
| `reportFormat` | string | No | The report format name or the report qualifiers as typed (see featreport); default table. (default: table) |
| `qualifiers` | string | No | Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sformat1 F, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes. |
| `files` | object | No | More files by name (file name -> text) that sequenceValue may name. |
| `sequenceValue` | string | No | The USA of the sequence input (-sequence): [format::]file[begin:end:r], file:id, file:field:value, file{a\|b}, @list, list::list, asis::TEXT, db:id; the text of sequence is the file 'sequence.dat'. |
| `definitions` | string | No | The text of an EMBOSS configuration file (.embossrc) that defines databases for sequenceValue. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `report` | string | The text of the report file. |
| `messages` | string[] | What the library printed (warnings and errors). |

## Errors

| Code | Message |
|------|---------|
| 1001 | The sequence cannot be empty |
| 1001 | Unable to read sequence ... (no sequence in the file) |
| 1001 | Unknown Overlap / rangetype / sense code |
| 1001 | Report format ... needs 7 tags (diffseq) |
| 1001 | An input on which the reference program never stops or dies by a fault, and the opposite sense on opposite strands (a fatal error of the library), are refused |

## Examples

### Example 1: A mRNA and an exon of an EMBL entry

**User Prompt:**
> Report the pairs of a mRNA feature and an exon feature of this EMBL entry in the excel format.

**Expected Tool Call:**
```json
{"tool": "twofeat", "arguments": {"sequence": "ID   X1; SV 1; linear; genomic DNA; STD; UNC; 60 BP.\nXX\nAC   X1;\nXX\nDE   test\nXX\nFH   Key             Location/Qualifiers\nFH\nFT   mRNA            1..40\nFT   exon            45..55\nSQ   Sequence 60 BP; 15 A; 15 C; 15 G; 15 T; 0 other;\n     acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt        60\n//\n", "aType": "mRNA", "bType": "exon", "reportFormat": "excel"}}
```

**Response:**
```json
{"report": "SeqName\tStart\tEnd\tScore\tStrand\tstarta\tenda\tstartb\tendb\nX1\t1\t55\t0.000\t+\t1\t40\t45\t55\n", "messages": []}
```

The feature that covers a pair runs from the first start to the last end and carries the positions of both features as tags.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
