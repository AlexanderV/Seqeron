# extractfeat

EMBOSS 6.6.0 'extractfeat': the features of the sequences of a file (any input format of seqret with features) that pass the source, type, sense, score, tag and value filters, as sequences with the sequence before and after them.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `extractfeat` |
| **Method ID** | `EmbossExtractfeatProgram.Extractfeat` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 extractfeat: the features of the sequences of a file (any input format of seqret that carries features: embl, genbank, refseq, refseqp, swiss, nbrf, gff2, gff3) that pass the filters (source, type, sense, score, tag, value; wildcards * and ?, several patterns separated by blanks, commas, semicolons or bars) are extracted as sequences, named SEQUENCE_START_END (with featInName also _TYPE) and described by [type] and the tags of describe, optionally with the sequence before and after the feature (before / after: negative values count from the other end of the feature; the missing sequence at the ends is padded with N or X) and with the parts of a join as one sequence (join). The sequences are returned in FASTA or in the output format of outputFormat (fasta unknown pearson, gcg gcg8, nbrf pir, codata, strider, staden, text plain raw, ig iguspto, gde, fitch, acedb, mase). A part of a join that names another entry (J00194.1:100..103) is read from the same file or database. Equals the installed program on 1 000 runs (every input format with features, every filter, before / after / join, 22 output formats, damaged files). A5-8c3h, B09 F258.

## Core Documentation Reference

- Source: [EmbossExtractfeatProgram.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/EmbossExtractfeatProgram.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_extractfeat_oracle.py](../../../../scripts/alignment/generate_extractfeat_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file with features (see the formats above), or a sequence as letters (no features: nothing is extracted). |
| `sequenceFormat` | string | No | The input format (-sformat); default: the program finds the format. |
| `before` | integer | No | Sequence before each feature (-before): bases / residues; negative: the start is that far before the END of the feature. (default: 0) |
| `after` | integer | No | Sequence after each feature (-after); negative: the end is that far after the START of the feature. (default: 0) |
| `join` | boolean | No | Write the parts of a join (CDS, mRNA, ...) as one sequence (-join). (default: false) |
| `featInName` | boolean | No | Add the feature type to the sequence name (-featinname). (default: false) |
| `describe` | string | No | Tag names to add with their values to the description (-describe): names separated by blanks, commas, semicolons or bars; wildcards. (default: ) |
| `source` | string | No | Sources of the features to extract (-source); default *. (default: *) |
| `type` | string | No | Types of the features to extract (-type); default *. (default: *) |
| `sense` | integer | No | Sense of the features (-sense): 0 any, 1 forward, -1 reverse. (default: 0) |
| `minScore` | number | No | Minimum score (-minscore); with maxScore both 0 the score is not tested. (default: 0) |
| `maxScore` | number | No | Maximum score (-maxscore). (default: 0) |
| `tag` | string | No | Tags a feature must have (-tag); default *. (default: *) |
| `value` | string | No | Values the tag must have (-value); default *. (default: *) |
| `outputFormat` | string | No | The sequence output format (-osformat); default fasta. |
| `qualifiers` | string | No | Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sformat1 F, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes. |
| `files` | object | No | More files by name (file name -> text) that sequenceValue may name: list files, other sequence files, the data files of the databases of definitions. |
| `sequenceValue` | string | No | The USA of the sequence input (-sequence): [format::]file[begin:end:r], file:id, file:field:value, file{a\|b}, @list, list::list, asis::TEXT, db:id; the text of sequence is the file 'sequence.dat'. |
| `definitions` | string | No | The text of an EMBOSS configuration file (.embossrc) that defines databases for sequenceValue: DBNAME name [ type: N format: fasta method: direct directory: "dir" filename: "*.fa" ], SET name value. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sequenceText` | string | The text of the sequence output file (FASTA or the output format). |
| `sequences` | integer | The number of sequences written. |
| `messages` | string[] | What the library printed (warnings and errors). |

## Errors

| Code | Message |
|------|---------|
| 1001 | The sequence cannot be empty |
| 1001 | Unable to read sequence ... (no sequence in the file) |
| 1001 | An input on which the reference program never stops or dies by a fault, and an output format this port does not write, are refused |

## Examples

### Example 1: The mRNA and exon features of an EMBL entry

**User Prompt:**
> Extract the mRNA and exon features of this EMBL entry, with the feature type in the sequence name.

**Expected Tool Call:**
```json
{"tool": "extractfeat", "arguments": {"sequence": "ID   X1; SV 1; linear; genomic DNA; STD; UNC; 60 BP.\nXX\nAC   X1;\nXX\nDE   test\nXX\nFH   Key             Location/Qualifiers\nFH\nFT   mRNA            1..40\nFT   exon            45..55\nSQ   Sequence 60 BP; 15 A; 15 C; 15 G; 15 T; 0 other;\n     acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt        60\n//\n", "type": "mRNA exon", "featInName": true}}
```

**Response:**
```json
{"sequenceText": ">X1_1_40_mRNA [mRNA] test\nacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt\n>X1_45_55_exon [exon] test\nacgtacgtacg\n", "sequences": 2, "messages": []}
```

The name of a feature sequence is the sequence name with the first and last position of the feature; the description is the type in brackets, the tags of describe and the description of the sequence.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
