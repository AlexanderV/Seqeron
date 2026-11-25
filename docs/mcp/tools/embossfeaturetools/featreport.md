# featreport

EMBOSS 6.6.0 'featreport': the first feature table of a feature file reported for a sequence in a report format.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `featreport` |
| **Method ID** | `EmbossFeatreportProgram.Featreport` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 featreport: the first feature table of a feature file (embl, genbank, refseq, refseqp, swiss, pir, gff2, gff3, dasgff, or any: the library tries every reader) is reported for a sequence in a report format (reportFormat: gff default, gff3, embl, genbank, pir, swiss, dasgff, debug, listfile, dbmotif, diffseq, draw, restrict, excel, feattable, motif, nametable, regions, seqtable, simple, srs, summary, table, tagseq, or the report qualifiers as typed: -rformat table -raccshow -rdesshow -rscoreshow N -rstrandshow N -rusashow -rmaxall 3 -rmaxseq 2). The formats gff, gff3, embl, genbank, pir, swiss and dasgff are written by the feature writers of the library (joins, tags, flags); the others list the start, end, strand and score of the features. The features come with an optional range (featureBegin, featureEnd). Equals the installed program on 1 000 runs (every feature format, every report format, ranges, damaged files). A5-8c3h, B09 F262. The report file name qualifiers of reportFormat (-rname -rextension -rdirectory: the default file name is the name of the sequence + .featreport, a directory is put in front of a name that has none, a directory that does not exist stops the program) and the qualifiers of the features input in qualifiers (-fformat2 -fopenfile2 -iquery2 -ioffset2 -fask2 -fbegin2 -fend2 -freverse2 -fcircular2; only -fbegin2 / -fend2 change the report, as in the program) are applied as the program applies them; the file name is returned in reportFile. A5-8c3h, B09 F278.

## Core Documentation Reference

- Source: [EmbossFeatreportProgram.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Emboss.FeatureTools/EmbossFeatreportProgram.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_featreport_oracle.py](../../../../scripts/alignment/generate_featreport_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file (any input format of seqret), or a sequence as letters. |
| `features` | string | Yes | The text of the feature file (see the formats above). |
| `sequenceFormat` | string | No | The input format of the sequence (-sformat); default: the program finds the format. |
| `reportFormat` | string | No | The report format name or the report qualifiers as typed (see above); default gff. (default: gff) |
| `qualifiers` | string | No | Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sformat1 F, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes. |
| `featureBegin` | integer | No | The first position of the features to report (-fbegin2); default 1. (default: 1) |
| `featureEnd` | integer | No | The last position of the features to report (-fend2); default the length of the table. |
| `files` | object | No | More files by name (file name -> text) that sequenceValue and featuresValue may name. |
| `sequenceValue` | string | No | The USA of the sequence input (-sequence): [format::]file[begin:end:r], file:id, file:field:value, file{a\|b}, @list, list::list, asis::TEXT, db:id; the text of sequence is the file 'sequence.dat'. |
| `featuresValue` | string | No | The value of -features: [format::]file[begin:end:r]; the text of features is the file 'features.dat'. |
| `definitions` | string | No | The text of an EMBOSS configuration file (.embossrc) that defines databases for sequenceValue: DBNAME definitions (methods direct, emboss and emblcd), RES / RESOURCE, SERVER, INCLUDE (the included file is an entry of files), IF / IFDEF / ELSE / ENDIF and SET blocks; the files of a database, and for the methods emboss (dbxfasta / dbxflat index files NAME.ent, .xid, .pxid ...) and emblcd (dbifasta / dbiflat files division.lkp, entrynam.idx, acnum.trg ...) its index files as latin-1 text, are entries of files by path. |
| `outputFile` | string | No | The name of the report file (-outfile); default: standard output, or - when reportFormat has -rname / -rextension / -rdirectory - the file named after the sequence (name.featreport) in that directory; it is returned in reportFile. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `report` | string | The text of the report file. |
| `messages` | string[] | What the library printed (warnings and errors). |

## Errors

| Code | Message |
|------|---------|
| 1001 | The sequence cannot be empty |
| 1001 | The features cannot be empty |
| 1001 | Bad value for '-features' and no more retries (the text is no feature file) |
| 1001 | Report format ... needs 7 tags (diffseq) / not supported for ... sequences |
| 1001 | An input on which the reference program never stops or dies by a fault is refused |

## Examples

### Example 1: The features of an EMBL entry as an Excel table

**User Prompt:**
> Report the features of this EMBL entry in the excel report format.

**Expected Tool Call:**
```json
{"tool": "featreport", "arguments": {"sequence": "ID   X1; SV 1; linear; genomic DNA; STD; UNC; 60 BP.\nXX\nAC   X1;\nXX\nDE   test\nXX\nFH   Key             Location/Qualifiers\nFH\nFT   mRNA            1..40\nFT   exon            45..55\nSQ   Sequence 60 BP; 15 A; 15 C; 15 G; 15 T; 0 other;\n     acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt        60\n//\n", "features": "ID   X1; SV 1; linear; genomic DNA; STD; UNC; 60 BP.\nXX\nAC   X1;\nXX\nDE   test\nXX\nFH   Key             Location/Qualifiers\nFH\nFT   mRNA            1..40\nFT   exon            45..55\nSQ   Sequence 60 BP; 15 A; 15 C; 15 G; 15 T; 0 other;\n     acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt        60\n//\n", "reportFormat": "excel"}}
```

**Response:**
```json
{"report": "SeqName\tStart\tEnd\tScore\tStrand\nX1\t1\t40\t0.000\t+\nX1\t45\t55\t0.000\t+\n", "messages": []}
```

The excel format has no header; the other formats print the run date and the command line.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
