# diffseq

EMBOSS 6.6.0 diffseq: the regions of difference between two nearly identical sequences (336 nucleotide / protein pairs, every -rformat: report and both GFF3 feature files byte for byte). The exact word matches of wordSize symbols (the engine and the minimal non-overlapping set of megamerger) are the regions of identity; the stretches between them are the differences (a difference with length 0 in one sequence is an insertion in the other); with globalDifferences the differences before the first and after the last match are reported too (a start or end that matches by itself is not counted). A nucleotide pair also gets the SNP count and its transitions / transversions. Returns the differences and the program's files: the report in the default format 'diffseq' (7 report tags; reportFormat takes any -rformat or the report qualifiers as typed, e.g. '-rformat table -raccshow'), and the GFF3 files -aoutfeat / -boutfeat. The sequence qualifiers of both inputs (-sbegin1 -send1 -sreverse1 -slower1 -supper1 -snucleotide1 ... and the same with 2, given in reportFormat as typed after the files, e.g. '-sbegin1 20 -send1 200 -sreverse2 -rformat table') are applied as the program applies them: the sequence is cut to the range (positions, lengths and the report numbers refer to the cut sequence, the header shows the original range), reversed sequences are reported on the reverse strand; 260 / 260 runs equal the installed program. The inputs may be the text of a GenBank or EMBL file: its features (type, location with complement / join / order, qualifiers) are read as the program reads them, follow -sbegin / -send / -sreverse, and the features that overlap a difference are listed with it (first_feature / second_feature, the 'Feature:' lines of the diffseq report); the feature files -aoutfeat / -boutfeat are written in any output feature format of the program (-offormat4 / -offormat5: gff, gff2, gff3, bed, embl, genbank, ddbj, refseq, refseqp, pir, swiss, dasgff, draw, debug and their aliases, or the 'format::' prefix of the file name; a format the program refuses for the kind of sequence - embl for protein, pir for nucleotide - gives an empty file; 1 100 runs equal the installed program); 870 runs equal the installed program (420 with the usual qualifiers, 450 with the tags the program checks with a special validator - bad values become notes, correctable ones are corrected), including its loss of the features of a GenBank entry whose SOURCE block is followed directly by FEATURES. The features of a sequence can come from a SEPARATE file: -ufo1 / -ufo2 ([format::]file), -fformat1, -fopenfile1 (and the same with 2) in reportFormat name files given in featureFiles (file name -> text), read in every input format of the EMBOSS feature readers (embl, genbank, refseq, genpept, refseqp, gff2, gff3, the gff3 of EMBOSS before 6.4, swiss, pir, dasgff; with no format every format is tried; featcopy / seqret traces of 7 800 files equal, 450 diffseq runs equal); a UFO switches the features of the sequence file off. Not read: SwissProt / GFF / GenPept files as the SEQUENCE input.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `diffseq` |
| **Method ID** | `GenomicAnalyzer.Diffseq` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 diffseq: the regions of difference between two nearly identical sequences (336 nucleotide / protein pairs, every -rformat: report and both GFF3 feature files byte for byte). The exact word matches of wordSize symbols (the engine and the minimal non-overlapping set of megamerger) are the regions of identity; the stretches between them are the differences (a difference with length 0 in one sequence is an insertion in the other); with globalDifferences the differences before the first and after the last match are reported too (a start or end that matches by itself is not counted). A nucleotide pair also gets the SNP count and its transitions / transversions. Returns the differences and the program's files: the report in the default format 'diffseq' (7 report tags; reportFormat takes any -rformat or the report qualifiers as typed, e.g. '-rformat table -raccshow'), and the GFF3 files -aoutfeat / -boutfeat. The sequence qualifiers of both inputs (-sbegin1 -send1 -sreverse1 -slower1 -supper1 -snucleotide1 ... and the same with 2, given in reportFormat as typed after the files, e.g. '-sbegin1 20 -send1 200 -sreverse2 -rformat table') are applied as the program applies them: the sequence is cut to the range (positions, lengths and the report numbers refer to the cut sequence, the header shows the original range), reversed sequences are reported on the reverse strand; 260 / 260 runs equal the installed program. The inputs may be the text of a GenBank or EMBL file: its features (type, location with complement / join / order, qualifiers) are read as the program reads them, follow -sbegin / -send / -sreverse, and the features that overlap a difference are listed with it (first_feature / second_feature, the 'Feature:' lines of the diffseq report); the feature files -aoutfeat / -boutfeat are written in any output feature format of the program (-offormat4 / -offormat5: gff, gff2, gff3, bed, embl, genbank, ddbj, refseq, refseqp, pir, swiss, dasgff, draw, debug and their aliases, or the 'format::' prefix of the file name; a format the program refuses for the kind of sequence - embl for protein, pir for nucleotide - gives an empty file; 1 100 runs equal the installed program); 870 runs equal the installed program (420 with the usual qualifiers, 450 with the tags the program checks with a special validator - bad values become notes, correctable ones are corrected), including its loss of the features of a GenBank entry whose SOURCE block is followed directly by FEATURES. The features of a sequence can come from a SEPARATE file: -ufo1 / -ufo2 ([format::]file), -fformat1, -fopenfile1 (and the same with 2) in reportFormat name files given in featureFiles (file name -> text), read in every input format of the EMBOSS feature readers (embl, genbank, refseq, genpept, refseqp, gff2, gff3, the gff3 of EMBOSS before 6.4, swiss, pir, dasgff; with no format every format is tried; featcopy / seqret traces of 7 800 files equal, 450 diffseq runs equal); a UFO switches the features of the sequence file off. Not read: SwissProt / GFF / GenPept files as the SEQUENCE input.

## Core Documentation Reference

- Algorithm: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence: nucleotide or protein letters, or the text of a FASTA, GenBank (LOCUS) or EMBL (ID) file, whose features are read. |
| `sequence2` | string | Yes | Second sequence (the same kind, or any of the file texts). |
| `wordSize` | integer | No | Word size (diffseq -wordsize, default 10, at least 2). (default: 10) |
| `globalDifferences` | boolean | No | Report the differences at the start and the end (diffseq -globaldifferences). (default: false) |
| `name1` | string | No | Name of the first sequence. (default: seqa) |
| `name2` | string | No | Name of the second sequence. (default: seqb) |
| `reportFormat` | string | No | The report format (-rformat name) or the qualifiers as typed after the files (e.g. '-rformat table -raccshow'; also the feature output qualifiers '-offormat4 embl -offormat5 gff2' / '-ofopenfile4 name', and the feature input qualifiers '-ufo1 f1.gff -fformat2 swiss -fopenfile2 f2.dat'); default format diffseq, feature files GFF3. |
| `featureFiles` | object | No | The feature files the inputs name with -ufo / -fformat / -fopenfile, by file name (for example {"f1.gff": "##gff-version 3\n..."}); a file that is not here does not exist. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `differences` | array | The regions of difference (1-based starts / ends, lengths; a length 0 is a gap). |
| `overlapStart1` | integer | Start of the overlap in the first sequence (-1: no match). |
| `overlapEnd1` | integer | End of the overlap in the first sequence. |
| `overlapStart2` | integer | Start of the overlap in the second sequence. |
| `overlapEnd2` | integer | End of the overlap in the second sequence. |
| `snps` | integer | SNP count (nucleotide pairs). |
| `transitions` | integer | Transitions among the SNPs. |
| `transversions` | integer | Transversions among the SNPs. |
| `report` | string | The report file text (-rformat, default diffseq). |
| `featureFile1` | string | The feature file of the first sequence (-aoutfeat): GFF3, or the format of -offormat4 (F220). |
| `featureFile2` | string | The feature file of the second sequence (-boutfeat): GFF3, or the format of -offormat5. |
| `featureFileName1` | string | The name of the file the first table is written to (the value of -aoutfeat without a format:: prefix, or -ofopenfile; default name.diffgff). |
| `featureFileName2` | string | The name of the file the second table is written to. |

## Errors

| Code | Message |
|------|---------|
| InvalidParams | Invalid arguments (an empty or unreadable sequence, a word size below 2, a report format the program refuses for the sequence type). |

## Examples

### Example 1: Compare two nearly identical sequences

**User Prompt:**
> Compare two nearly identical sequences

**Expected Tool Call:**
```json
{"tool": "diffseq", "arguments": {"sequence1": "ACGTACGTACGTAAGGCCTTAACCGGTT", "sequence2": "ACGTACGTACGTAAGGCATTAACCGGTT", "wordSize": 8}}
```

**Response:**
```json
{"differences": [{"start1": 18, "end1": 18, "length1": 1, "start2": 18, "end2": 18, "length2": 1}], "snps": 1, "transitions": 1, "transversions": 0}
```

## See Also

- [megamerger](megamerger.md) - merge by word matches
- [seqmatchall](seqmatchall.md) - all-against-all word matches
