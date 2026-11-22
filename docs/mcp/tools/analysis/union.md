# union

EMBOSS 6.6.0 union: the sequences of a file (any input format of seqret) are joined into ONE sequence that carries the name, description, accession and type of the first; with findOverlap the bases by which a sequence repeats the end of the one before (the longest end of the previous sequence that is the start of this one) are joined once, and with overlapFile the overlap of every sequence is returned (one number per sequence, the first is 0). The result is written in FASTA or in the output format of outputFormat (fasta unknown pearson, gcg gcg8, nbrf pir, codata, strider, staden, text plain raw, ig iguspto, gde, fitch, acedb, mase, ncbi, gifasta, nib (nibble), dasdna, das, asn1, embl (em, emblnew, emold, emblold), genbank (gb, ddbj), refseq, genpept, refseqp, swiss (sw, swissprot, uniprot, swissnew, swnew, swold, swissold, swissprotold, swissprotnew), gff2, gff3 (gff), fastq (fastq-sanger, fastq-illumina, fastq-solexa: the quality values of the sequences are written), sam, scf, and the sequence-set formats msf, clustal (aln), selex, phylip, phylipnon (phylip3), hennig86, mega, meganon, nexus (paup), nexusnon (paupnon), jackknifer, jackknifernon, treecon (written when the output is closed), experiment, debug and bam (BGZF blocks deflated as zlib 1.3 does)). With feature the features of the sequences are joined into the features of the result (the positions moved by the length of the sequences before; with source one feature of the type source, tag origid = the name of the sequence, for every sequence joined) and returned as featureText in the feature format of featureFormat (gff3 default, gff2, embl, genbank, refseq, bed, dasgff, debug; an unknown name is written as gff3 with the library's errors). Equals the installed program on 200 runs (fragments that overlap by random lengths, the reader inputs, damaged files, every output format, the format-for-sequence-type errors). A5-8c3h, B09 F253. The joined sequence can be written to a file of its own (outputSingle with outputDirectory / outputExtension: -ossingle2 -osdirectory2 -osextension2; returned in singleFiles), as seqret does; equals the installed program on 1 200 runs (every output format including experiment, debug and bam, the output qualifiers, directories that do not exist). A5-8c3h, B09 F277.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `union` |
| **Method ID** | `EmbossUnionProgram.Union` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 union: the sequences of a file (any input format of seqret) are joined into ONE sequence that carries the name, description, accession and type of the first; with findOverlap the bases by which a sequence repeats the end of the one before (the longest end of the previous sequence that is the start of this one) are joined once, and with overlapFile the overlap of every sequence is returned (one number per sequence, the first is 0). The result is written in FASTA or in the output format of outputFormat (fasta unknown pearson, gcg gcg8, nbrf pir, codata, strider, staden, text plain raw, ig iguspto, gde, fitch, acedb, mase, ncbi, gifasta, nib (nibble), dasdna, das, asn1, embl (em, emblnew, emold, emblold), genbank (gb, ddbj), refseq, genpept, refseqp, swiss (sw, swissprot, uniprot, swissnew, swnew, swold, swissold, swissprotold, swissprotnew), gff2, gff3 (gff), fastq (fastq-sanger, fastq-illumina, fastq-solexa: the quality values of the sequences are written), sam, scf, and the sequence-set formats msf, clustal (aln), selex, phylip, phylipnon (phylip3), hennig86, mega, meganon, nexus (paup), nexusnon (paupnon), jackknifer, jackknifernon, treecon (written when the output is closed), experiment, debug and bam (BGZF blocks deflated as zlib 1.3 does)). With feature the features of the sequences are joined into the features of the result (the positions moved by the length of the sequences before; with source one feature of the type source, tag origid = the name of the sequence, for every sequence joined) and returned as featureText in the feature format of featureFormat (gff3 default, gff2, embl, genbank, refseq, bed, dasgff, debug; an unknown name is written as gff3 with the library's errors). Equals the installed program on 200 runs (fragments that overlap by random lengths, the reader inputs, damaged files, every output format, the format-for-sequence-type errors). A5-8c3h, B09 F253. The joined sequence can be written to a file of its own (outputSingle with outputDirectory / outputExtension: -ossingle2 -osdirectory2 -osextension2; returned in singleFiles), as seqret does; equals the installed program on 1 200 runs (every output format including experiment, debug and bam, the output qualifiers, directories that do not exist). A5-8c3h, B09 F277.

## Core Documentation Reference

- Algorithm: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file (any input format of seqret), or a sequence as letters. |
| `sequenceFormat` | string | No | The input format (-sformat); default: the program finds the format. |
| `findOverlap` | boolean | No | Look for the overlaps of consecutive sequences when joining (-findoverlap). (default: false) |
| `overlapFile` | boolean | No | Return the overlap of every sequence (-overlapfile). (default: false) |
| `outputFormat` | string | No | The sequence output format (-osformat); default fasta. (default: fasta) |
| `qualifiers` | string | No | Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sformat1 F, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes. |
| `files` | object | No | More files by name (file name -> text) that sequence may name with a 'format::' prefix. |
| `sequenceValue` | string | No | The USA of the sequence input (-sequence) as the program reads it, instead of the file the text is put in: [format::]file[begin:end:r] (file = the name of an entry of files; the text of sequence is the file 'sequence.dat'), file:id / file:field:query / file{a\|b} (id, acc, sv, gi, des, key, org; wildcards * and ?), @list / list::list (a list file of files; # comments; every entry has its own range and format), asis::SEQUENCE; the range qualifiers -sbegin1 -send1 -sreverse1 of qualifiers are applied to the sequences and to their features. |
| `feature` | boolean | No | Join the features of the sequences (-feature); the features are returned as featureText. (default: false) |
| `source` | boolean | No | Create a source feature for every sequence joined (-source; with feature). (default: false) |
| `featureFormat` | string | No | The format of the feature output (-offormat2); default gff3 (gff). |
| `outputSingle` | boolean | No | Write the joined sequence to a file of its own, name.extension in lower case (-ossingle2); the file is returned in singleFiles. (default: false) |
| `outputDirectory` | string | No | The directory of the file of outputSingle (-osdirectory2). |
| `outputExtension` | string | No | The extension of the file of outputSingle (-osextension2); default: the output format as typed. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sequenceText` | string | The text of the sequence output file. |
| `overlapText` | string | The overlap of every sequence, one number per line (empty without overlapFile). |
| `messages` | array | What the library printed on its error stream (Warning: ..., Error: ...), in order. |

## Errors

| Code | Message |
|------|---------|
| InvalidParams | Invalid arguments (an empty sequence, a text no sequence reader reads, an input of an alignment or binary format this port does not read). |

## Examples

### Example 1: Join overlapping fragments

**User Prompt:**
> Join these fragments and count the overlaps

**Expected Tool Call:**
```json
{"tool": "union", "arguments": {"sequence": ">f1\nACGTACGTAC\n>f2\nCGTACGTACGGG\n>f3\nGGGTTTAAA\n", "findOverlap": true, "overlapFile": true}}
```

**Response:**
```json
{"sequenceText": ">f1\nACGTACGTACGGGTTTAAA\n", "overlapText": "0\n9\n3\n", "messages": []}
```

## See Also

- [seqret](seqret.md) - read and write sequences
- [maskfeat](maskfeat.md) - mask the regions of features
- [descseq](descseq.md) - alter the name or description of a sequence
