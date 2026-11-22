# seqret

EMBOSS 6.6.0 seqret: the sequences of a file are read by the sequence readers of the EMBOSS library (ajseqread.c) and written in FASTA (with the GFF3 or another feature format of the features the file carries). Input formats read byte for byte as the installed program reads them (200 inputs of each, with the library's warnings and errors): fasta (ncbi, pearson, gifasta, dbid), gcg (gcg8), embl (em), swiss (sw, swissprot, uniprot; the DE description of RecName / AltName / SubName / Includes / Contains / Flags is assembled as the library does), nbrf (pir), genbank (gb, ddbj, refseq), refseqp, genpept, codata, strider, ig, igstrict, iguspto, staden, textonly (plain, asis, unknown), fastq, fastq-sanger, fastq-illumina, fastq-solexa, sam, gff2, gff3 (gff; the sequence after ##FASTA / the ##DNA lines), raw. With no format every format the program tries is tried in turn (a text no format reads is read as plain text). The features of embl, genbank, refseq, refseqp, swiss, nbrf, gff2 and gff3 files are returned with featureFormat when feature is true. Not read yet (alignment and binary formats: clustal, msf, phylip, nexus, mega, selex, stockholm, mase, fitch, treecon, hennig86, jackknifer, gde, acedb, ace, pdb, abi, bam, nibble, scf, das, ensembl): the tool refuses an input of such a format. A5-8c3g, B09 F248.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `seqret` |
| **Method ID** | `EmbossSeqretPrograms.Seqret` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 seqret: the sequences of a file are read by the sequence readers of the EMBOSS library (ajseqread.c) and written in FASTA (with the GFF3 or another feature format of the features the file carries). Input formats read byte for byte as the installed program reads them (200 inputs of each, with the library's warnings and errors): fasta (ncbi, pearson, gifasta, dbid), gcg (gcg8), embl (em), swiss (sw, swissprot, uniprot; the DE description of RecName / AltName / SubName / Includes / Contains / Flags is assembled as the library does), nbrf (pir), genbank (gb, ddbj, refseq), refseqp, genpept, codata, strider, ig, igstrict, iguspto, staden, textonly (plain, asis, unknown), fastq, fastq-sanger, fastq-illumina, fastq-solexa, sam, gff2, gff3 (gff; the sequence after ##FASTA / the ##DNA lines), raw. With no format every format the program tries is tried in turn (a text no format reads is read as plain text). The features of embl, genbank, refseq, refseqp, swiss, nbrf, gff2 and gff3 files are returned with featureFormat when feature is true. Not read yet (alignment and binary formats: clustal, msf, phylip, nexus, mega, selex, stockholm, mase, fitch, treecon, hennig86, jackknifer, gde, acedb, ace, pdb, abi, bam, nibble, scf, das, ensembl): the tool refuses an input of such a format. A5-8c3g, B09 F248.

## Core Documentation Reference

- Algorithm: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file (any input format above), or a sequence as letters. A value 'asis::ACGT' is the sequence itself. |
| `sequenceFormat` | string | No | The input format (-sformat): the names above, any case, or an EDAM term; default: the program finds the format. |
| `feature` | boolean | No | Read and write the features of the sequences (-feature). (default: false) |
| `firstOnly` | boolean | No | Read one sequence and stop (-firstonly). (default: false) |
| `featureFormat` | string | No | The feature output format (-offormat2): gff3 (default), gff, gff2, bed, embl, em, genbank, gb, ddbj, refseq, refseqp, pir, nbrf, swiss, sw, swissprot, uniprot, dasgff, draw, debug. |
| `qualifiers` | string | No | Sequence input qualifiers as typed after the sequence on the command line: -slower1 -supper1 -snucleotide1 -sprotein1 (also -sformat1 F); the range qualifiers are not applied by this tool yet. |
| `files` | object | No | More files by name (file name -> text) that sequence may name with a 'format::' prefix. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sequenceText` | string | The text of the sequence output file (FASTA). |
| `featureText` | string | The text of the feature output file (empty when no feature was written). |
| `sequences` | integer | The number of sequences written. |
| `messages` | array | What the library printed on its error stream (Warning: ..., Error: ...), in order. |

## Errors

| Code | Message |
|------|---------|
| InvalidParams | Invalid arguments (an empty sequence, a text no sequence reader reads - 'Unable to read sequence' -, an input of an alignment or binary format this port does not read, an input on which the reference program dies by a fault). |

## Examples

### Example 1: A Swiss-Prot entry with its description and features

**User Prompt:**
> Read this Swiss-Prot entry with its features

**Expected Tool Call:**
```json
{"tool": "seqret", "arguments": {"sequence": "ID   ABC_HUMAN   Reviewed;   10 AA.\nAC   P12345;\nDE   RecName: Full=Kinase alpha; Short=KA;\nDE   AltName: Full=Other name;\nFT   DOMAIN       1      5       Kinase domain.\nSQ   SEQUENCE   10 AA;  1000 MW;  0 CRC64;\n     MKVLAAGIVG\n//\n", "sequenceFormat": "swiss", "feature": true}}
```

**Response:**
```json
{"sequenceText": ">ABC_HUMAN P12345 Kinase alpha (KA) (Other name)\nMKVLAAGIVG\n", "featureText": "##gff-version 3\n##sequence-region ABC_HUMAN 1 10\n#!Date 2026-11-22\n#!Type Protein\n#!Source-version EMBOSS 6.6.0.0\nABC_HUMAN\tSWISSPROT\tpolypeptide_domain\t1\t5\t.\t.\t.\tID=ABC_HUMAN.1;note=Kinase domain\n", "sequences": 1, "messages": []}
```

## See Also

- [featcopy](featcopy.md) - read and write feature tables
- [feattext](feattext.md) - original text of a feature table
- [diffseq](diffseq.md) - differences between two sequences
