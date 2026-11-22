# seqret

EMBOSS 6.6.0 seqret: the sequences of a file are read by the sequence readers of the EMBOSS library (ajseqread.c) and written in FASTA or, with outputFormat, in another sequence format (with the GFF3 or another feature format of the features the file carries). Input formats read byte for byte as the installed program reads them (200 inputs of each of 66 names, with the library's warnings and errors): fasta (ncbi, pearson, gifasta, dbid), gcg (gcg8), embl (em), swiss (sw, swissprot, uniprot; the DE description of RecName / AltName / SubName / Includes / Contains / Flags is assembled as the library does), nbrf (pir), genbank (gb, ddbj, refseq), refseqp, genpept, codata, strider, ig, igstrict, iguspto, staden, textonly (plain, asis, unknown), fastq, fastq-sanger, fastq-illumina, fastq-solexa, sam, gff2, gff3 (gff; the sequence after ##FASTA / the ##DNA lines), gde, raw, and the alignment formats clustal (aln), msf, fitch, mase, biomart, treecon, acedb, phylip (interleaved; falls back to phylipnon as the library does), phylipnon, nexus (paup), mega, selex, stockholm (pfam), jackknifer, hennig86 (the library loops over memory on a sequence line: refused), the entries of experiment, pdb (the ATOM residues), pdbseq (SEQRES), pdbnuc, pdbnucseq (the members of the alignment are returned one after the other). ace / consed (the library reads the AS line and then dies by a segmentation fault) are refused, as are the few inputs on which the library reads freed or foreign memory (a repeated taxon label of a nexus file, ...). With no format every format the program tries is tried in turn (a text no format reads is read as plain text). The features of embl, genbank, refseq, refseqp, swiss, nbrf, gff2 and gff3 files are returned with featureFormat when feature is true. Output: 19 format names (fasta unknown pearson, gcg gcg8, nbrf pir, codata, strider, staden, text plain raw, ig, gde, fitch, acedb, mase) written byte for byte as the installed program writes them (200 runs of each: the format as -osformat2 or the 'format::' prefix of the output name, outputSingle with outputDirectory / outputExtension, unknown formats, the format-for-sequence-type errors). Not read yet (the binary abi, bam, nibble, scf, das, ensembl and the CO contig lines of embl) and not written yet (embl, genbank, swiss, gff, fastq, sam and the alignment-set formats): the tool refuses such an input or format. A5-8c3g / A5-8c3h, B09 F248 / F251 - F252.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `seqret` |
| **Method ID** | `EmbossSeqretPrograms.Seqret` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 seqret: the sequences of a file are read by the sequence readers of the EMBOSS library (ajseqread.c) and written in FASTA or, with outputFormat, in another sequence format (with the GFF3 or another feature format of the features the file carries). Input formats read byte for byte as the installed program reads them (200 inputs of each of 66 names, with the library's warnings and errors): fasta (ncbi, pearson, gifasta, dbid), gcg (gcg8), embl (em), swiss (sw, swissprot, uniprot; the DE description of RecName / AltName / SubName / Includes / Contains / Flags is assembled as the library does), nbrf (pir), genbank (gb, ddbj, refseq), refseqp, genpept, codata, strider, ig, igstrict, iguspto, staden, textonly (plain, asis, unknown), fastq, fastq-sanger, fastq-illumina, fastq-solexa, sam, gff2, gff3 (gff; the sequence after ##FASTA / the ##DNA lines), gde, raw, and the alignment formats clustal (aln), msf, fitch, mase, biomart, treecon, acedb, phylip (interleaved; falls back to phylipnon as the library does), phylipnon, nexus (paup), mega, selex, stockholm (pfam), jackknifer, hennig86 (the library loops over memory on a sequence line: refused), the entries of experiment, pdb (the ATOM residues), pdbseq (SEQRES), pdbnuc, pdbnucseq (the members of the alignment are returned one after the other). ace / consed (the library reads the AS line and then dies by a segmentation fault) are refused, as are the few inputs on which the library reads freed or foreign memory (a repeated taxon label of a nexus file, ...). With no format every format the program tries is tried in turn (a text no format reads is read as plain text). The features of embl, genbank, refseq, refseqp, swiss, nbrf, gff2 and gff3 files are returned with featureFormat when feature is true. Output: 19 format names (fasta unknown pearson, gcg gcg8, nbrf pir, codata, strider, staden, text plain raw, ig, gde, fitch, acedb, mase) written byte for byte as the installed program writes them (200 runs of each: the format as -osformat2 or the 'format::' prefix of the output name, outputSingle with outputDirectory / outputExtension, unknown formats, the format-for-sequence-type errors). Not read yet (the binary abi, bam, nibble, scf, das, ensembl and the CO contig lines of embl) and not written yet (embl, genbank, swiss, gff, fastq, sam and the alignment-set formats): the tool refuses such an input or format. A5-8c3g / A5-8c3h, B09 F248 / F251 - F252.

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
| `qualifiers` | string | No | Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1 (the range and the strand: applied to the sequences and to their features), -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sformat1 F, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes. |
| `files` | object | No | More files by name (file name -> text) that sequence may name with a 'format::' prefix. |
| `outputFormat` | string | No | The sequence output format (-osformat2): fasta (default; unknown, pearson), gcg (gcg8), nbrf (pir), codata, strider, staden, text (plain, raw), ig, gde, fitch, acedb, mase. iguspto is accepted by the library but its writer dies on every sequence that has no full description (all but swiss and iguspto input), so such a run is refused. Other output formats of the library (embl, genbank, swiss, msf, clustal, phylip, nexus, ...) are not written yet. |
| `outputSingle` | boolean | No | Write every sequence to a file of its own, name.extension in lower case (-ossingle2); the files are returned in singleFiles. (default: false) |
| `outputDirectory` | string | No | The directory of the files of outputSingle (-osdirectory2). |
| `outputExtension` | string | No | The extension of the files of outputSingle (-osextension2); default: the output format as typed. |
| `sequenceValue` | string | No | The USA of the sequence input (-sequence) as the program reads it, instead of the file the text is put in: [format::]file[begin:end:r] (file = the name of an entry of files; the text of sequence is the file 'sequence.dat'), file:id / file:field:query / file{a\|b} (id, acc, sv, gi, des, key, org; wildcards * and ?), @list / list::list (a list file of files; # comments; every entry has its own range and format), asis::SEQUENCE; the range qualifiers -sbegin1 -send1 -sreverse1 of qualifiers are applied to the sequences and to their features. |
| `definitions` | string | No | The text of an EMBOSS configuration file (.embossrc) that defines databases for sequenceValue: DBNAME name [ type: N format: fasta method: direct directory: "dir" filename: "*.fa" ], SET name value; the data files are in files. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sequenceText` | string | The text of the sequence output file (FASTA). |
| `featureText` | string | The text of the feature output file (empty when no feature was written). |
| `sequences` | integer | The number of sequences written. |
| `messages` | array | What the library printed on its error stream (Warning: ..., Error: ...), in order. |
| `singleFiles` | object | With outputSingle: the file of every sequence by path (directory/name.extension); null otherwise. |

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
