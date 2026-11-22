# coderet

EMBOSS 6.6.0 'coderet': the coding sequences, messenger RNAs, translations and non-coding regions of the sequences of a file with feature tables.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `coderet` |
| **Method ID** | `EmbossCoderetProgram.Coderet` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 coderet: from the nucleotide sequences of a file with feature tables (embl, genbank, refseq, gff2, gff3, ...) the coding sequences (the CDS features, joins as one sequence; files NAME.cds, sequences ACCESSION_cds_N), the messenger RNAs (mRNA features; NAME.mrna, _mrna_N), the translations (the /translation tag of the CDS features; NAME.prot, _pro_N) and the non-coding regions (what the CDS, mRNA and exon features leave; NAME.noncoding, _noncoding_START) are written, NAME being the first sequence in lower case. The tool returns the log of the program (the counts of every sequence), the files by name and the messages; cds, mrna, translation and nonCoding select the outputs and outputFormat is the sequence output format (fasta unknown pearson, gcg gcg8, nbrf pir, codata, strider, staden, text plain raw, ig iguspto, gde, fitch, acedb, mase, ncbi, gifasta, nib (nibble), dasdna, das, asn1, embl (em, emblnew, emold, emblold), genbank (gb, ddbj), refseq, genpept, refseqp, swiss (sw, swissprot, uniprot, swissnew, swnew, swold, swissold, swissprotold, swissprotnew), gff2, gff3 (gff), fastq (fastq-sanger, fastq-illumina, fastq-solexa: the quality values of the sequences are written), sam, scf, and the sequence-set formats msf, clustal (aln), selex, phylip, phylipnon (phylip3), hennig86, mega, meganon, nexus (paup), nexusnon (paupnon), jackknifer, jackknifernon, treecon (written when the output is closed), experiment, debug and bam (BGZF blocks deflated as zlib 1.3 does)). Equals the installed program on 1 000 runs (EMBL entries with joins, complements, orders, partial ends and translations, the reader inputs with features, damaged files, one to four outputs in 29 format names). A5-8c3h, B09 F259.

## Core Documentation Reference

- Source: [EmbossCoderetProgram.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/EmbossCoderetProgram.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_coderet_oracle.py](../../../../scripts/alignment/generate_coderet_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The text of the sequence file with features (nucleotide sequences), see the formats above. |
| `sequenceFormat` | string | No | The input format (-sformat); default: the program finds the format. |
| `cds` | boolean | No | Write the coding sequences (-cdsoutseq). (default: true) |
| `mrna` | boolean | No | Write the messenger RNAs (-mrnaoutseq). (default: true) |
| `translation` | boolean | No | Write the translations (-translationoutseq). (default: true) |
| `nonCoding` | boolean | No | Write the non-coding regions (-restoutseq). (default: true) |
| `outputFormat` | string | No | The sequence output format of the four outputs (-osformat3 .. -osformat6); default fasta. |
| `qualifiers` | string | No | Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1, -sformat1 F, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes. |
| `files` | object | No | More files by name (file name -> text) that sequenceValue may name. |
| `sequenceValue` | string | No | The USA of the sequence input (-seqall): [format::]file[begin:end:r], file:id, file:field:value, file{a\|b}, @list, list::list, db:id; the text of sequence is the file 'sequence.dat'. |
| `definitions` | string | No | The text of an EMBOSS configuration file (.embossrc) that defines databases for sequenceValue. |
| `outputSingle` | boolean | No | Write every sequence of the four outputs to a file of its own, name.extension in lower case (-ossingle3 .. -ossingle6); the files are in files. (default: false) |
| `outputDirectory` | string | No | The directory of the output files (-osdirectory3 .. -osdirectory6). |
| `outputExtension` | string | No | The extension of the output files (-osextension3 .. -osextension6); default cds, mrna, prot and noncoding. |
| `outputName` | string | No | The base name of the output files (-osname3 .. -osname6); default: the name of the first sequence in lower case. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `log` | string | The log of the program: the counts of CDS, mRNA, non-coding features and translations of every sequence. |
| `files` | object | The sequence files by name. |
| `messages` | string[] | What the library printed (warnings and errors). |

## Errors

| Code | Message |
|------|---------|
| 1001 | The sequence cannot be empty |
| 1001 | Unable to read sequence ... (no sequence in the file) and Sequence is not nucleic (a protein file) |
| 1001 | An input on which the reference program never stops or dies by a fault, and an output format this port does not write, are refused |

## Examples

### Example 1: The coding sequence and translation of an EMBL entry

**User Prompt:**
> Get the CDS and its translation from this EMBL entry.

**Expected Tool Call:**
```json
{"tool": "coderet", "arguments": {"sequence": "ID   X1; SV 1; linear; genomic DNA; STD; UNC; 60 BP.\nXX\nAC   X1;\nXX\nDE   test\nXX\nFH   Key             Location/Qualifiers\nFH\nFT   CDS             join(3..12,20..29)\nFT                   /translation=\"MKV\"\nSQ   Sequence 60 BP; 15 A; 15 C; 15 G; 15 T; 0 other;\n     acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt        60\n//\n", "mrna": false, "nonCoding": false}}
```

**Response:**
```json
{"log": "   CDS Trans Total Sequence\n ===== ===== ===== ========\n     1     1     1 X1\n", "files": {"x1.cds": ">x1_cds_1\ngtacgtacgttacgtacgta\n", "x1.prot": ">x1_pro_1\nMKV\n"}, "messages": []}
```

The sequences are named ACCESSION_cds_N and ACCESSION_pro_N in lower case; the CDS of a join is its parts one after the other.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
