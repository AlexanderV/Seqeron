# parse_coordinate_alignment

Reads the alignments of a BED, PSL / PSLX (nucleotides, and a protein aligned to translated DNA with its CDS feature), chain, MAF, Mauve xmfa, SAM (header metadata and targets, MD tag, hard / soft clips, qualities, tags), Exonerate vulgar / cigar, GCG MSF, EMBOSS pair / srspair / simple (needle, water, stretcher), HHsearch hhr or BLAST / FASTA tabular (-outfmt 7, -m 8CB / 8CC) file as Bio.Align.parse does and returns them as JSON objects {"sequences", "coordinates", "attrs"} (see write_coordinate_alignment): sequences with their lengths and the residues a PSLX file defines, the block coordinates (a decreasing coordinate = reverse strand), score / matches / misMatches / repMatches / nCount / thickStart / thickEnd / itemRgb / annotations.id as the format carries them. Text Biopython refuses (inconsistent columns, block counts, starts and ends) is refused with the same condition; the file-level values (the iterator's metadata, the SAM targets, the parser warnings) are returned in metadataJson / targetsJson / warnings. VERIFICATION: equal to Biopython 1.88 on 2 331 / 2 331 committed texts (incl. EMBOSS 6.6.0, HH-suite 3.3.0, BLAST+ and fasta36 outputs).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `parse_coordinate_alignment` |
| **Method ID** | `AlignmentCoordinateFormats.Parse` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Reads the alignments of a BED, PSL / PSLX (nucleotides, and a protein aligned to translated DNA with its CDS feature), chain, MAF, Mauve xmfa, SAM (header metadata and targets, MD tag, hard / soft clips, qualities, tags), Exonerate vulgar / cigar, GCG MSF, EMBOSS pair / srspair / simple (needle, water, stretcher), HHsearch hhr or BLAST / FASTA tabular (-outfmt 7, -m 8CB / 8CC) file as Bio.Align.parse does and returns them as JSON objects {"sequences", "coordinates", "attrs"} (see write_coordinate_alignment): sequences with their lengths and the residues a PSLX file defines, the block coordinates (a decreasing coordinate = reverse strand), score / matches / misMatches / repMatches / nCount / thickStart / thickEnd / itemRgb / annotations.id as the format carries them. Text Biopython refuses (inconsistent columns, block counts, starts and ends) is refused with the same condition; the file-level values (the iterator's metadata, the SAM targets, the parser warnings) are returned in metadataJson / targetsJson / warnings. VERIFICATION: equal to Biopython 1.88 on 2 331 / 2 331 committed texts (incl. EMBOSS 6.6.0, HH-suite 3.3.0, BLAST+ and fasta36 outputs).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/AlignmentCoordinateFormats.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/AlignmentCoordinateFormats.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md) (section Biopython coordinate formats)
- Oracle: `scripts/alignment/generate_biopython_coord_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/biopython_coord_oracle.tsv`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `text` | string | Yes | The text of the file. |
| `format` | string | No | Format: bed, psl, chain, maf, mauve, sam, exonerate, msf, emboss, hhr or tabular. (default: psl) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The text of the file (write) or the JSON of the parsed alignments (parse). |
| `format` | string | The format name. |
| `count` | integer | The number of alignments. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Unknown coordinate format |

## See Also

- [write_alignment](write_alignment.md)
