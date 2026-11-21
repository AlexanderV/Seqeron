# write_coordinate_alignment

Writes alignments given by their block coordinates in a Biopython 1.88 Bio.Align coordinate format, byte for byte as Bio.Align.write does: bed (3 to 12 columns, bedN), psl (header, mask lower / upper for repeat matches, wildcard N), chain, maf (any number of sequences, score / pass / empty-region annotations, quality and i lines of the records, file header metadata), mauve (xmfa: header from the metadata and identifiers, one block per alignment). An alignment is a JSON object {"sequences": [{"id": "chr1", "length": 1000, "parts": {"0": "ACGT..."} or null}, {query}], "coordinates": [[target positions], [query positions]], "attrs": {"score": 42, "matches": 10, ...}} (target first, query second; a position that decreases along the alignment is a sequence on the reverse strand; parts = the known residues by start position, null = unknown contents; PSL counts matches / mismatches from the residues when they are known, otherwise every aligned base is a match). VERIFICATION: written text equal to Biopython 1.88 on 620 / 620 committed cases (scripts/alignment/generate_biopython_coord_oracle.py).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `write_coordinate_alignment` |
| **Method ID** | `AlignmentCoordinateFormats.Write` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Writes alignments given by their block coordinates in a Biopython 1.88 Bio.Align coordinate format, byte for byte as Bio.Align.write does: bed (3 to 12 columns, bedN), psl (header, mask lower / upper for repeat matches, wildcard N), chain, maf (any number of sequences, score / pass / empty-region annotations, quality and i lines of the records, file header metadata), mauve (xmfa: header from the metadata and identifiers, one block per alignment). An alignment is a JSON object {"sequences": [{"id": "chr1", "length": 1000, "parts": {"0": "ACGT..."} or null}, {query}], "coordinates": [[target positions], [query positions]], "attrs": {"score": 42, "matches": 10, ...}} (target first, query second; a position that decreases along the alignment is a sequence on the reverse strand; parts = the known residues by start position, null = unknown contents; PSL counts matches / mismatches from the residues when they are known, otherwise every aligned base is a match). VERIFICATION: written text equal to Biopython 1.88 on 620 / 620 committed cases (scripts/alignment/generate_biopython_coord_oracle.py).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/AlignmentCoordinateFormats.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/AlignmentCoordinateFormats.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md) (section Biopython coordinate formats)
- Oracle: `scripts/alignment/generate_biopython_coord_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/biopython_coord_oracle.tsv`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `alignmentsJson` | array<string> | Yes | The alignments, each a JSON object as described (a JSON array of objects). |
| `format` | string | No | Format: bed, psl, chain, maf or mauve. (default: psl) |
| `bedN` | integer | No | BED: number of columns, 3 to 12. (default: 12) |
| `header` | boolean | No | PSL: write the five header lines. (default: true) |
| `mask` | string | No | PSL: null, 'lower' or 'upper': repeat-masked residues of the target are counted as repMatches. |
| `wildcard` | string | No | PSL: the wildcard residue counted in nCount (default N; empty = none). |
| `identifiers` | array<string> | No | Mauve: the identifiers of the sequences (the file names, or the entry numbers of a combined file as strings: '0', '1', ...); sequences are numbered by their index. |
| `metadataJson` | string | No | MAF / Mauve: the metadata of the file header as a JSON object (MAF: {"MAF Version": "1", "Scoring": ..., "Program": ..., "Comments": [...], "name", "description", "mafDot", "visibility", "speciesOrder": [...]}, default {"MAF Version": "1"}; Mauve: {"FormatVersion": "Mauve1", "File": "genome.fa" (one combined file), "BackboneFile": ...}). |

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
