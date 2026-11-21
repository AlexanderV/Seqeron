# write_coordinate_alignment

Writes pairwise alignments given by their block coordinates in a Biopython 1.88 Bio.Align coordinate format, byte for byte as Bio.Align.write does: bed (3 to 12 columns, bedN), psl (header, mask lower / upper for repeat matches, wildcard N), chain. An alignment is a JSON object {"sequences": [{"id": "chr1", "length": 1000, "parts": {"0": "ACGT..."} or null}, {query}], "coordinates": [[target positions], [query positions]], "attrs": {"score": 42, "matches": 10, ...}} (target first, query second; a position that decreases along the alignment is a sequence on the reverse strand; parts = the known residues by start position, null = unknown contents; PSL counts matches / mismatches from the residues when they are known, otherwise every aligned base is a match). VERIFICATION: written text equal to Biopython 1.88 on 390 / 390 committed cases (scripts/alignment/generate_biopython_coord_oracle.py).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `write_coordinate_alignment` |
| **Method ID** | `AlignmentCoordinateFormats.Write` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Writes pairwise alignments given by their block coordinates in a Biopython 1.88 Bio.Align coordinate format, byte for byte as Bio.Align.write does: bed (3 to 12 columns, bedN), psl (header, mask lower / upper for repeat matches, wildcard N), chain. An alignment is a JSON object {"sequences": [{"id": "chr1", "length": 1000, "parts": {"0": "ACGT..."} or null}, {query}], "coordinates": [[target positions], [query positions]], "attrs": {"score": 42, "matches": 10, ...}} (target first, query second; a position that decreases along the alignment is a sequence on the reverse strand; parts = the known residues by start position, null = unknown contents; PSL counts matches / mismatches from the residues when they are known, otherwise every aligned base is a match). VERIFICATION: written text equal to Biopython 1.88 on 390 / 390 committed cases (scripts/alignment/generate_biopython_coord_oracle.py).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/AlignmentCoordinateFormats.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/AlignmentCoordinateFormats.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md) (section Biopython coordinate formats)
- Oracle: `scripts/alignment/generate_biopython_coord_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/biopython_coord_oracle.tsv`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `alignmentsJson` | array<string> | Yes | The alignments, each a JSON object as described (a JSON array of objects). |
| `format` | string | No | Format: bed, psl or chain. (default: psl) |
| `bedN` | integer | No | BED: number of columns, 3 to 12. (default: 12) |
| `header` | boolean | No | PSL: write the five header lines. (default: true) |
| `mask` | string | No | PSL: null, 'lower' or 'upper': repeat-masked residues of the target are counted as repMatches. |
| `wildcard` | string | No | PSL: the wildcard residue counted in nCount (default N; empty = none). |

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
