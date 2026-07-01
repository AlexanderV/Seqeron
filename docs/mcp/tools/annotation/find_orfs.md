# find_orfs

Find all open reading frames (ORFs) in a DNA sequence across forward and reverse strands.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_orfs` |
| **Method ID** | `GenomeAnnotator.FindOrfs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans the (up to) six reading frames — three forward, and three on the reverse complement when
`searchBothStrands` is true — for open reading frames. An ORF runs from a start codon
(`ATG`/`GTG`/`TTG` when `requireStartCodon` is true) to the first in-frame stop codon
(`TAA`/`TAG`/`TGA`). Only ORFs whose translated length is at least `minLength` amino acids are
returned. Each ORF carries its 0-based forward-strand `start`, exclusive `end` (which includes the
stop codon), reading `frame` (1..3), `isReverseComplement` flag, the ORF nucleotide `sequence`, and
the translated `proteinSequence` (ending with `*` for the stop codon).

## Core Documentation Reference

- Source: [GenomeAnnotator.cs#L104](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs#L104)
- Algorithm: [ORF_Detection.md](../../../algorithms/Annotation/ORF_Detection.md)

## Input Schema

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `dnaSequence` | string | Yes | — | DNA sequence to search |
| `minLength` | integer | No | 100 | Minimum ORF length in amino acids |
| `searchBothStrands` | boolean | No | true | Whether to also search the reverse complement |
| `requireStartCodon` | boolean | No | true | Whether to require a start codon (ATG/GTG/TTG) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `orfs` | array | `{ start, end, frame, isReverseComplement, sequence, proteinSequence }` per ORF |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Simple ATG..TAA ORF

`ATGAAAAAAAAATAA` (start `ATG`, three `AAA` codons, stop `TAA`), `minLength=1`, forward strand only →
one ORF spanning positions 0..15, protein `MKKK*`.

**Response:**
```json
{ "orfs": [ { "start": 0, "end": 15, "frame": 1, "isReverseComplement": false, "sequence": "ATGAAAAAAAAATAA", "proteinSequence": "MKKK*" } ] }
```

### Example 2: Below minimum length

The same sequence with `minLength=10` yields no ORF (the ORF is only 4 aa including the stop).

**Response:**
```json
{ "orfs": [] }
```

## Performance

- **Time Complexity:** O(n) per frame
- **Space Complexity:** O(k) for the returned ORFs

## See Also

- [longest_orfs_per_frame](longest_orfs_per_frame.md) - Longest ORF per reading frame
- [predict_genes](predict_genes.md) - ORF-based gene prediction
- [find_ribosome_binding_sites](find_ribosome_binding_sites.md) - Shine–Dalgarno motifs upstream of ORFs
