# assess_completeness

Assess assembly completeness with nucleotide marker genes (BUSCO-like).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `assess_completeness` |
| **Method ID** | `GenomeAssemblyAnalyzer.AssessCompleteness` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Each marker and its reverse complement are seeded with exact 28-bp ACGT words (BLAST+ megablast
`word_size`; a marker shorter than 28 bp is seeded by its whole length) and every seed cluster is
extended by Smith-Waterman-Gotoh local alignment with the BLAST+ blastn scoring (+2/−3, gap 5 + 2·k).
A hit has identity = identical / alignment columns and coverage = aligned marker span / marker length.

Classification follows compleasm (`record_1st_gene_label`, aligned-fraction mode):
a hit is valid when coverage ≥ `minLengthFraction` **or** identity ≥ `identityThreshold`; no valid hit →
Missing; complete hits have coverage ≥ `coverageThreshold`; none → Fragmented; one → Single; several →
Duplicated when they lie on different sequences or form more than one cluster of overlapping regions on
one sequence (abutting regions merge, as compleasm `get_region_clusters` on half-open coordinates).

Percentages are BUSCO's, all over n: C% = (S + D)/n, D% = D/n, and `buscoNotation` is BUSCO's
`C:x%[S:x%,D:x%],F:x%,M:x%,n:N` (Python `round(·, 1)`).

This is not BUSCO itself: BUSCO/compleasm search lineage protein HMM profiles through gene predictors
or miniprot; here the caller's nucleotide markers are aligned directly, so introns/frameshifts are not
modelled and the 28-bp exact seed is the only significance filter (a marker diverged so that no exact
28-bp word survives is Missing, as with megablast).

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L727](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L727)
- Evidence: [ASSEMBLY-COMPLETE-001](../../../Evidence/ASSEMBLY-COMPLETE-001-Evidence.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `assembly` | array | Yes | Assembled sequences `{ id, sequence }` (case-insensitive). |
| `markerGenes` | array | Yes | Marker genes `{ id, sequence }` (nucleotide, non-empty). |
| `identityThreshold` | number | No | compleasm `min_identity` in `[0, 1]` (default 0.9). |
| `coverageThreshold` | number | No | compleasm `min_complete` in `[0, 1]` (default 0.9). |
| `minLengthFraction` | number | No | compleasm `min_length_percent` in `[0, 1]` (default 0.6). |

## Output Schema

Fields: `totalGenes` (n), `complete` (C = S + D), `completeSingleCopy` (S), `completeDuplicated` (D),
`fragmented` (F), `missing` (M), `completenessPercent` (C/n·100), `duplicationPercent` (D/n·100),
`singleCopyPercent`, `fragmentedPercent`, `missingPercent`, `buscoNotation`.

## Errors

| Code | Message |
|------|---------|
| 1001 | Assembly cannot be null |
| 1002 | Marker genes cannot be null |
| 1003 | Identity threshold must be in [0, 1] |
| 1004 | Coverage threshold must be in [0, 1] |
| 1005 | Minimum length fraction must be in [0, 1] |

An empty or null marker sequence, or a null assembly sequence, is rejected by the library (ArgumentException).

## Example

Assembly `s1 = s2 = A40C40`; markers `m1 = A40C40`, `m2 = (GT)40` →
`totalGenes 2, completeDuplicated 1, missing 1, completenessPercent 50, duplicationPercent 50,
buscoNotation "C:50.0%[S:0.0%,D:50.0%],F:0.0%,M:50.0%,n:2"`.

## References

- [GenomeAssemblyAnalyzer.AssessCompleteness](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L727)
- Simão FA et al. (2015) BUSCO. Bioinformatics 31:3210; Manni M et al. (2021) BUSCO v5. Mol Biol Evol 38:4647.
- Huang N, Li H (2023) compleasm. Bioinformatics 39:btad595.
