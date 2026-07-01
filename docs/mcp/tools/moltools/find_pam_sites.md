# find_pam_sites

Find all CRISPR PAM sites (both strands) in a sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `find_pam_sites` |
| **Method ID** | `CrisprDesigner.FindPamSites` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans both strands for the PAM of the chosen CRISPR system (IUPAC-aware, e.g. `NGG`, `NNGRRT`, `TTTV`) and returns each match with its forward-strand position, matched PAM, adjacent guide/target window, target start, and strand. Sites whose guide window would fall outside the sequence are skipped.

## Core Documentation Reference

- Source: [CrisprDesigner.cs#L51](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs#L51)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to scan (non-empty). |
| `system_type` | enum | No | CRISPR system (default `SpCas9`). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sites` | array | PAM sites, each `{position, pamSequence, targetSequence, targetStart, isForwardStrand}`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: SpCas9

`ATATATATATATATATATATTGG` → forward PAM `TGG` at position 20; its 20-nt target is the preceding `ATAT…AT`.

### Example 2: `TGG` (too short) → `sites: []`.

## See Also

- [crispr_system_info](crispr_system_info.md), [design_guide_rnas](design_guide_rnas.md)
