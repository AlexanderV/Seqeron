# create_mirna

Build a MiRna record (T→U normalized) with seed metadata.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `create_mirna` |
| **Method ID** | `MiRnaAnalyzer.CreateMiRna` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Constructs a `MiRna` record from a name and nucleotide sequence. The sequence is upper-cased and
T→U normalised (so DNA input becomes RNA), and the positions 2–8 seed is extracted with
`seedStart = 1`, `seedEnd = 7`.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L107](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L107)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `name` | string | Yes | miRNA name / identifier |
| `sequence` | string | Yes | miRNA nucleotide sequence |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `miRna.name` | string | miRNA name |
| `miRna.sequence` | string | Upper-cased, T→U normalised sequence |
| `miRna.seedSequence` | string | Seed region (positions 2–8) |
| `miRna.seedStart` | integer | Always 1 |
| `miRna.seedEnd` | integer | Always 7 |

## Errors

| Code | Message |
|------|---------|
| 1001 | Name cannot be null or empty |
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: DNA input normalised to RNA

**Response:**
```json
{ "miRna": { "name": "let-7a", "sequence": "UGAGGUAGUAGGUUGUAUAGUU", "seedSequence": "GAGGUAG", "seedStart": 1, "seedEnd": 7 } }
```

### Example 2: RNA input preserved

**Response:**
```json
{ "miRna": { "name": "miR-21", "sequence": "UAGCUUAUCAGACUGAUGUUGA", "seedSequence": "AGCUUAU", "seedStart": 1, "seedEnd": 7 } }
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(n)

## See Also

- [mirna_seed_sequence](mirna_seed_sequence.md) — seed extraction only
- [find_mirna_target_sites](find_mirna_target_sites.md) — target-site search using the record
