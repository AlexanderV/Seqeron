# standardize_repeat_motif

MISA repeat-type class and Krait standard motif of a repeat unit.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `standardize_repeat_motif` |
| **Method ID** | `RepeatFinder.GetStandardMotif` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Standardizes a repeat unit two ways. `canonicalClass` — the MISA repeat-type class "considering sequence
complementary" (`misa.pl` `.statistics`, `RepeatFinder.GetCanonicalMotifClass`): `X/Y` with X, Y the
lexicographically smallest rotations of the motif and of its reverse complement, smaller first (CA → `AC/GT`).
`standardMotif` — the Krait standard motif (Du et al. 2018, lmdu/krait `motif.py` `StandardMotif.standard`): the
smallest member of the motif's equivalence set under Krait's base order A < T < C < G at `level` 0 (motif itself),
1 (+ rotations), 2 (+ reverse-complement rotations; default, the MISA grouping), 3 (+ complement; Krait GUI
default) or 4 (+ reverse). Motif must be A/C/G/T only (case-insensitive). `tandem_repeat_summary` reports the
per-class / per-standard-motif SSR counts.

## Core Documentation Reference

- Source: [RepeatFinder.cs#L740](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs#L740)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `motif` | string | Yes | Repeat unit of A/C/G/T (case-insensitive) |
| `level` | integer | No | Krait standardization level (default 2) |

## Output Schema

`canonicalClass` (MISA class), `standardMotif` (Krait), `level`

## Errors

| Code | Message |
|------|---------|
| 1001 | Motif cannot be null or empty |
| 1002 | Motif contains a non-ACGT symbol |
| 1003 | level must be 0-4 |

## Examples

### Example 1: Dinucleotide CA

**Input:** `{"motif": "CA"}`

**Output:**

```json
{"canonicalClass": "AC/GT", "standardMotif": "AC", "level": 2}
```

### Example 2: Krait GUI level 3

**Input:** `{"motif": "CTG", "level": 3}`

**Output:**

```json
{"canonicalClass": "AGC/CTG", "standardMotif": "ACG", "level": 3}
```

### Example 3: Level 0 = motif itself

**Input:** `{"motif": "ACAT", "level": 0}`

**Output:**

```json
{"canonicalClass": "ACAT/ATGT", "standardMotif": "ACAT", "level": 0}
```

## Performance

- O(m²) for a motif of length m.

## See Also

- [tandem_repeat_summary](tandem_repeat_summary.md)
- [find_microsatellites](find_microsatellites.md)
