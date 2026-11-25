# find_repetitive_elements

Find tandem and inverted repeats in a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_repetitive_elements` |
| **Method ID** | `GenomeAnnotator.FindRepetitiveElements` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans a DNA sequence for two classes of exact repeats, delegating to the canonical `RepeatFinder`:

- **Tandem repeats** — `RepeatFinder.FindMicrosatellites(seq, 1, 60, minCopies)`: every maximal perfect run of
  period 1–60 is reported once, at its left end, by its primitive A/C/G/T unit and its complete copies
  (`end − start = copies × unit`, copies ≥ `minCopies`, span ≥ `minRepeatLength`). Rotations of a run, non-primitive
  units and runs of `N` are not reported (Kolpakov & Kucherov 1999; MISA / pytrf convention).
- **Inverted repeats** (form W·G·W̄ᴿ) — `RepeatFinder.FindInvertedRepeats(seq, max(minRepeatLength, 2), 100,
  minLoopLength: 0, maxArmLength: max(100, arm minimum))`: the maximal exact stems, i.e. the result set of EMBOSS
  `palindrome -minpallen m -maxpallen 100 -gaplimit 100 -nummismatches 0 -overlap Y` (EMBOSS / IUPACpal default
  limits). Only A·T / C·G pair; sub-stems of a longer stem are not reported. For inverted repeats `minRepeatLength`
  is the minimum **arm** length. The `sequence` is `leftArm + rightArm` (loop 0) or `leftArm[loop]rightArm`.

Each hit reports a 0-based inclusive `start`, an exclusive `end`, a `type` of `tandem_repeat` or
`inverted_repeat`, and the repeat `sequence`. Tandem repeats come first (by period, then start), then inverted
repeats (by left-arm start, then right-arm start).

## Core Documentation Reference

- Source: [GenomeAnnotator.cs#L2059](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs#L2059)

## Input Schema

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `dnaSequence` | string | Yes | — | DNA sequence to scan |
| `minRepeatLength` | integer | No | 10 | Minimum tandem-repeat span and minimum inverted-repeat arm length (nt) |
| `minCopies` | integer | No | 2 | Minimum number of tandem copies (≥2) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `repeats` | array | `{ start, end, type, sequence }` per repeat element |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Tandem repeat

`ATTCGATTCGATTCG` (ATTCG × 3), `minRepeatLength=5` → one `tandem_repeat` spanning [0, 15).

**Response:**
```json
{ "repeats": [ { "start": 0, "end": 15, "type": "tandem_repeat", "sequence": "ATTCGATTCGATTCG" } ] }
```

### Example 2: Inverted repeat

`GAATTC` (the EcoRI palindrome) → one `inverted_repeat` spanning [0, 6).

### Example 3: Maximal stem only

`AGGATCCT`, `minRepeatLength=3` → one `inverted_repeat` `AGGATCCT` [0, 8) (EMBOSS palindrome); its sub-stems
`AGG[AT]CCT` and `GGATCC` are not reported.

## Performance

- **Time Complexity:** O(n · 60) tandem scan; O(n · 100 · arm) inverted-repeat scan (loop ≤ 100)
- **Space Complexity:** O(k)

## See Also

- [find_promoter_motifs](find_promoter_motifs.md) - Promoter motifs
- [codon_usage](codon_usage.md) - Codon usage counts
