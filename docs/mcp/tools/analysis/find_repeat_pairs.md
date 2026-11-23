# find_repeat_pairs

Maximal exact repeated pairs (MUMmer repeat-match list) on the forward strand and / or as reverse complements.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_repeat_pairs` |
| **Method ID** | `GenomicAnalyzer.FindRepeatPairs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Every maximal exact repeated pair of the sequence (Gusfield 1997 §7.12; MUMmer `repeat-match`, REPuter): two copies that cannot be extended on either side, reported as the 0-based forward-strand starts of both copies, the length and the spacing. The list is the union of the canonical direct and reverse-complement repeat finders with every length and spacing filter lifted, so it is the unfiltered `repeat-match` output (verified against MUMmer 3.23). Use `find_direct_repeats` / `find_reverse_complement_repeats` for the spacing and length filters and `find_repeats` for the repeated STRINGS with all their positions.

## Core Documentation Reference

- Source: [GenomicAnalyzer.Repeats.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.Repeats.cs)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence: A/C/G/T plus IUPAC codes such as N (case-insensitive); non-A/C/G/T symbols never repeat. |
| `minLength` | integer | No | Minimum repeat length (repeat-match -n; default 5, >= 2). (default: 5) |
| `strands` | string | No | Which strands (case-insensitive): Both (default) \| Forward \| ReverseComplement. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | Pairs: `{ firstPosition, secondPosition, length, spacing, isReverseComplement }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty / invalid nucleotide sequence |
| 1002 | minLength must be at least 2 / unknown strands |

## Examples

### Example 1: Direct pair

**Input:** `{ "sequence": "ACGTTTACGAAAT", "minLength": 3, "strands": "Forward" }`

**Response:** `{ "items": [ { "firstPosition": 0, "secondPosition": 6, "length": 3, "spacing": 3 } ] }` (`repeat-match -n 3` prints `1 7 3`, 1-based).

## See Also

- [find_repeats](find_repeats.md)
- [find_direct_repeats](find_direct_repeats.md)
- [find_reverse_complement_repeats](find_reverse_complement_repeats.md)
