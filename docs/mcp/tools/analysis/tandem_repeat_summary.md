# tandem_repeat_summary

Aggregate statistics across all microsatellites in a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `tandem_repeat_summary` |
| **Method ID** | `RepeatFinder.GetTandemRepeatSummary` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Summarizes all microsatellites (STRs, unit length 1–6) found in a DNA sequence:
the total number of repeats, total repeat bases, the percentage of the sequence
covered by repeats (union of spans, so ≤ 100), per-type counts (mono/di/tri/tetra),
the longest repeat and the most frequent repeat unit.

## Core Documentation Reference

- Source: [RepeatFinder.cs#L871](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs#L871)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `minRepeats` | integer | No | Minimum repeats (default 3, ≥ 2) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `totalRepeats` | integer | Number of microsatellites found |
| `totalRepeatBases` | integer | Sum of repeat spans (may overlap) |
| `percentageOfSequence` | number | Percent of the sequence covered (0–100) |
| `mononucleotideRepeats` | integer | Count of mononucleotide STRs |
| `dinucleotideRepeats` | integer | Count of dinucleotide STRs |
| `trinucleotideRepeats` | integer | Count of trinucleotide STRs |
| `tetranucleotideRepeats` | integer | Count of tetranucleotide STRs |
| `longestRepeat` | object/null | The longest microsatellite (or null) |
| `mostFrequentUnit` | string/null | The most frequent repeat unit (or null) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |

## Examples

### Example 1: Single trinucleotide STR

**User Prompt:**
> Summarize tandem repeats in "CAGCAGCAG".

**Expected Tool Call:**
```json
{
  "tool": "tandem_repeat_summary",
  "arguments": { "sequence": "CAGCAGCAG", "minRepeats": 3 }
}
```

**Response:**
```json
{ "totalRepeats": 1, "totalRepeatBases": 9, "percentageOfSequence": 100.0, "mononucleotideRepeats": 0, "dinucleotideRepeats": 0, "trinucleotideRepeats": 1, "tetranucleotideRepeats": 0, "mostFrequentUnit": "CAG" }
```

### Example 2: No STRs

**User Prompt:**
> Tandem repeat summary of "ACGT".

**Expected Tool Call:**
```json
{
  "tool": "tandem_repeat_summary",
  "arguments": { "sequence": "ACGT", "minRepeats": 3 }
}
```

**Response:**
```json
{ "totalRepeats": 0, "totalRepeatBases": 0, "percentageOfSequence": 0.0, "mononucleotideRepeats": 0, "dinucleotideRepeats": 0, "trinucleotideRepeats": 0, "tetranucleotideRepeats": 0 }
```

## Performance

- **Time Complexity:** O(n · 6) STR scan + O(k log k) grouping.
- **Space Complexity:** O(number of STRs).

## See Also

- [find_microsatellites](find_microsatellites.md) — the per-STR list
- [find_tandem_repeats](find_tandem_repeats.md) — general tandem repeats
