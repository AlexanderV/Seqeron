# find_microsatellites

Short Tandem Repeats (STRs / microsatellites) in a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_microsatellites` |
| **Method ID** | `RepeatFinder.FindMicrosatellites` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds **microsatellites** (Short Tandem Repeats): units of length
`minUnitLength..maxUnitLength` repeated at least `minRepeats` times consecutively.
Each hit is classified by repeat type (Mononucleotide, Dinucleotide, …). Redundant
units — those that are themselves a shorter unit repeated (e.g. `AA`) — are skipped.
Per unit length each maximal perfect run is reported once, at its left end, with the
number of complete copies (a trailing partial copy is not counted, and rotations of
the same run — e.g. `TA` inside `ATATATA` — are not re-reported); units containing
non-ACGT symbols (e.g. `N`) are never reported. Positions are 0-based.

## Core Documentation Reference

- Source: [RepeatFinder.cs#L85](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs#L85)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `minUnitLength` | integer | No | Minimum unit length (default 1, ≥ 1) |
| `maxUnitLength` | integer | No | Maximum unit length (default 6) |
| `minRepeats` | integer | No | Minimum repeats (default 3, ≥ 2) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | `{ position, repeatUnit, repeatCount, totalLength, repeatType }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Invalid unit-length or minRepeats bounds |

## Examples

### Example 1: Dinucleotide (CA)₄

**User Prompt:**
> Find microsatellites in "CACACACA" (unit 2–6, ≥ 3 repeats).

**Expected Tool Call:**
```json
{
  "tool": "find_microsatellites",
  "arguments": { "sequence": "CACACACA", "minUnitLength": 2, "maxUnitLength": 6, "minRepeats": 3 }
}
```

**Response:**
```json
{ "items": [ { "position": 0, "repeatUnit": "CA", "repeatCount": 4, "totalLength": 8, "repeatType": "Dinucleotide" } ] }
```

### Example 2: Trinucleotide (CAG)₃

**User Prompt:**
> Trinucleotide STRs in "CAGCAGCAG".

**Expected Tool Call:**
```json
{
  "tool": "find_microsatellites",
  "arguments": { "sequence": "CAGCAGCAG", "minUnitLength": 3, "maxUnitLength": 6, "minRepeats": 3 }
}
```

**Response:**
```json
{ "items": [ { "position": 0, "repeatUnit": "CAG", "repeatCount": 3, "totalLength": 9, "repeatType": "Trinucleotide" } ] }
```

## Performance

- **Time Complexity:** O(n · (maxUnitLength − minUnitLength + 1)) character comparisons.
- **Space Complexity:** O(number of STRs).

## See Also

- [find_tandem_repeats](find_tandem_repeats.md) — general tandem repeats
- [tandem_repeat_summary](tandem_repeat_summary.md) — aggregate STR statistics
