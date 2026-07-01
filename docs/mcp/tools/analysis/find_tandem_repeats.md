# find_tandem_repeats

Consecutive repeating units (tandem repeats) in a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_tandem_repeats` |
| **Method ID** | `GenomicAnalyzer.FindTandemRepeats` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds **tandem repeats** — runs of a unit repeated consecutively (e.g. `ATGATGATG`) —
whose unit length is ≥ `minUnitLength` and which repeat ≥ `minRepetitions` times. Each
result reports the repeat unit, its 0-based start position, the number of consecutive
copies, and the total span (`unit.Length × repetitions`). The scan skips past each
tandem it reports.

## Core Documentation Reference

- Source: [GenomicAnalyzer.cs#L115](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs#L115)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `minUnitLength` | integer | No | Minimum repeat-unit length (default 2, ≥ 1) |
| `minRepetitions` | integer | No | Minimum consecutive copies (default 2, ≥ 2) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | `{ unit, position, repetitions, totalLength }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |
| 1002 | Minimum repeat-unit length must be at least 1 |
| 1002 | minRepetitions must be ≥ 2 |

## Examples

### Example 1: Trinucleotide tandem (ATG ×3)

**User Prompt:**
> Find tandem 3-mers in "ATGATGATG".

**Expected Tool Call:**
```json
{
  "tool": "find_tandem_repeats",
  "arguments": { "sequence": "ATGATGATG", "minUnitLength": 3, "minRepetitions": 2 }
}
```

**Response:**
```json
{ "items": [ { "unit": "ATG", "position": 0, "repetitions": 3, "totalLength": 9 } ] }
```

### Example 2: Mononucleotide run (unit length 1)

**User Prompt:**
> Find tandem repeats in "AAAAA" with unit length ≥ 1.

**Expected Tool Call:**
```json
{
  "tool": "find_tandem_repeats",
  "arguments": { "sequence": "AAAAA", "minUnitLength": 1, "minRepetitions": 2 }
}
```

**Response:**
```json
{ "items": [ { "unit": "A", "position": 0, "repetitions": 5, "totalLength": 5 }, { "unit": "AA", "position": 0, "repetitions": 2, "totalLength": 4 } ] }
```
With `minUnitLength = 1` the scan also considers unit length 2, so the run is
reported both as A×5 and as the AA×2 tiling.

## Performance

- **Time Complexity:** O(n²) over unit lengths and start positions.
- **Space Complexity:** O(number of tandems).

## See Also

- [find_microsatellites](find_microsatellites.md) — STRs with repeat-type classification
- [find_repeats](find_repeats.md) — any repeated substring
- [tandem_repeat_summary](tandem_repeat_summary.md) — aggregate statistics
